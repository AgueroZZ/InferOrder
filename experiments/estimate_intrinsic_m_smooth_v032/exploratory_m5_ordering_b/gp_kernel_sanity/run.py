"""GP covariance replacements and matched-prior CAVI on the saved ordering-B data."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[key]='1'
import sys,json,time,hashlib,warnings
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.special import logsumexp,xlogy
from scipy.optimize import minimize
from scipy.stats import spearmanr,norm
import GPy
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
sys.path.insert(0,str(BASE/'collapsed_comparison'))
import compare
JITTER=1e-6

def kernel(grid,kind,ell=1.,amp=1.):
    dist=np.abs(grid[:,None]-grid[None,:])/ell
    if kind=='rbf': C=np.exp(-.5*dist**2)
    elif kind=='matern32': C=(1+np.sqrt(3)*dist)*np.exp(-np.sqrt(3)*dist)
    else: raise ValueError(kind)
    return amp*(C+JITTER*np.eye(len(grid)))

def sparse_covariance(C,neighbors):
    """Ordered nearest-predecessor conditionals give an SPD banded precision."""
    K=len(C);B=np.eye(K);v=np.zeros(K)
    for i in range(K):
        idx=np.arange(max(0,i-neighbors),i)
        if len(idx):
            beta=np.linalg.solve(C[np.ix_(idx,idx)],C[idx,i]);B[i,idx]=-beta
            v[i]=C[i,i]-C[i,idx]@beta
        else:v[i]=C[i,i]
    assert np.min(v)>0
    Q=(B.T/v)@B
    return cho_solve(cho_factor(Q,lower=True),np.eye(K)),Q

def corr(a,b):return float(abs(spearmanr(a,b).statistic))
def save(tag,result,arrays):
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT/(tag+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(ROOT/(tag+'.npz'),**arrays)
    print(tag,{k:result.get(k) for k in ('rho','rho_vs_bgplvm','iterations','seconds','converged','elbo')},flush=True)

def swap(case,kind,neighbors=0):
    ref=json.loads((BASE/'collapsed_comparison'/f'{case}_reference.json').read_text())
    K=len(ref['Q']);grid=np.linspace(-np.sqrt(3),np.sqrt(3),K)
    C=kernel(grid,kind);exact=C.copy()
    if neighbors:C,Q=sparse_covariance(C,neighbors)
    else:Q=cho_solve(cho_factor(C,lower=True),np.eye(K))
    ref['Q']=((Q+Q.T)/2).tolist();ref['rank']=K;ref['logdet_Q']=-np.linalg.slogdet(C)[1]
    ref['initial']['lambda']=[1.]*len(ref['initial']['lambda'])
    model=compare.Model(ref);R=np.array(ref['initial']['R']);start=time.perf_counter()
    fit=compare.cavi(model,R,maxiter=4000)
    final=fit['R']@grid
    result=dict(case=case,method='Q_swap',kernel=kind,neighbors=neighbors,rho=corr(ref['truth'],final),
       iterations=fit['iterations'],seconds=time.perf_counter()-start,converged=fit['converged'],elbo=fit['value'],
       covariance_relative_error=float(np.linalg.norm(C-exact)/np.linalg.norm(exact)),
       precision_density=float(np.mean(abs(Q)>1e-10)),min_eigenvalue=float(np.linalg.eigvalsh(C).min()),
       initialization='Native MPCurve hard PCA assignments; kernel amplitude 1; native noise and pi',
       grid=grid.tolist(),history=fit['history'])
    save(f'{case}_swap_{kind}_m{neighbors}',result,dict(final=final,initial=R@grid,R=fit['R'],sigma=fit['sigma'],lam=fit['lam'],pi=fit['pi'],covariance=C,precision=Q))

class GridGP:
    """Mean-field categorical positions and Gaussian GP values, with a proper prior."""
    def __init__(self,Y,grid,kind):
        self.Y=Y;self.grid=grid;self.kind=kind;self.n,self.d=Y.shape;self.K=len(grid)
        self.pi=np.exp(-.5*grid**2);self.pi/=self.pi.sum();self.y2=np.sum(Y*Y)
    def evaluate(self,R,theta):
        amp,ell,noise=np.exp(theta);C=kernel(self.grid,self.kind,ell,amp)
        # Covariance form avoids forming a poorly conditioned RBF precision.
        L=np.linalg.cholesky(C);counts=R.sum(0);b=R.T@self.Y/noise
        H=np.eye(self.K)+(L.T*(counts/noise))@L
        ch=cho_factor(H,lower=True);S=L@cho_solve(ch,L.T);mu=S@b
        logdetH=2*np.log(np.diag(ch[0])).sum()
        value=-.5*(self.n*self.d*np.log(2*np.pi*noise)+self.y2/noise)+.5*np.sum(b*mu)-.5*self.d*logdetH
        value+=np.sum(R.sum(0)*np.log(self.pi))-np.sum(xlogy(R,R))
        return float(value),mu,S,C
    def fit(self,R,theta,maxiter=2000):
        start=time.perf_counter();value,mu,S,C=self.evaluate(R,theta);history=[];opt_calls=0
        for it in range(1,maxiter+1):
            noise=np.exp(theta[2]);score=self.Y@mu.T/noise-.5*(np.sum(mu*mu,axis=1)+self.d*np.diag(S))/noise+np.log(self.pi)
            R=np.exp(score-logsumexp(score,axis=1,keepdims=True))
            # Conditional maximization of just three shared hyperparameters.
            if it==1 or it%5==0:
                opt=minimize(lambda t:-self.evaluate(R,t)[0]/(self.n*self.d),theta,method='L-BFGS-B',
                    bounds=[(-8,5),(np.log(.03),np.log(10)),(-12,4)],options=dict(maxiter=12,ftol=1e-10,gtol=1e-7))
                opt_calls+=opt.nfev
                if -opt.fun*self.n*self.d>=self.evaluate(R,theta)[0]-1e-7:theta=opt.x
            new,mu,S,C=self.evaluate(R,theta);history.append(dict(iteration=it,elbo=new,seconds=time.perf_counter()-start))
            delta=new-value;value=new
            if it%5==0 and delta>=-1e-7 and delta/(self.n*self.d)<1e-6:break
        return dict(R=R,theta=theta,mu=mu,S=S,history=history,iterations=it,converged=it<maxiter,
                    seconds=time.perf_counter()-start,elbo=value,hyper_evaluations=opt_calls)

def aligned(case,K=100,kind='rbf'):
    raw=np.loadtxt(BASE/'external_methods/inputs'/f'{case}_X.csv',delimiter=',',skiprows=1)
    Y=raw-raw.mean(0);positions=np.genfromtxt(BASE/'external_methods/inputs'/f'{case}_positions.csv',delimiter=',',names=True)
    np.random.seed(20260929)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore');gpy=GPy.models.BayesianGPLVM(Y,1,num_inducing=50)
    initial=gpy.X.mean.values[:,0].copy();variance=gpy.X.variance.values[:,0].copy()
    grid=np.linspace(-4,4,K);edges=np.r_[-np.inf,(grid[1:]+grid[:-1])/2,np.inf]
    R=np.diff(norm.cdf((edges[None,:]-initial[:,None])/np.sqrt(variance[:,None])),axis=1)
    R=np.maximum(R,1e-300);R/=R.sum(1,keepdims=True)
    theta=np.log([float(gpy.kern.variance[0]),float(gpy.kern.lengthscale[0]),float(gpy.likelihood.variance[0])])
    gp=GridGP(Y,grid,kind);fit=gp.fit(R,theta)
    bgpath=next((BASE/'external_methods/gpy_results').glob(f'*_{case}_BayesianGPLVM_pca_m50_s20260929.json'))
    bg=json.loads(bgpath.read_text());final=fit['R']@grid
    save(f'{case}_aligned_{kind}_K{K}',dict(case=case,method='aligned_grid_gp',kernel=kind,K=K,
      rho=corr(positions['truth'],final),rho_vs_bgplvm=corr(bg['final'],final),initial_rho=corr(positions['truth'],initial),
      amplitude=float(np.exp(fit['theta'][0])),lengthscale=float(np.exp(fit['theta'][1])),noise=float(np.exp(fit['theta'][2])),
      edge_mass=float(fit['R'][:,[0,-1]].sum()/len(R)),bgplvm_reference=str(bgpath.relative_to(BASE)),
      **{k:fit[k] for k in ['history','iterations','converged','seconds','elbo','hyper_evaluations']}),
      dict(R=fit['R'],final=final,initial=initial,grid=grid,mean=fit['mu'],covariance=fit['S'],theta=fit['theta']))

def verify():
    grid=np.linspace(-2,2,15);rng=np.random.default_rng(17);Y=rng.normal(size=(20,3));idx=rng.integers(0,15,20)
    R=np.eye(15)[idx];gp=GridGP(Y,grid,'rbf');theta=np.log([1.3,.7,.2]);val,mu,S,C=gp.evaluate(R,theta)
    Ky=C[np.ix_(idx,idx)]+.2*np.eye(len(idx));ch=cho_factor(Ky,lower=True)
    exact_mu=C[:,idx]@cho_solve(ch,Y);exact_S=C-C[:,idx]@cho_solve(ch,C[idx,:])
    loglik=-.5*(len(Y)*3*np.log(2*np.pi)+3*2*np.log(np.diag(ch[0])).sum()+np.sum(Y*cho_solve(ch,Y)))
    checks=dict(mean_error=float(np.max(abs(mu-exact_mu))),covariance_error=float(np.max(abs(S-exact_S))),
      evidence_error=abs(val-loglik-np.log(gp.pi[idx]).sum()))
    for kind in ['rbf','matern32']:
        k=GPy.kern.RBF(1,variance=1.3,lengthscale=.7) if kind=='rbf' else GPy.kern.Matern32(1,variance=1.3,lengthscale=.7)
        checks[kind+'_kernel_error']=float(np.max(abs(kernel(grid,kind,.7,1.3)-k.K(grid[:,None])-1.3*JITTER*np.eye(len(grid)))))
    assert max(checks.values())<1e-8,checks
    (ROOT/'verification.json').write_text(json.dumps(checks,indent=2));print('VERIFIED',checks,flush=True)

if __name__=='__main__':
    verify()
    for case in ('B_noise0','B_noise05','B_noise1'):
        for kind in ('rbf','matern32'):
            for neighbors in (0,5):swap(case,kind,neighbors)
        aligned(case)
    aligned('B_noise05',200)
    aligned('B_noise05',100,'matern32')
