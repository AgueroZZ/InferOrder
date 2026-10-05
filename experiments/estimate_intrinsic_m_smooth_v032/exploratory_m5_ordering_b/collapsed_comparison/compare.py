"""Profile the Gaussian curve block of the frozen single-ordering MPCurve ELBO."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import argparse, hashlib, json, pathlib, time
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from scipy.special import logsumexp, xlogy
from scipy.stats import spearmanr
ROOT=pathlib.Path(__file__).resolve().parent

def rho(truth,R):
    return float(abs(spearmanr(truth,R@np.linspace(0,1,R.shape[1])).statistic))

class Model:
    def __init__(self,ref):
        self.encoding='logits';self.ref=ref;self.X=np.array(ref['X']);self.Q=np.array(ref['Q'])
        self.n,self.d=self.X.shape;self.K=self.Q.shape[0]
        self.rank=ref['rank'];self.logdetQ=ref['logdet_Q']
        self.sumx2=(self.X**2).sum(0)
        self.eigenvalues,self.V=np.linalg.eigh(self.Q)
        self.eigenvalues[:self.K-self.rank]=0.
    def posterior(self,R,sigma,lam,clip=False):
        Nk=R.sum(0);counts=np.maximum(Nk,1e-8) if clip else Nk
        GX=self.X.T@R;means=[];covs=[];logdets=[];rough=[]
        # Work in RW2 eigen-coordinates: exact null modes and stable large-lambda algebra.
        for j in range(self.d):
            A=(self.V.T*(counts/sigma[j]))@self.V+np.diag(lam[j]*self.eigenvalues)
            chol=cho_factor(A,lower=True,check_finite=False)
            mc=cho_solve(chol,self.V.T@(GX[j]/sigma[j]),check_finite=False)
            Sc=cho_solve(chol,np.eye(self.K),check_finite=False)
            means.append(self.V@mc);covs.append(self.V@Sc@self.V.T)
            rough.append(np.sum(self.eigenvalues*(mc*mc+np.diag(Sc))))
            logdets.append(2*np.log(np.diag(chol[0])).sum())
        return np.array(means),np.array(covs),np.array(logdets),np.array(rough)
    def quantities(self,R,sigma,lam,pi,post=None):
        m,S,logdetA,roughness=self.posterior(R,sigma,lam) if post is None else post
        variance=np.diagonal(S,axis1=1,axis2=2);Nk=R.sum(0);GX=self.X.T@R
        residual=self.sumx2-2*np.sum(m*GX,axis=1)+np.sum((m*m+variance)*Nk,axis=1)
        likelihood=-.5*(self.n*self.d*np.log(2*np.pi)+self.n*np.log(sigma).sum()+(residual/sigma).sum())
        prior=.5*np.sum(self.rank*(np.log(lam)-np.log(2*np.pi))+self.logdetQ-lam*roughness)
        entropyU=.5*np.sum(self.K*np.log(2*np.pi*np.e)-logdetA)
        categorical=np.sum(Nk*np.log(np.maximum(pi,1e-300)))-np.sum(xlogy(R,R))
        score=self.X@(m/sigma[:,None])-.5*np.sum((m*m+variance)/sigma[:,None],axis=0)
        return likelihood+prior+entropyU+categorical,score,residual,roughness,(m,S,logdetA,roughness)
    def profile_value_gradient(self,z,fixed=False):
        logits=z[:self.n*self.K].reshape(self.n,self.K)
        if self.encoding=='sqrt':
            denom=np.sum(logits*logits,axis=1,keepdims=True)
            R=logits*logits/denom;logR=np.log(np.maximum(R,1e-300))
        else:
            logR=logits-logsumexp(logits,axis=1,keepdims=True);R=np.exp(logR)
        init=self.ref['initial']
        sigma=np.array(init['sigma2']) if fixed else np.exp(z[-self.d:])
        lam=np.array(init['lambda']) if fixed else np.exp(z[-2*self.d:-self.d])
        pi=np.array(init['pi']) if fixed else R.mean(0)
        value,score,resid,rough,post=self.quantities(R,sigma,lam,pi)
        gradR=score+np.log(np.maximum(pi,1e-300))-logR
        centered=gradR-np.sum(R*gradR,axis=1,keepdims=True)
        gradlogits=2*logits/denom*centered if self.encoding=='sqrt' else R*centered
        grad=gradlogits.ravel()
        if not fixed:
            grad=np.r_[grad,.5*(self.rank-lam*rough),-.5*self.n+.5*resid/sigma]
        # Independent determinant expression for the same profile objective.
        m,S,logdetA,_=post; b=(self.X.T@R)/sigma[:,None]
        profile=(-.5*(self.n*self.d*np.log(2*np.pi)+self.n*np.log(sigma).sum()+(self.sumx2/sigma).sum())
                 +.5*np.sum(self.rank*np.log(lam)+self.logdetQ+(self.K-self.rank)*np.log(2*np.pi)
                             +np.sum(b*m,axis=1)-logdetA)
                 +np.sum(R.sum(0)*np.log(np.maximum(pi,1e-300)))-np.sum(xlogy(R,R)))
        return value,grad,dict(R=R,sigma=sigma,lam=lam,pi=pi,post=post,profile_difference=float(profile-value))

def cavi(model,R0,fixed=False,maxiter=2000,tol=1e-6,trace=True,state=None):
    init=model.ref['initial'];R=R0.copy();sigma=np.array(init['sigma2']);lam=np.array(init['lambda']);pi=np.array(init['pi'])
    if state is not None:
        sigma=state['sigma'].copy();lam=state['lam'].copy();pi=state['pi'].copy()
    start=time.perf_counter();post=model.posterior(R,sigma,lam,clip=True)
    value=model.quantities(R,sigma,lam,pi,post)[0]
    history=[dict(iteration=0,seconds=0.,elbo=float(value))];converged=False;states=[]
    for it in range(1,maxiter+1):
        score=model.quantities(R,sigma,lam,pi,post)[1]+np.log(np.maximum(pi,1e-300))
        R=np.exp(score-logsumexp(score,axis=1,keepdims=True))
        if not fixed:
            pi=np.maximum(R.mean(0),np.finfo(float).eps);pi/=pi.sum()
            _,_,resid,rough,_=model.quantities(R,sigma,lam,pi,post)
            sigma=np.clip(resid/model.n,1e-10,1e10)
            lam=np.clip(model.rank/np.maximum(rough,1e-12),1e-10,1e10)
        post=model.posterior(R,sigma,lam,clip=True)
        new=model.quantities(R,sigma,lam,pi,post)[0]
        if it<=3:states.append(dict(R=R.copy(),sigma=sigma.copy(),lam=lam.copy(),value=float(new),post=post))
        if trace:history.append(dict(iteration=it,seconds=time.perf_counter()-start,elbo=float(new)))
        delta=new-value;value=new
        if delta>=0 and delta/(model.n*model.d)<tol:
            converged=True;break
    return dict(R=R,sigma=sigma,lam=lam,pi=pi,value=float(value),iterations=it,converged=converged,
                elapsed=time.perf_counter()-start,history=history,states=states)

def validate(model):
    init=model.ref['initial'];R=np.array(init['R']);sigma=np.array(init['sigma2']);lam=np.array(init['lambda']);pi=np.array(init['pi'])
    value=model.quantities(R,sigma,lam,pi)[0]
    errs={'initial_elbo':abs(value-init['elbo'])}
    replay=cavi(model,R,maxiter=3,tol=0)
    for i,(got,expected) in enumerate(zip(replay['states'],model.ref['stages'])):
        errs[f'R_step{i+1}']=float(np.max(abs(got['R']-expected['R'])))
        errs[f'elbo_step{i+1}']=abs(got['value']-expected['elbo'])
    rng=np.random.default_rng(123);soft=.99*R+.01/model.K
    z=np.r_[(np.sqrt(soft) if model.encoding=='sqrt' else np.log(soft)).ravel(),np.log(lam),np.log(sigma)]
    value,grad,state=model.profile_value_gradient(z)
    errs['profile_expression']=abs(state['profile_difference'])
    # Directional derivatives test the large categorical block and both hyperparameter blocks separately.
    for group,slc in [('positions',slice(0,model.n*model.K)),('lambda',slice(-2*model.d,-model.d)),('sigma',slice(-model.d,None))]:
        v=np.zeros_like(z);v[slc]=rng.normal(size=v[slc].shape);v/=np.linalg.norm(v)
        h=1e-4;fd=(model.profile_value_gradient(z+h*v)[0]-model.profile_value_gradient(z-h*v)[0])/(2*h)
        errs[f'gradient_{group}']=float(abs(fd-grad@v)/max(1,abs(fd)))
    if max(errs.values())>1e-5:raise AssertionError(errs)
    return errs

def run(case,epsilon=.01,fixed=False,repeats=3):
    ref=json.loads((ROOT/f'{case}_reference.json').read_text());model=Model(ref)
    validation=validate(model);print('VALIDATED',case,validation,flush=True)
    R0=(1-epsilon)*np.array(ref['initial']['R'])+epsilon/model.K
    modes='fixed' if fixed else 'adaptive';tag=f'{case}_{modes}_eps{epsilon:g}'
    outputs=[]
    for method in ('cavi_python','collapsed_lbfgs','collapsed_sqrt'):
        model.encoding='sqrt' if method=='collapsed_sqrt' else 'logits'
        validation=validate(model)
        timings=[];result=None
        for repeat in range(repeats):
            if method=='cavi_python':
                result=cavi(model,R0,fixed=fixed)
                timings.append(result['elapsed'])
            else:
                z0=(np.sqrt(R0) if model.encoding=='sqrt' else np.log(R0)).ravel()
                if not fixed:z0=np.r_[z0,np.log(ref['initial']['lambda']),np.log(ref['initial']['sigma2'])]
                start=time.perf_counter();history=[];evaluations=0;last={};invalid_trials=0
                def objective(z):
                    nonlocal evaluations,last,invalid_trials
                    evaluations+=1
                    try:
                        value,grad,state=model.profile_value_gradient(z,fixed=fixed)
                    except np.linalg.LinAlgError:
                        # Reject numerically singular trial points without changing the prior.
                        invalid_trials+=1
                        return 1e100,np.zeros_like(z)
                    last=dict(z=z.copy(),value=value,grad=grad,state=state)
                    if evaluations==1:history.append(dict(iteration=0,seconds=time.perf_counter()-start,elbo=float(value)))
                    return -value/(model.n*model.d),-grad/(model.n*model.d)
                def callback(z):
                    history.append(dict(iteration=len(history),seconds=time.perf_counter()-start,elbo=float(last['value'])))
                bounds=[(-40.,40.)]*(model.n*model.K)
                if not fixed:bounds += [(np.log(1e-10),np.log(1e10))]*(2*model.d)
                opt=minimize(objective,z0,method='L-BFGS-B',jac=True,bounds=bounds,callback=callback,
                             options=dict(maxiter=3000,maxfun=6000,ftol=1e-10,gtol=1e-8,maxls=40,maxcor=15))
                elapsed=time.perf_counter()-start;timings.append(elapsed)
                value,grad,state=model.profile_value_gradient(opt.x,fixed=fixed)
                result=dict(R=state['R'],sigma=state['sigma'],lam=state['lam'],pi=state['pi'],value=float(value),
                 iterations=int(opt.nit),evaluations=int(opt.nfev),converged=bool(opt.success),status=str(opt.message),
                 elapsed=elapsed,history=history,gradient_inf=float(np.max(abs(grad))),invalid_trials=invalid_trials)
            print(tag,method,'repeat',repeat+1,'rho',rho(ref['truth'],result['R']),'ELBO',result['value'],
                  'seconds',timings[-1],'iterations',result['iterations'],flush=True)
        diagnostic=cavi(model,result['R'],fixed=fixed,maxiter=1,tol=0,state=result)
        result['one_cavi_sweep_gain']=diagnostic['value']-result['value']
        record={k:v for k,v in result.items() if k not in ('R','sigma','lam','pi','states')}
        record.update(case=case,mode=modes,epsilon=epsilon,method=method,rho=rho(ref['truth'],result['R']),
          timing_repeats=timings,median_seconds=float(np.median(timings)),initial_elbo=float(result['history'][0]['elbo']),
          validation=validation,reference_sha256=hashlib.sha256((ROOT/f'{case}_reference.json').read_bytes()).hexdigest(),
          script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest())
        np.savez_compressed(ROOT/f'{tag}_{method}.npz',R=result['R'],sigma=result['sigma'],lam=result['lam'],pi=result['pi'])
        (ROOT/f'{tag}_{method}.json').write_text(json.dumps(record,indent=2)+'\n');outputs.append(record)
    return outputs

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('case');parser.add_argument('--epsilon',type=float,default=.01)
    parser.add_argument('--fixed',action='store_true');parser.add_argument('--repeats',type=int,default=3)
    args=parser.parse_args();run(args.case,args.epsilon,args.fixed,args.repeats)
