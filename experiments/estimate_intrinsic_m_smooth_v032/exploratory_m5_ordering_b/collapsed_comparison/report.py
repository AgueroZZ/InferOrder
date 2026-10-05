import sys
from pathlib import Path
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root))
import compare as c
import json,csv,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import rankdata
records=[json.loads(p.read_text()) for p in sorted(root.glob('*_eps*.json'))]
rows=[];audit=[]
for r in records:
 target=json.loads((root/f"{r['case']}_{r['mode']}_eps{r['epsilon']:g}_cavi_python.json").read_text())['value']-.0036
 hits=[h['seconds'] for h in r['history'] if h['elbo']>=target]
 row={k:r.get(k) for k in ['case','mode','epsilon','method','rho','value','iterations','evaluations','median_seconds','converged','one_cavi_sweep_gain']}
 row['seconds_to_cavi_target']=hits[0] if hits else None;rows.append(row)
 if r['mode']=='adaptive' and r['epsilon']==.01 and r['method']!='collapsed_lbfgs':
  z=np.load(root/f"{r['case']}_adaptive_eps0.01_{r['method']}.npz")
  audit.append(dict(case=r['case'],method=r['method'],R=z['R'].tolist(),sigma2=z['sigma'].tolist(),lambda_=z['lam'].tolist(),pi=z['pi'].tolist(),expected=r['value']))
(root/'endpoint_inputs.json').write_text(json.dumps(audit))
with (root/'summary.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
checks={}
fig,axes=plt.subplots(3,3,figsize=(11,10),sharex=True,sharey=True)
for col,case in enumerate(['B_noise0','B_noise05','B_noise1']):
 ref=json.loads((root/f'{case}_reference.json').read_text());warm=json.loads((root/f'{case}_warm_reference.json').read_text())
 py=np.load(root/f'{case}_adaptive_eps0.01_cavi_python.npz')
 pr=json.loads((root/f'{case}_adaptive_eps0.01_cavi_python.json').read_text())
 checks[case]=dict(max_R_error=float(np.max(abs(py['R']-warm['R']))),elbo_error=abs(pr['value']-warm['elbo']),native_rho=c.rho(ref['truth'],np.array(ref['native_final']['R'])),native_elbo=ref['native_final']['elbo'],native_seconds=float(np.median(ref['native_timing_seconds'])))
 truth=np.array(ref['truth']);init=rankdata(ref['pca'])/len(truth)
 for row,method in enumerate(['native','cavi_python','collapsed_sqrt']):
  R=np.array(ref['native_final']['R']) if method=='native' else np.load(root/f'{case}_adaptive_eps0.01_{method}.npz')['R']
  est=R@np.linspace(0,1,R.shape[1]); est=rankdata(est)/len(truth)
  if np.corrcoef(est,truth)[0,1]<0:est=1-est
  initial=init if np.corrcoef(init,truth)[0,1]>=0 else 1-init
  ax=axes[row,col];ax.scatter(truth,initial,s=8,c='#aaaaaa',alpha=.5,label='Initial PCA');ax.scatter(truth,est,s=9,c='#0072B2',alpha=.65,label='Final')
  ax.set_title(f"{['No noise','Half noise','Original noise'][col]} | |Spearman rho| = {c.rho(truth,R):.3f}")
  if col==0:ax.set_ylabel(['Native default CAVI','Matched-start CAVI','Collapsed + L-BFGS'][row]+'\nEstimated rank / n')
  if row==2:ax.set_xlabel('True latent position')
axes[0,0].legend(fontsize=8);fig.suptitle('Ordering B: profiling the curve block does not unfold the PCA ordering',fontsize=13);fig.tight_layout();fig.savefig(root/'ordering_comparison.png',dpi=170);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(11,4))
for method,label,color in [('cavi_python','CAVI','#0072B2'),('collapsed_sqrt','Collapsed (square coordinates)','#D55E00'),('collapsed_lbfgs','Collapsed (logits; stalled)','#999999')]:
 r=json.loads((root/f'B_noise05_adaptive_eps0.01_{method}.json').read_text());h=r['history']
 for ax in axes:ax.plot([max(x['seconds'],1e-5) for x in h],[x['elbo'] for x in h],label=label,color=color);ax.set_xscale('log');ax.set_xlabel('Elapsed seconds (single thread)');ax.set_ylabel('MPCurve ELBO')
axes[0].set_title('Half noise: full optimization trace');axes[1].set_title('Half noise: terminal objective region');axes[1].set_ylim(-1800,-1400);axes[0].legend(fontsize=8);fig.tight_layout();fig.savefig(root/'convergence_half_noise.png',dpi=170);plt.close(fig)
(root/'reference_checks.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2));print(json.dumps([r for r in rows if r['mode']=='adaptive' and r['epsilon']==.01],indent=2))
