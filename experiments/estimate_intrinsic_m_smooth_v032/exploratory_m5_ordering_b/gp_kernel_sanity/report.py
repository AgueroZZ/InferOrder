"""Summarize recovery, agreement, numerical validation, and approximation error."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,csv,hashlib,sys
import numpy as np
import scipy,GPy
from scipy.stats import rankdata,spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
rows=[];refs={};checks={}
for case in ['B_noise0','B_noise05','B_noise1']:
 ref=json.loads((BASE/'collapsed_comparison'/f'{case}_reference.json').read_text())
 bgpath=next((BASE/'external_methods/gpy_results').glob(f'*_{case}_BayesianGPLVM_pca_m50_s20260929.json'))
 bg=json.loads(bgpath.read_text());refs[case]=(ref,bg)
 rows.append(dict(case=case,method='Native RW2',rho=abs(spearmanr(ref['truth'],np.array(ref['native_final']['R'])@np.arange(50)).statistic),agreement_bgplvm=abs(spearmanr(bg['final'],np.array(ref['native_final']['R'])@np.arange(50)).statistic)))
 rows.append(dict(case=case,method='Bayesian GPLVM (50 inducing)',rho=bg['final_rho'],agreement_bgplvm=1.))
 for f in sorted(ROOT.glob(f'{case}_*.json')):
  if 'audit' in f.name:continue
  r=json.loads(f.read_text());z=np.load(f.with_suffix('.npz'));final=z['final']
  rows.append(dict(case=case,method=f.stem[len(case)+1:],rho=r['rho'],agreement_bgplvm=abs(spearmanr(bg['final'],final).statistic)))
  history=r['history'];checks[f.stem]=dict(min_elbo_increment=float(np.min(np.diff([h['elbo'] for h in history]))) if len(history)>1 else None,
     converged=r['converged'],initial_final_rho=abs(spearmanr(z['initial'],final).statistic))
with (ROOT/'summary.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['case','method','rho','agreement_bgplvm']);w.writeheader();w.writerows(rows)
(ROOT/'fit_checks.json').write_text(json.dumps(checks,indent=2))
fig,axes=plt.subplots(2,3,figsize=(12,7.5),sharex=True,sharey=True)
for col,case in enumerate(refs):
 ref,bg=refs[case];truth=np.array(ref['truth'])
 for row,tag in enumerate([f'{case}_swap_rbf_m0',f'{case}_aligned_rbf_K100']):
  z=np.load(ROOT/(tag+'.npz'));r=json.loads((ROOT/(tag+'.json')).read_text());ax=axes[row,col]
  for values,label,color,size in [(z['initial'],'CAVI initial PCA','#999999',9),(bg['final'],'Bayesian GPLVM','#0072B2',12),(z['final'],'GP-kernel CAVI','#D55E00',11)]:
   rank=rankdata(values)/len(values)
   orientation=spearmanr(truth,bg['final']).statistic * spearmanr(bg['final'],values).statistic
   if orientation<0:rank=1-rank
   ax.scatter(truth,rank,s=size,color=color,alpha=.55,label=label)
  ax.set_title(f"{['No noise','Half noise','Original noise'][col]}\nTruth rho: CAVI {r['rho']:.3f}, BGPLVM {bg['final_rho']:.3f}")
  if row==1:ax.set_xlabel('True latent position')
  if col==0:ax.set_ylabel(['Replace Q only','Align GP model settings'][row]+'\nEstimated rank / n')
axes[0,0].legend(fontsize=8);fig.suptitle('Ordering B: GP covariance in CAVI versus Bayesian GPLVM',fontsize=14);fig.tight_layout();fig.savefig(ROOT/'ordering_comparison.png',dpi=170);plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(12,3.8))
for col,case in enumerate(refs):
 ref,bg=refs[case];z=np.load(ROOT/f'{case}_aligned_rbf_K100.npz');x=rankdata(bg['final'])/300;y=rankdata(z['final'])/300
 if spearmanr(x,y).statistic<0:y=1-y
 ax=axes[col];points=ax.scatter(x,y,c=ref['truth'],cmap='viridis',vmin=0,vmax=1,s=12);ax.plot([0,1],[0,1],color='#aaaaaa',linestyle='--');ax.set_xlabel('Bayesian GPLVM rank / n');ax.set_ylabel('GP-kernel CAVI rank / n');ax.set_title(f"{['No noise','Half noise','Original noise'][col]} | agreement {abs(spearmanr(x,y).statistic):.3f}")
fig.tight_layout(rect=(0,0,.91,1));cax=fig.add_axes([.93,.2,.015,.6]);fig.colorbar(points,cax=cax,label='True latent position');fig.savefig(ROOT/'ordering_agreement.png',dpi=170);plt.close(fig)
provenance=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,GPy=GPy.__version__,
 seed=20260929,jitter_relative_to_kernel_amplitude=1e-6,package_commit='15f2b0bbe5dfa61cd46da5160b2bc251e75a0475',
 inputs={str(f.relative_to(BASE)):hashlib.sha256(f.read_bytes()).hexdigest() for f in (BASE/'external_methods/inputs').glob('B_*csv')},
 comparator_files=[str(next((BASE/'external_methods/gpy_results').glob(f'*_{case}_BayesianGPLVM_pca_m50_s20260929.json')).relative_to(BASE)) for case in refs])
(ROOT/'provenance.json').write_text(json.dumps(provenance,indent=2));print(json.dumps(rows,indent=2));print('Minimum ELBO increment',min(x['min_elbo_increment'] for x in checks.values()))
