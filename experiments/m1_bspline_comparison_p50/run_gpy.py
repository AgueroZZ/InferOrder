"""GPy comparison with explicit common PC1 inputs and resumable fit artifacts."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['MPLCONFIGDIR'] = '/tmp/mpcurve-gpy-matplotlib'
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
import warnings
import numpy as np
import scipy
from scipy.stats import spearmanr, kendalltau
import GPy

ROOT = Path(__file__).resolve().parent

def latent(model, bayesian):
    return np.asarray(model.X.mean.values if bayesian else model.X.values).ravel().copy()

def run(case, method):
    path = ROOT / 'results' / f'{case["id"]}_{method}.json'
    if path.exists():
        return
    seed = int(case['seed']) + 100000
    np.random.seed(seed)
    input_path = ROOT / 'inputs' / f'{case["id"]}_X.csv'
    Y = np.loadtxt(input_path, delimiter=',', skiprows=1)
    positions = np.genfromtxt(ROOT / 'inputs' / f'{case["id"]}_positions.csv',
                             delimiter=',', names=True)
    X = positions['pca'][:, None].copy()
    bayesian = method == 'BayesianGPLVM'
    result = dict(id=case['id'], method=method, seed=seed,
                  input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  versions=dict(GPy=GPy.__version__, numpy=np.__version__, scipy=scipy.__version__))
    started = time.monotonic()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            model = (GPy.models.BayesianGPLVM(Y=Y, input_dim=1, X=X, num_inducing=50)
                     if bayesian else GPy.models.GPLVM(Y=Y, input_dim=1, X=X))
            initial = latent(model, bayesian)
            initial_parameters = model.param_array.copy()
            assert np.max(np.abs(initial - positions['pca'])) < 1e-10
            blocks = []
            for block in range(3):
                model.optimize(optimizer='lbfgsb', max_iters=2000, messages=False,
                               bfgs_factor=1e7, gtol=1e-5)
                opt = model.optimization_runs[-1]
                blocks.append(dict(block=block+1, status=opt.status,
                                   evaluations=int(opt.funct_eval), objective=float(model.log_likelihood())))
                if not opt.status.startswith('Maximum'):
                    break
            final = latent(model, bayesian)
            result.update(initial=initial.tolist(), final=final.tolist(), blocks=blocks,
                          status=blocks[-1]['status'], converged=blocks[-1]['status'].startswith('Converged'),
                          rho=float(abs(spearmanr(positions['truth'], final).statistic)),
                          tau=float(abs(kendalltau(positions['truth'], final).statistic)),
                          noise_variance=float(model.likelihood.variance[0]))
            np.savez_compressed(path.with_suffix('.npz'), initial_parameters=initial_parameters,
                                final_parameters=model.param_array.copy(), initial=initial, final=final,
                                parameter_names=np.asarray(model.parameter_names()))
        except Exception as error:
            result.update(status=f'error: {error}', converged=False, rho=None, tau=None)
    result['elapsed_seconds'] = time.monotonic() - started
    result['warnings'] = sorted(set(str(w.message) for w in caught))
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(case['id'], method, result['status'], result['rho'],
          round(result['elapsed_seconds'], 2), 'seconds', flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--replications', type=int, nargs='+')
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    with (ROOT / 'manifest.csv').open() as handle:
        cases = list(csv.DictReader(handle))
    jobs = [(case, method) for case in cases
            if args.replications is None or int(case['replication']) in args.replications
            for method in ('GPLVM', 'BayesianGPLVM')]
    if args.workers == 1:
        for job in jobs:
            run(*job)
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(run, *job) for job in jobs]
            for future in futures:
                future.result()
