from pathlib import Path
import json,numpy as np
root=Path(__file__).resolve().parent
for f in root.glob('*_swap_*.npz'):
 z=np.load(f);(root/(f.stem+'_audit_input.json')).write_text(json.dumps(dict(Q=z['precision'].tolist(),R=z['R'].tolist())))
