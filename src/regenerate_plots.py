from pathlib import Path
import json
import numpy as np
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import analysis
out=ROOT/'results/refined_20261104'
runs={}
for kind in ['circular','lowe','uniform']:
    full=np.load(out/f'{kind}_chain.npz')['chain']
    chain=full[len(full)//2:]
    flat=np.load(out/f'{kind}_samples.npy')
    diag=json.loads((out/f'{kind}_diagnostics.json').read_text())
    post=analysis.Posterior(kind,diag['pmax'])
    runs[kind]=(flat,chain,diag,post)
analysis.make_outputs(runs,out,137010)
print('Regenerated titled figures from saved posterior samples.')
