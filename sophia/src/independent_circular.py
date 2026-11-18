"""Independent circular ensemble from a separate high-period starting region."""
from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import analysis as a

out=ROOT/'results/independent_circular_20261104'
out.mkdir(exist_ok=False)
post=a.Posterior('circular')
rng=np.random.default_rng(941370)
center=np.array([np.log(650),.65,np.log(.014),.004,.462,.60,.12])
spread=np.array([.3,.15,.04,.008,.04,.05,.05])
pos=[]
while len(pos)<48:
    x=center+spread*rng.normal(size=7)
    if np.isfinite(post(x)):pos.append(x)
sampler=a.emcee.EnsembleSampler(48,7,post,moves=[(a.emcee.moves.StretchMove(),.8),(a.emcee.moves.DEMove(),.2)])
sampler.random_state=np.random.RandomState(941370).get_state()
start=time.time()
for i in range(6):
    sampler.run_mcmc(pos if i==0 else None,10000,progress=False)
    np.savez_compressed(out/'circular_chain.npz',chain=sampler.get_chain(),log_probability=sampler.get_log_prob(),acceptance=sampler.acceptance_fraction,random_state=np.array(sampler.random_state,dtype=object))
    print('Independent circular',sampler.iteration,flush=True)
chain=sampler.get_chain(discard=30000)
flat=chain.reshape(-1,7)
np.save(out/'circular_samples.npy',flat)
diag=a.diagnostics(chain,sampler.acceptance_fraction)
diag.update(seed=941370,steps=60000,walkers=48,seconds=time.time()-start,start_center=center.tolist(),physical_period_quantiles=a.quantiles(np.exp(flat[:,0])),independent_of_refinement=True)
(out/'diagnostics.json').write_text(json.dumps(diag,indent=2))
print(json.dumps(diag,indent=2),flush=True)
