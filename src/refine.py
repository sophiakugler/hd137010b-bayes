"""Extend trusted local checkpoints; preserve the pilot run in a separate folder."""
from pathlib import Path
import sys,json,time,argparse
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import analysis as a

parser=argparse.ArgumentParser()
parser.add_argument('--source',default='results/reproducible_20261104')
parser.add_argument('--output',default='results/refined_20261104')
parser.add_argument('--min-steps',type=int,default=60000)
parser.add_argument('--max-steps',type=int,default=120000)
args=parser.parse_args()
source=ROOT/args.source
out=ROOT/args.output
out.mkdir(exist_ok=False)
manifest=json.loads((source/'manifest.json').read_text())
manifest.update(continuation_from=str(source.relative_to(ROOT)),min_total_steps=args.min_steps,max_total_steps=args.max_steps,stopping='Retained half must exceed 50 estimated autocorrelation times in every coordinate and descriptive split Rhat < 1.05. Stop at requested maximum and disclose failure if not met.',moves='80% StretchMove + 20% DEMove during continuation')
(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
runs={}
for kind in ['circular','lowe','uniform']:
    # allow_pickle is only for RNG state in our own trusted generated checkpoint.
    data=np.load(source/f'{kind}_chain.npz',allow_pickle=True)
    old=data['chain'];oldlog=data['log_probability'];rngstate=tuple(data['random_state'])
    post=a.Posterior(kind)
    initial_diag=a.diagnostics(old[len(old)//2:],data['acceptance'])
    if (len(old)>=args.min_steps and initial_diag['convergence_screen_passed']) or len(old)>=args.max_steps:
        import shutil
        shutil.copy2(source/f'{kind}_chain.npz',out/f'{kind}_chain.npz')
        diag=json.loads((source/f'{kind}_diagnostics.json').read_text())
        retained=old[len(old)//2:]
        flat=retained.reshape(-1,old.shape[-1])
        np.save(out/f'{kind}_samples.npy',flat)
        (out/f'{kind}_diagnostics.json').write_text(json.dumps(diag,indent=2))
        runs[kind]=(flat,retained,diag,post)
        print(kind,'retained existing checkpoint',len(old),flush=True)
        continue
    sampler=a.emcee.EnsembleSampler(old.shape[1],old.shape[2],post,moves=[(a.emcee.moves.StretchMove(),.8),(a.emcee.moves.DEMove(),.2)])
    sampler.random_state=rngstate
    state=a.emcee.State(old[-1],log_prob=oldlog[-1],random_state=rngstate)
    start=time.time()
    for block in range((args.max_steps-len(old)+9999)//10000):
        sampler.run_mcmc(state if block==0 else None,min(10000,args.max_steps-len(old)-sampler.iteration),progress=False)
        full=np.concatenate([old,sampler.get_chain()],axis=0)
        logp=np.concatenate([oldlog,sampler.get_log_prob()],axis=0)
        retained=full[len(full)//2:]
        diag=a.diagnostics(retained,sampler.acceptance_fraction)
        diag.update(seed=json.loads((source/f'{kind}_diagnostics.json').read_text())['seed'],steps=len(full),walkers=full.shape[1],seconds_continuation=time.time()-start,sigma_flux=post.sigma,pmax=2000,acceptance_scope='Continuation only')
        np.savez_compressed(out/f'{kind}_chain.npz',chain=full,log_probability=logp,acceptance=sampler.acceptance_fraction,random_state=np.array(sampler.random_state,dtype=object))
        (out/f'{kind}_diagnostics.json').write_text(json.dumps(diag,indent=2))
        print(kind,len(full),'max tau',round(max(diag['tau'])),'max descriptive Rhat',round(max(diag['split_rhat_descriptive']),3),'screen',diag['convergence_screen_passed'],flush=True)
        if len(full)>=args.min_steps and diag['convergence_screen_passed']:break
    flat=retained.reshape(-1,full.shape[-1])
    np.save(out/f'{kind}_samples.npy',flat)
    runs[kind]=(flat,retained,diag,post)
    del data,old,oldlog,full,logp,sampler
a.make_outputs(runs,out,137010)
print('Refinement complete',flush=True)
