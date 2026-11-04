"""Reproducible HD 137010 b prior experiment; no legacy samples are overwritten.

Run from repository root: .venv/Scripts/python src/analysis.py --steps 20000
The three models differ only in their explicitly conditional eccentricity prior.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np
import scipy
from scipy.stats import ks_2samp, wasserstein_distance
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
import batman
import emcee

ROOT = Path(__file__).resolve().parents[1]
TREF = 2458055.0969
G = 6.67430e-11
LABELS = ['ln P', 'b', 'ln radius ratio', 't0 offset [d]', 'log10 density', 'u1', 'u2', 'sqrt(e) cos(w)', 'sqrt(e) sin(w)']

def scaled_axis(period, density):
    return (G * density * 1000 * (period * 86400)**2 / (3*np.pi))**(1/3)

def unpack(x):
    p, b, r, dt, lr, u1, u2 = x[:7]
    e = np.sum(np.asarray(x[7:9])**2) if len(x) == 9 else 0.
    w = np.arctan2(x[8], x[7]) if len(x) == 9 else np.pi/2
    return np.exp(p), b, np.exp(r), TREF+dt, 10**lr, u1, u2, e, w

def log_prior(x, kind, pmax=2000.):
    if not np.all(np.isfinite(x)):
        return -np.inf
    p,b,r,t0,rho,u1,u2,e,w = unpack(x)
    if not (90 < p < pmax and 0 <= b < 1 and .005 < r < .03 and abs(t0-TREF)<.1 and 1<rho<5):
        return -np.inf
    if not (0 <= e < .8 and u1>0 and u1+u2<1 and u1+2*u2>0):
        return -np.inf
    a = scaled_axis(p,rho)
    if a*(1-e) <= 1+r:
        return -np.inf
    # Density measure is d(log10 rho); log-period sampling supplies Jacobian P.
    # Thus a physical P^(-5/3) density becomes P^(-2/3) in ln P.
    lp = -2/3*np.log(p) - .5*((x[4]-.462)/.041)**2
    lp -= .5*((u1-.60)/.07)**2 + .5*((u2-.12)/.07)**2
    # Flat ln(radius ratio). Cartesian sqrt(e) coordinates give flat de dw.
    if kind == 'lowe':
        lp -= e/.06
    return lp

class Transit:
    def __init__(self, t):
        self.t = np.asarray(t)
        self.params = batman.TransitParams()
        q = self.params
        q.t0=TREF; q.per=355.; q.rp=.01368; q.a=270.; q.inc=90.
        q.ecc=0.; q.w=90.; q.u=[.56,.10]; q.limb_dark='quadratic'
        self.engine = batman.TransitModel(q,self.t,supersample_factor=10,exp_time=29.4/1440)

    def __call__(self,x):
        p,b,r,t0,rho,u1,u2,e,w=unpack(x)
        a=scaled_axis(p,rho)
        ci=b/a*(1+e*np.sin(w))/(1-e*e)
        if not 0 <= ci <= 1:
            return None
        q=self.params
        q.t0=t0; q.per=p; q.rp=r; q.a=a; q.inc=np.degrees(np.arccos(ci))
        q.ecc=e; q.w=np.degrees(w); q.u=[u1,u2]
        return self.engine.light_curve(q)

class Posterior:
    def __init__(self,kind,pmax=2000.):
        d=np.loadtxt(ROOT/'data/apjladf06ft3_mrt.txt',skiprows=15)
        t=d[:,0]+2454833
        mask=np.isfinite(d[:,2]) & (abs(t-TREF)<.8)
        self.t=t[mask]; self.flux=d[mask,2]
        oot=abs(self.t-TREF)>.3
        self.sigma=np.std(self.flux[oot],ddof=1)
        self.kind=kind; self.pmax=pmax
        self.model=Transit(self.t)
        # Marginalize a constant baseline with an improper flat prior.
        # With fixed sigma, the missing determinant is parameter-independent.
    def __call__(self,x):
        lp=log_prior(x,self.kind,self.pmax)
        if not np.isfinite(lp): return -np.inf
        m=self.model(x)
        if m is None or not np.all(np.isfinite(m)): return -np.inf
        res=self.flux-m
        res-=np.mean(res)
        return lp-.5*np.sum((res/self.sigma)**2)

def duration(x):
    p,b,r,t0,rho,u1,u2,e,w=unpack(x)
    a=scaled_axis(p,rho)
    ci=b/a*(1+e*np.sin(w))/(1-e*e)
    arg=np.sqrt(max(0,(1+r)**2-b*b))/(a*np.sqrt(1-ci*ci))
    return 24*p/np.pi*np.arcsin(np.clip(arg,0,1))*np.sqrt(1-e*e)/(1+e*np.sin(w))

def prior_draws(kind,n,rng,pmax=2000.):
    # Inverse CDF for physical period density P^-5/3.
    u=rng.random(n)
    p=(90**(-2/3)+u*(pmax**(-2/3)-90**(-2/3)))**(-1.5)
    base=np.column_stack([np.log(p),rng.uniform(0,1,n),rng.uniform(np.log(.005),np.log(.03),n),rng.uniform(-.1,.1,n),rng.normal(.462,.041,n),rng.normal(.60,.07,n),rng.normal(.12,.07,n)])
    if kind!='circular':
        e=( -.06*np.log1p(-rng.random(n)*(1-np.exp(-.8/.06))) if kind=='lowe' else rng.uniform(0,.8,n))
        w=rng.uniform(-np.pi,np.pi,n)
        base=np.column_stack([base,np.sqrt(e)*np.cos(w),np.sqrt(e)*np.sin(w)])
    return np.array([x for x in base if np.isfinite(log_prior(x,kind,pmax))])

def diagnostics(chain,acceptance):
    n=chain.shape[0]
    try:
        tau=emcee.autocorr.integrated_time(chain,tol=0)
    except Exception:
        tau=np.full(chain.shape[-1],np.nan)
    # Split-chain Rhat over ensemble walkers is descriptive only: walkers interact.
    h=n//2
    split=np.concatenate([chain[:h],chain[-h:]],axis=1)
    within=np.var(split,axis=0,ddof=1).mean(axis=0)
    between=h*np.var(split.mean(axis=0),axis=0,ddof=1)
    rhat=np.sqrt(((h-1)/h*within+between/h)/within)
    enough=bool(np.all(np.isfinite(tau)) and np.all(n>50*tau) and np.all(rhat<1.05))
    return dict(tau=tau.tolist(),approx_ess=(n*chain.shape[1]/tau).tolist(),split_rhat_descriptive=rhat.tolist(),retained_steps=n,acceptance_mean=float(np.mean(acceptance)),convergence_screen_passed=enough,warning='Ensemble walkers interact; this Rhat is not an independent-chain convergence proof. Run independent ensembles and longer chains for final inference.')

def sample(kind,steps,walkers,seed,pmax,out):
    rng=np.random.default_rng(seed)
    post=Posterior(kind,pmax)
    center=np.array([np.log(350),.25,np.log(.0137),0.,.462,.56,.10]+([0.,0.] if kind!='circular' else []))
    spread=np.array([.15,.15,.03,.004,.03,.04,.04]+([.2,.2] if kind!='circular' else []))
    pos=[]
    while len(pos)<walkers:
        x=center+spread*rng.normal(size=len(center))
        if np.isfinite(post(x)): pos.append(x)
    sampler=emcee.EnsembleSampler(walkers,len(center),post)
    sampler.random_state=np.random.RandomState(seed).get_state()
    start=time.time()
    # Checkpoint complete history in blocks, without thinning away diagnostics.
    for end in range(0,steps,1000):
        state=sampler.run_mcmc(pos if end==0 else None,min(1000,steps-end),progress=False)
        np.savez_compressed(out/f'{kind}_chain.npz',chain=sampler.get_chain(),log_probability=sampler.get_log_prob(),acceptance=sampler.acceptance_fraction,random_state=np.array(sampler.random_state,dtype=object))
        print(f'{kind}: {sampler.iteration}/{steps} steps; {time.time()-start:.0f}s',flush=True)
    chain=sampler.get_chain(discard=steps//2)
    flat=chain.reshape(-1,len(center))
    np.save(out/f'{kind}_samples.npy',flat)
    diag=diagnostics(chain,sampler.acceptance_fraction)
    diag.update(seed=seed,steps=steps,walkers=walkers,seconds=time.time()-start,sigma_flux=post.sigma,pmax=pmax)
    (out/f'{kind}_diagnostics.json').write_text(json.dumps(diag,indent=2))
    return flat,chain,diag,post

def quantiles(x): return [float(v) for v in np.percentile(x,[16,50,84])]

def derived(flat,rng):
    # Draw mass independently of transit density, derive radius consistently from M/rho.
    mass=10**rng.normal(-.139,.010,len(flat))
    rho=10**flat[:,4]
    radius=(mass*1.408/rho)**(1/3)
    teff=rng.normal(4770,90,len(flat))
    luminosity=radius**2*(teff/5772)**4
    p=np.exp(flat[:,0]); e=np.sum(flat[:,7:9]**2,axis=1) if flat.shape[1]==9 else np.zeros(len(flat))
    a=(mass*(p/365.25)**2)**(1/3)
    s=luminosity/a**2
    return dict(period_days=p,eccentricity=e,semi_major_axis_AU=a,radius_Earth=np.exp(flat[:,2])*radius*109.076,
                flux_at_a_Earth=s,mean_flux_Earth=s/np.sqrt(1-e*e),periastron_AU=a*(1-e),apastron_AU=a*(1+e),
                periastron_flux_Earth=s/(1-e)**2,apastron_flux_Earth=s/(1+e)**2,flux_ratio=((1+e)/(1-e))**2)

def savefig(out,name):
    plt.tight_layout()
    plt.savefig(out/f'{name}.png',dpi=180)
    plt.savefig(out/f'{name}.svg')
    plt.close()

def make_outputs(runs,out,seed):
    rng=np.random.default_rng(seed)
    colors=['#2864a0','#269a81','#cc6b39']
    summaries={}; physical={}
    for kind,(flat,chain,diag,post) in runs.items():
        d=derived(flat,rng); physical[kind]=d
        np.savez_compressed(out/f'{kind}_derived.npz',**d)
        summaries[kind]={k:quantiles(v) for k,v in d.items()}
        summaries[kind]['period_near_lower_bound_fraction']=float(np.mean(d['period_days']<99))
        summaries[kind]['period_near_upper_bound_fraction']=float(np.mean(d['period_days']>.9*diag['pmax']))
        summaries[kind]['convergence_screen_passed']=diag['convergence_screen_passed']
        fig,axes=plt.subplots(flat.shape[1],1,figsize=(10,12),sharex=True)
        stride=max(1,len(chain)//1000)
        for j,ax in enumerate(axes):
            ax.plot(np.arange(0,len(chain),stride),chain[::stride,:,j],alpha=.15,lw=.4); ax.set_ylabel(LABELS[j],fontsize=8)
        axes[-1].set_xlabel('Retained MCMC step (burn-in removed; display subsampled)')
        savefig(out,f'{kind}_traces')
    # Distribution comparison, using all retained samples in numerical summaries.
    fig,axs=plt.subplots(1,3,figsize=(13,3.8))
    names={'circular':'Circular','lowe':'Low eccentricity','uniform':'Broad eccentricity'}
    for (kind,d),color in zip(physical.items(),colors):
        for ax,key in zip(axs,['period_days','semi_major_axis_AU','mean_flux_Earth']):
            vals=d[key]; bins=np.geomspace(max(vals.min()*.9,1e-4),vals.max()*1.1,80)
            ax.hist(vals,bins=bins,density=True,histtype='step',color=color,label=names[kind],lw=1.5)
            ax.set_xscale('log')
    for ax,label in zip(axs,['Period [d]','Semi-major axis [AU]','Orbit-mean flux [Earth]']):
        ax.set_xlabel(label); ax.set_ylabel('Density')
        ax.xaxis.set_minor_formatter(NullFormatter())
    for ax,ticks in zip(axs,[[100,300,1000],[.4,.7,1,2],[.03,.1,.3,1,3]]):
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:g}'))
    axs[0].legend(); savefig(out,'prior_sensitivity')
    plt.figure(figsize=(8,4))
    for kind,color in zip(runs,colors):
        prior=prior_draws(kind,15000,rng)
        dur=np.array([duration(x) for x in prior])
        plt.hist(dur,bins=np.linspace(0,40,100),weights=np.ones(len(dur))/len(dur)/(40/99),histtype='step',label=kind,color=color)
    plt.axvline(9.76,color='black',ls='--',label='Published duration')
    plt.xlabel('Approximate first-to-fourth contact duration [h]');plt.ylabel('Density (display restricted to 0–40 h)');plt.legend()
    savefig(out,'prior_predictive_duration')
    # Compare circular ECDF on complete support, not a truncated histogram.
    import csv
    with (ROOT/'data/posteriors.csv').open() as f:
        ref=np.array([float(row['P']) for row in csv.DictReader(f)])
    ours=physical['circular']['period_days']
    plt.figure(figsize=(8,4))
    for vals,label in [(ref,'Venner supplied samples'),(ours,'New circular run')]:
        v=np.sort(vals);plt.plot(v,np.arange(1,len(v)+1)/len(v),label=label)
    plt.xscale('log');plt.xlabel('Period [d]');plt.ylabel('Empirical cumulative probability');plt.legend()
    savefig(out,'circular_validation')
    validation=dict(reference_quantiles=quantiles(ref),new_quantiles=quantiles(ours),ks_distance=float(ks_2samp(ref,ours).statistic),wasserstein_days=float(wasserstein_distance(ref,ours)),warning='Distances are descriptive; do not interpret iid KS p-values for correlated MCMC samples. This is an approximate reproduction, not identical noise/detrending/selection modelling.')
    (out/'validation.json').write_text(json.dumps(validation,indent=2))
    for kind,(flat,chain,diag,post) in runs.items():
        take=rng.choice(len(flat),min(300,len(flat)),replace=False)
        curves=np.array([post.model(flat[i]) for i in take])
        curves+=np.mean(post.flux-curves,axis=1)[:,None]
        med=np.median(curves,axis=0)
        fig,(ax,ar)=plt.subplots(2,1,figsize=(9,6),sharex=True,gridspec_kw={'height_ratios':[2,1]})
        hours=(post.t-TREF)*24
        lo,hi=np.percentile(curves,[2.5,97.5],axis=0)
        ax.plot(hours,(post.flux-1)*1e6,'.',color='.5'); ax.fill_between(hours,(lo-1)*1e6,(hi-1)*1e6,alpha=.3);ax.plot(hours,(med-1)*1e6)
        ax.set_ylabel('Flux relative to unity [ppm]');ax.set_title(kind+' posterior model band (95%; excludes measurement noise)')
        residual=(post.flux-med)/post.sigma
        ar.plot(hours,residual,'.');ar.axhline(0,color='k',lw=.5);ar.set_ylabel('Residual / sigma');ar.set_xlabel('Time from reference transit [h]')
        savefig(out,f'{kind}_model_check')
        # Posterior replicated discrepancy including measurement noise.
        observed=np.sum(((post.flux[None,:]-curves)/post.sigma)**2,axis=1)
        noise=rng.normal(size=curves.shape); noise-=noise.mean(axis=1)[:,None]
        replicated=np.sum(noise**2,axis=1)
        check=dict(posterior_predictive_chi_square_tail=float(np.mean(replicated>=observed)),median_residual_lag1=float(np.corrcoef(residual[:-1],residual[1:])[0,1]),rms_ppm=float(np.std(post.flux-med)*1e6),warning='Fixed white-noise sigma estimated from local out-of-transit flux; no red-noise or detrending uncertainty model.')
        (out/f'{kind}_model_check.json').write_text(json.dumps(check,indent=2))
    if 'uniform' in runs:
        flat=runs['uniform'][0]
        e=physical['uniform']['eccentricity']; p=physical['uniform']['period_days']
        plt.figure(figsize=(8,4));sel=rng.choice(len(flat),min(12000,len(flat)),replace=False)
        plt.scatter(p[sel],e[sel],c=flat[sel,1],s=2,alpha=.25,cmap='viridis');plt.xscale('log');plt.colorbar(label='Impact parameter b');plt.xlabel('Period [d]');plt.ylabel('Eccentricity');savefig(out,'period_eccentricity')
        post=runs['uniform'][3];grid=Transit(np.linspace(TREF-.35,TREF+.35,450))
        # Select real high-posterior samples from distinct period quantile bands.
        scores=np.array([post(x) for x in flat[sel]])
        fig,ax=plt.subplots(figsize=(9,4));ax.plot((post.t-TREF)*24,(post.flux-1)*1e6,'.',color='.7',label='K2')
        examples=[]
        for q in [.1,.5,.9]:
            target=np.quantile(p,q); candidates=np.where(abs(np.log(p[sel]/target))<.08)[0]
            if len(candidates)==0: continue
            idx=sel[candidates[np.argmax(scores[candidates])]];x=flat[idx]
            model=grid(x)+np.mean(post.flux-post.model(x))
            ax.plot((grid.t-TREF)*24,(model-1)*1e6,label=f'P={p[idx]:.0f} d, e={e[idx]:.2f}, b={x[1]:.2f}')
            examples.append(dict(zip(['P','b','radius_ratio','t0','rho','u1','u2','e','omega_radians'],map(float,unpack(x)))))
        ax.set_xlabel('Time from reference transit [h]');ax.set_ylabel('Flux [ppm]');ax.legend(fontsize=8);savefig(out,'representative_transits')
        (out/'representative_transits.json').write_text(json.dumps(examples,indent=2))
    (out/'summary.json').write_text(json.dumps(summaries,indent=2))
    rows=['model,quantity,q16,q50,q84']
    for kind,d in physical.items():
        rows += [f'{kind},{key},'+','.join(map(str,quantiles(v))) for key,v in d.items()]
    (out/'summary.csv').write_text('\n'.join(rows)+'\n')
    return summaries

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--steps',type=int,default=20000);ap.add_argument('--walkers',type=int,default=48);ap.add_argument('--seed',type=int,default=137010);ap.add_argument('--pmax',type=float,default=2000);ap.add_argument('--output',default='results/reproducible_v1');args=ap.parse_args()
    out=ROOT/args.output
    out.mkdir(parents=True,exist_ok=False)
    manifest=dict(arguments=vars(args),parameter_order=LABELS,time_reference_BJD=TREF,period_bounds_days=[90,args.pmax],eccentricity_max=.8,conditioning='Priors specified for the observed transiting candidate; no additional e,omega transit-probability weight. Not an intrinsic population model.',versions={k:v for k,v in [('numpy',np.__version__),('scipy',scipy.__version__),('emcee',emcee.__version__),('batman',batman.__version__)]},inputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'data/apjladf06ft3_mrt.txt',ROOT/'data/posteriors.csv',Path(__file__)]})
    manifest['parameter_order_by_model']={k:LABELS[:7] if k=='circular' else LABELS for k in ['circular','lowe','uniform']}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    runs={k:sample(k,args.steps,args.walkers,args.seed+i,args.pmax,out) for i,k in enumerate(['circular','lowe','uniform'])}
    make_outputs(runs,out,args.seed)
    print('Finished:',out,flush=True)

if __name__=='__main__': main()
