"""Sample the circular-orbit single-transit posterior with emcee.

    python run_circular.py --period-prior loguniform_kipping

Runs in chunks until the chain is longer than 50 integrated autocorrelation
times and tau has changed by < 5 % (Foreman-Mackey et al. 2013 heuristics),
or until --max-steps. Chains (thinned x10) are stored in results/circular_<prior>.npz.
"""

import argparse
from pathlib import Path

import emcee
import numpy as np

import transit as tr

parser = argparse.ArgumentParser()
parser.add_argument("--period-prior", default="loguniform_kipping", choices=list(tr.PERIOD_PRIORS))
parser.add_argument("--p-max", type=float, default=tr.P_MAX)
parser.add_argument("--nwalkers", type=int, default=50)
parser.add_argument("--max-steps", type=int, default=300_000)
parser.add_argument("--chunk", type=int, default=10_000)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

rng = np.random.default_rng(args.seed)
t, f, p_min = tr.load_data()
model = tr.CircularTransitModel(t, f, p_min, period_prior=args.period_prior, p_max=args.p_max)

# Initialise walkers spread along the P-b degeneracy found by the MAP fit:
# at fixed duration and density, P ~ 300 d / (1 - b^2)^(3/2).
ndim = len(tr.PARAM_NAMES)
b0 = rng.uniform(0, 0.7, args.nwalkers)
p0 = np.column_stack([
    rng.normal(3222.0972, 0.002, args.nwalkers),
    np.log(300 / (1 - b0**2) ** 1.5) + rng.normal(0, 0.05, args.nwalkers),
    b0,
    rng.normal(np.log(0.0134), 0.01, args.nwalkers),
    rng.normal(0.57, 0.02, args.nwalkers),
    rng.normal(0.11, 0.02, args.nwalkers),
    rng.normal(0.46, 0.01, args.nwalkers),
    rng.normal(-0.139, 0.005, args.nwalkers),
    rng.normal(1.0, 2e-6, args.nwalkers),
    rng.normal(np.log(3.2e-5), 0.02, args.nwalkers),
])
assert all(np.isfinite(model(p)) for p in p0)

moves = [(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)]
sampler = emcee.EnsembleSampler(args.nwalkers, ndim, model, moves=moves)

state, old_tau, tau_history = p0, np.inf, []
while sampler.iteration < args.max_steps:
    state = sampler.run_mcmc(state, args.chunk, progress=False)
    tau = sampler.get_autocorr_time(tol=0)
    n = sampler.iteration
    tau_history.append([n, *tau])
    print(f"step {n:>7d}  tau_max = {tau.max():7.1f} ({tr.PARAM_NAMES[tau.argmax()]})  "
          f"N/tau = {n / tau.max():6.1f}  acc = {sampler.acceptance_fraction.mean():.2f}", flush=True)
    # tau estimates themselves fluctuate by a few %, so require 5 % stability
    if n > 50 * tau.max() and np.all(np.abs(old_tau - tau) / tau < 0.05):
        print("converged", flush=True)
        break
    old_tau = tau

# Store every THIN-th step: tau >~ 100, so this loses no information but keeps files small.
THIN = 10
out = Path(__file__).resolve().parent / "results"
out.mkdir(exist_ok=True)
np.savez_compressed(
    out / (f"circular_{args.period_prior}"
           + ("" if args.p_max == tr.P_MAX else f"_pmax{args.p_max:.0f}") + ".npz"),
    chain=sampler.get_chain(thin=THIN),            # (steps / THIN, walkers, ndim)
    log_prob=sampler.get_log_prob(thin=THIN),
    thin=THIN,
    acceptance=sampler.acceptance_fraction,
    tau_history=np.array(tau_history),
    param_names=tr.PARAM_NAMES,
    p_max=args.p_max,
)
