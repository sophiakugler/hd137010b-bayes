"""Sample the light-curve shape posterior (flat in ln v) with emcee.

    python run_shape.py

The physical posterior under any period / density / eccentricity prior is then
obtained with `transit.shape_to_physical` (importance weights). Same convergence
rule and output format as run_circular.py; output: results/shape.npz.
"""

import argparse
from pathlib import Path

import emcee
import numpy as np

import transit as tr

parser = argparse.ArgumentParser()
parser.add_argument("--nwalkers", type=int, default=50)
parser.add_argument("--max-steps", type=int, default=300_000)
parser.add_argument("--min-steps", type=int, default=100_000)
parser.add_argument("--chunk", type=int, default=10_000)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

rng = np.random.default_rng(args.seed)
t, f, _ = tr.load_data()
model = tr.ShapeTransitModel(t, f)

ndim = len(tr.SHAPE_PARAM_NAMES)
b0 = rng.uniform(0, 0.7, args.nwalkers)
p0 = np.column_stack([
    rng.normal(3222.0972, 0.002, args.nwalkers),
    np.log(0.8 * np.sqrt(1 - b0**2)) + rng.normal(0, 0.02, args.nwalkers),
    b0,
    rng.normal(np.log(0.0134), 0.01, args.nwalkers),
    rng.normal(0.57, 0.02, args.nwalkers),
    rng.normal(0.11, 0.02, args.nwalkers),
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
    print(f"step {n:>7d}  tau_max = {tau.max():7.1f} ({tr.SHAPE_PARAM_NAMES[tau.argmax()]})  "
          f"N/tau = {n / tau.max():6.1f}  acc = {sampler.acceptance_fraction.mean():.2f}", flush=True)
    # many samples are wanted here because the reweighting to long-P tails uses them all
    if n >= args.min_steps and n > 50 * tau.max() and np.all(np.abs(old_tau - tau) / tau < 0.05):
        print("converged", flush=True)
        break
    old_tau = tau

THIN = 10
out = Path(__file__).resolve().parent / "results"
out.mkdir(exist_ok=True)
np.savez_compressed(
    out / "shape.npz",
    chain=sampler.get_chain(thin=THIN),
    log_prob=sampler.get_log_prob(thin=THIN),
    thin=THIN,
    acceptance=sampler.acceptance_fraction,
    tau_history=np.array(tau_history),
    param_names=tr.SHAPE_PARAM_NAMES,
)
