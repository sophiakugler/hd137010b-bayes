"""Nested sampling (dynesty) of the circular model in the original parametrisation.

    python run_dynesty.py --period-prior kipping [--e-prior beta_kipping13]
    python run_dynesty.py --shape         # shape-space model (flat ln v), cf. run_shape.py

Independent cross-check of the emcee runs and of the reweighted eccentric posteriors
(and gives the evidence ln Z).
Output: results/dynesty_<prior>[_<e-prior>].npz with equally weighted posterior samples.

Eccentric runs: p(e) p(w) is sampled through the prior transform, while the transit-
probability factor (1 + e sin w)/(1 - e^2) and the periastron cut enter the likelihood.
The stored ln Z is therefore relative to that unnormalised prior; the normalised evidence
is ln Z + ln C_P - ln C (C_P = integral of P^x dlnP, C from transit.log_prior_norm).
"""

import argparse
from pathlib import Path

import dynesty
import numpy as np
from scipy import stats

import transit as tr

parser = argparse.ArgumentParser()
parser.add_argument("--period-prior", default="kipping", choices=list(tr.PERIOD_PRIORS))
parser.add_argument("--e-prior", default="circular", choices=["circular", *tr.E_PRIORS])
parser.add_argument("--shape", action="store_true", help="sample the shape model instead")
parser.add_argument("--p-max", type=float, default=tr.P_MAX)
parser.add_argument("--nlive", type=int, default=1000)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

t, f, p_min = tr.load_data()
eccentric = args.e_prior != "circular"
Model = tr.EccentricTransitModel if eccentric else tr.CircularTransitModel
model = Model(t, f, p_min, period_prior=args.period_prior, p_max=args.p_max)
s = tr.PERIOD_PRIORS[args.period_prior]        # p(lnP) ~ P^s on [p_min, p_max]
ndim = len(tr.PARAM_NAMES) + 2 * eccentric
shape_model = tr.ShapeTransitModel(t, f)


def shape_prior_transform(u):
    x = np.empty_like(u)
    x[0] = tr.T0_RANGE[0] + u[0] * (tr.T0_RANGE[1] - tr.T0_RANGE[0])
    x[1] = tr.LN_V_RANGE[0] + u[1] * (tr.LN_V_RANGE[1] - tr.LN_V_RANGE[0])
    x[3] = np.log(1e-3) + u[3] * np.log(100)
    x[2] = u[2] * (1 + np.exp(x[3]))
    x[4] = stats.norm.ppf(u[4], *tr.U1_PRIOR)
    x[5] = stats.norm.ppf(u[5], *tr.U2_PRIOR)
    x[6] = 0.99 + 0.02 * u[6]
    x[7] = -15 + 10 * u[7]
    return x


def prior_transform(u):
    x = np.empty_like(u)
    x[0] = tr.T0_RANGE[0] + u[0] * (tr.T0_RANGE[1] - tr.T0_RANGE[0])
    lo, hi = p_min**s, args.p_max**s
    x[1] = np.log(lo + u[1] * (hi - lo)) / s if s != 0 else np.log(p_min) + u[1] * np.log(args.p_max / p_min)
    x[3] = np.log(1e-3) + u[3] * np.log(100)               # ln k, k in (1e-3, 0.1)
    x[2] = u[2] * (1 + np.exp(x[3]))                        # b ~ U(0, 1+k)
    x[4] = stats.norm.ppf(u[4], *tr.U1_PRIOR)
    x[5] = stats.norm.ppf(u[5], *tr.U2_PRIOR)
    x[6] = stats.norm.ppf(u[6], *tr.LOG_RHO_PRIOR)
    x[7] = stats.norm.ppf(u[7], *tr.LOG_M_PRIOR)
    x[8] = 0.99 + 0.02 * u[8]
    x[9] = -15 + 10 * u[9]
    if eccentric:
        x[10] = tr.E_PRIORS[args.e_prior].ppf(u[10])
        x[11] = 2 * np.pi * u[11]
    return x


def log_like(x):
    # the physical limb-darkening bounds truncate the prior; impose them here
    u1, u2 = x[4], x[5]
    if not (u1 > 0 and u1 + u2 < 1 and u1 + 2 * u2 > 0):
        return -1e300
    if args.shape:
        return shape_model.log_likelihood(x)
    if not eccentric:
        return model.log_likelihood(x)
    P, e, w = np.exp(x[1]), x[10], x[11]
    lw = tr.log_orbit_weight(P, e, w, x[6], args.period_prior, p_min, args.p_max) - s * x[1]
    if not np.isfinite(lw):
        return -1e300
    try:
        return model.log_likelihood(x) + lw
    except RuntimeError:  # batman's Kepler solver fails for e -> 1: use the shape model there
        v = tr.a_over_rstar(P, x[6]) / P * tr.g_factor(e, w)
        return shape_model.log_likelihood(np.r_[x[0], np.log(v), x[2:6], x[8:10]]) + lw


if args.shape:
    prior_transform, ndim = shape_prior_transform, len(tr.SHAPE_PARAM_NAMES)
sampler = dynesty.NestedSampler(log_like, prior_transform, ndim, nlive=args.nlive,
                                sample="rslice", rstate=np.random.default_rng(args.seed))
sampler.run_nested(dlogz=0.01, print_progress=False)
res = sampler.results
weights = np.exp(res.logwt - res.logz[-1])
samples = dynesty.utils.resample_equal(res.samples, weights / weights.sum(), np.random.default_rng(args.seed))
print(f"ln Z = {res.logz[-1]:.3f} ± {res.logzerr[-1]:.3f}, n_eff = {res.eff:.1f}%, "
      f"n_samples = {len(samples)}, n_like = {np.sum(res.ncall)}")

out = Path(__file__).resolve().parent / "results"
np.savez_compressed(
    out / ("dynesty_shape.npz" if args.shape else
           f"dynesty_{args.period_prior}" + (f"_{args.e_prior}" if eccentric else "")
           + ("" if args.p_max == tr.P_MAX else f"_pmax{args.p_max:.0f}") + ".npz"),
    samples=samples, logz=res.logz[-1], logzerr=res.logzerr[-1],
    param_names=tr.SHAPE_PARAM_NAMES if args.shape else tr.PARAM_NAMES + ["e", "w"] * eccentric,
)
