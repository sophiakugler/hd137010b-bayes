# Project audit and reproducible extension — 4 November 2026

## Starting state

Branch `sophia`, commit `6f630cf`; working tree clean when inspected. No branch switch, commit, push, reset, or deletion was performed. Existing notebooks and posterior arrays were retained.

* `01`: data exploration exists, with stored outputs.
* `02`: exposure-integrated illustrative circular transit model exists.
* `03`: five-parameter circular MCMC exists, but the on-disk source uses flat stellar density and period priors, fixed limb darkening, and has no saved chain export or Venner comparison. It differs from the more developed workflow described in the earlier chat.
* `04`: empty notebook shell; `05`: zero-byte file.
* `06`: exploratory TESS Sector 91 downloads, plots, and hand-picked transit predictions. This is not an injection/recovery analysis or a non-detection likelihood. `lightkurve` is imported but absent from the declared project dependencies.
* Five posterior arrays exist. The two broad-eccentricity arrays are identical. Their filenames, dimensions, and approximate distributions suggest circular 7-column and eccentric 9-column parameter sets, but generating code, exact priors, chain ordering, burn-in, seeds, and diagnostics are missing. They are not used as verified inference in this extension.
* README links an Overleaf report. No local `.tex`, `.bib`, or report PDF was present. Overleaf contents were not inspected or overwritten.

## Added workflow

`src/analysis.py` fits the same published detrended K2 photometry with batman, integrating 29.4-minute exposures at 10 substeps. All three runs share their likelihood and nuisance priors. Circular, truncated exponential eccentricity (scale 0.06), and uniform eccentricity on [0,0.8) differ only in eccentricity assumptions. These are explicitly priors conditional on studying an observed transiting candidate, not an inference from an intrinsic planet population.

The physical period density is P^(-5/3) on [90,2000] days. The lower bound exceeds the K2 observing baseline; shorter-period eccentric solutions are excluded by design, not ruled out by this local likelihood. The upper bound follows the discovery analysis but needs a sensitivity extension for broad eccentricity. Density is Gaussian in log10(rho), mean 0.462 and SD 0.041. Radius ratio has a bounded log-uniform prior; impact parameter is uniform on [0,1); this excludes grazing transits. Limb-darkening coefficients have Gaussian priors with physical constraints. Midtime is represented relative to the published epoch for numerical stability. A constant flux baseline is marginalized. The white-noise scale is estimated from local out-of-transit points and held fixed.

The seven sampled coordinates are ln(P), b, ln(Rp/Rstar), t0 offset, log10(rho), u1, u2. Eccentric runs append sqrt(e)cos(omega), sqrt(e)sin(omega); their constant Jacobian gives a uniform density in e and omega before the exponential weight. Log-period sampling includes the required P Jacobian; do not add a second 1/P factor. The eccentric impact parameter relation is used when converting to inclination.

Results include complete chains and log posterior values, unthinned retained samples, random state checkpoints, seeds, input hashes, package versions, autocorrelation estimates, approximate ESS, descriptive split Rhat, traces, circular ECDF validation, prior-predictive duration distributions, posterior model and residual checks, period/eccentricity degeneracy, actual representative posterior transit curves, derived physical samples, and machine-readable quantiles. Checkpoints preserve all completed blocks; automatic continuation from checkpoints is not implemented.

## Scientific limits and completion criteria

1. Inspect `*_diagnostics.json`: short chains that fail the convergence screen are exploratory. Ensemble walkers are dependent, so the descriptive split Rhat is not a formal independent-chain diagnostic. Repeat runs from dispersed starts and extend until tail quantiles and autocorrelation estimates stabilize. Never infer convergence from a plausible-looking histogram alone.
2. Circular validation uses the supplied `P` column directly and compares full ECDFs. KS distance and Wasserstein distance are descriptive, without iid p-values. Different noise, nuisance, baseline, and selection choices prevent claiming exact replication. Inspect differences as well as matching medians.
3. Selection: P^(-5/3) is a long-period single-transit baseline convention. A full intrinsic-population analysis requires geometric transit conditioning, the exact observing window and detection efficiency, and explicit period population assumptions. Do not append a geometric factor to these conditional priors without redefining the model. Repeat the eccentric experiment with transit-conditioned population priors and alternative period limits before claiming population-level robustness.
4. Prior-predictive duration uses the usual constant-velocity eccentric-duration approximation; high-eccentricity cases are illustrative. The photometric likelihood itself uses batman. The displayed duration histogram is restricted to 0–40 hours.
5. Derived stellar mass is drawn in log10(M), mean -0.139, SD 0.010, and stellar radius is derived consistently from mass and sampled density. Teff is drawn independently (4770 +/- 90 K). The joint stellar isochrone covariance is unavailable, so derived radius/luminosity uncertainties are approximate. These are propagated Monte Carlo uncertainties, not additional transit constraints.
6. Flux at a is L/a^2; orbital mean flux is L/[a^2 sqrt(1-e^2)]. Periastron/apastron flux uses a(1-/+e). These do not establish surface temperature, liquid water, or habitability. No new habitable-zone probability is claimed.
7. TESS remains exploratory. A non-detection must not eliminate posterior samples until actual cadence coverage and injection/recovery efficiency have been evaluated. The original TESS notebook is retained.
8. `data/posteriors.csv` was already tracked at the starting commit. The extension does not publish it. Confirm permission before sharing a repository/archive containing the supplied samples.

## Reproduce

From the repository root, use the existing environment:

```
.venv\Scripts\python.exe tests\test_analysis.py
.venv\Scripts\python.exe src\analysis.py --steps 20000 --output results\new_run
```

Use a new output directory for each run (existing output directories are refused). Default is 48 walkers and seed 137010. The full 3-model run can take several minutes. Results are not automatically committed. Use `git diff` and review before staging.

## Source

Venner et al. (2026), https://doi.org/10.3847/2041-8213/adf06f ; open methods: https://arxiv.org/html/2601.19870v1 . This extension is a controlled course-project prior experiment, not the authors' full discovery pipeline.

## Executed analysis status

The user requested code and analysis only; optional report/presentation scaffolds were removed. The extra independent circular run was interrupted; its checkpoint is retained and excluded from final results.

| Model | Steps per walker | Convergence screen | Period q16 / q50 / q84 [d] |
|---|---:|---|---|
| circular | 60000 | True | 302.9 / 379.0 / 703.6 |
| lowe | 60000 | True | 288.1 / 384.6 / 703.3 |
| uniform | 70000 | True | 151.7 / 332.1 / 793.9 |

Passing the screen is not independent-ensemble validation. Independent eccentricity runs, broader period-support checks, and noise/selection sensitivity remain scientific follow-up work.

Circular ECDF KS distance against supplied samples: 0.100. The period tail differs; this is not an exact replication.
