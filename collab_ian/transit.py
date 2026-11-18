"""Single-transit model for HD 137010 b (Phase 1: circular-orbit baseline).

Reproduces the setup of Venner et al. (2026), Sec. 3.4:
batman quadratic-limb-darkening transit, supersampled to the 29.4-min K2
long cadence, fitted to the detrended K2 C15 photometry.

Parameter vector (all times in BKJD = BJD - 2454833):

    T0        transit mid-time [d]                     uniform
    lnP       ln orbital period [d]                    see `PERIOD_PRIORS`
    b         impact parameter                         uniform [0, 1+k)
    lnk       ln(Rp/R*)                                log-uniform
    u1, u2    quadratic limb darkening                 Gaussian (Claret 2018)
    log_rho   log10 stellar density [g cm^-3]          Gaussian (isochrones)
    log_M     log10 stellar mass [Msun]                Gaussian (isochrones)
    f0        out-of-transit flux level                uniform
    ln_sigma  ln white-noise level per point           log-uniform

log_M does not affect the light curve; it is only carried along (as in the
paper) to convert P into an absolute semi-major axis.
"""

import warnings
from pathlib import Path

import batman
import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA_FILE = HERE / "data" / "apjladf06ft3_mrt.txt"

# --- constants (cgs) ---
G = 6.674e-8
MSUN = 1.989e33
RSUN = 6.957e10
REARTH = 6.371e8
AU = 1.496e13
DAY = 86400.0

# --- K2 long cadence ---
EXP_TIME = 29.4 / 60 / 24  # d
SUPERSAMPLE = 10

# --- stellar / limb-darkening priors from Venner et al. (2026) ---
LOG_RHO_PRIOR = (0.462, 0.041)   # log10 g cm^-3
LOG_M_PRIOR = (-0.139, 0.010)    # log10 Msun
U1_PRIOR = (0.60, 0.07)
U2_PRIOR = (0.12, 0.07)
LUM = (0.232, 0.022)             # Lsun (symmetrised Table 1 value)
TEFF = (4770.0, 90.0)            # K

T0_GUESS = 2458055.0969 - 2454833.0
T0_RANGE = (T0_GUESS - 0.2, T0_GUESS + 0.2)
P_MAX = 2000.0  # d, upper cut as in the paper

PARAM_NAMES = ["T0", "lnP", "b", "lnk", "u1", "u2", "log_rho", "log_M", "f0", "ln_sigma"]

# Period priors. Kipping (2018): for an intrinsic period distribution p(P) ~ P^alpha,
# the prior for a planet seen to transit exactly once (P > baseline) is
#     p(P | single transit) ~ P^alpha * P^(-2/3) [transit prob.] * P^(-1) [window]
#                           = P^(alpha - 5/3).
# Since we sample lnP, the density used is p(lnP) = p(P) P ~ P^(alpha - 2/3).
# Value = exponent x in p(lnP) ~ P^x.
PERIOD_PRIORS = {
    "kipping": -2 / 3,                # alpha = 0    (uniform in P) -> p(P) ~ P^-5/3
    "uniform_a_kipping": -1.0,        # alpha = -1/3 (uniform in a ~ P^(2/3)) -> p(P) ~ P^-2
    "loguniform_kipping": -5 / 3,     # alpha = -1   (log-uniform)  -> p(P) ~ P^-8/3  [paper, per authors' samples]
    "loguniform": 0.0,                # log-uniform, no single-transit correction
}

# Eccentricity priors (omega uniform in all cases). The geometric transit probability
# additionally scales as (1 + e sin w) / (1 - e^2) (Barnes 2007; Kipping 2014), which
# must multiply any of these when conditioning on "the planet transits".
E_PRIORS = {
    "uniform": stats.uniform(0, 1),
    "beta_kipping13": stats.beta(0.867, 3.03),            # Kipping (2013), RV planets
    "rayleigh_vaneylen19": stats.rayleigh(scale=0.049),   # Van Eylen+ (2019), multis
    # Kipping+ (2025): mean e of temperate Earth-sized planets <~ 0.15. Their posterior is not
    # tabulated, so we take the most permissive member: a Rayleigh with mean e = 0.15.
    "rayleigh_mean015": stats.rayleigh(scale=0.15 / np.sqrt(np.pi / 2)),
}

# Kopparapu+ (2013) habitable-zone limits for HD 137010 as quoted by Venner+ (2026), Sec. 4.2,
# in units of Earth's insolation: (outer, inner) edges.
HZ_CONSERVATIVE = (0.30, 1.00)
HZ_OPTIMISTIC = (0.27, 1.60)


def transit_prob_factor(e, w):
    """Eccentricity part of the geometric transit probability, (1 + e sin w)/(1 - e^2)."""
    return (1 + e * np.sin(w)) / (1 - e**2)


def g_factor(e, w):
    """Transit-speed boost relative to a circular orbit: (1 + e sin w)/sqrt(1 - e^2).

    For fixed rho*, T14 ~ P^(1/3) / g, so a duration fitted as circular maps to
    P_true = P_circ * g^3 (photo-eccentric effect, Dawson & Johnson 2012).
    """
    return (1 + e * np.sin(w)) / np.sqrt(1 - e**2)


def load_data(window=1.5):
    """Detrended K2 flux within +-window days of the transit."""
    t, _, f = np.loadtxt(DATA_FILE, skiprows=15, unpack=True)
    # Shortest period that does not put a second transit inside the K2
    # baseline (ignoring data gaps). Irrelevant in practice: the ~10 h
    # duration already forces P >~ 200 d for a circular orbit.
    p_min = max(t.max() - T0_GUESS, T0_GUESS - t.min())
    m = np.abs(t - T0_GUESS) < window
    return t[m], f[m], p_min


def a_over_rstar(P, log_rho):
    """Kepler's third law: a/R* = (G rho P^2 / 3 pi)^(1/3)."""
    rho = 10.0**log_rho
    return (G * rho * (P * DAY) ** 2 / (3 * np.pi)) ** (1 / 3)


class CircularTransitModel:
    def __init__(self, t, f, p_min, period_prior="loguniform_kipping", p_max=P_MAX):
        if period_prior not in PERIOD_PRIORS:
            raise ValueError(f"period_prior must be one of {list(PERIOD_PRIORS)}")
        self.t, self.f = t, f
        self.p_min, self.p_max = p_min, p_max
        self.period_prior = period_prior

        self.params = batman.TransitParams()
        self.params.t0 = T0_GUESS
        self.params.per = 355.0
        self.params.rp = 0.0137
        self.params.a = 270.0
        self.params.inc = 90.0
        self.params.ecc = 0.0
        self.params.w = 90.0
        self.params.limb_dark = "quadratic"
        self.params.u = [0.6, 0.12]
        self.batman = batman.TransitModel(
            self.params, t, supersample_factor=SUPERSAMPLE, exp_time=EXP_TIME
        )

    def light_curve(self, theta, t=None):
        T0, lnP, b, lnk, u1, u2, log_rho, _, f0, _ = theta
        P, k = np.exp(lnP), np.exp(lnk)
        aR = a_over_rstar(P, log_rho)
        p = self.params
        p.t0, p.per, p.rp, p.a = T0, P, k, aR
        p.inc = np.degrees(np.arccos(b / aR))
        p.u = [u1, u2]
        model = self.batman
        if t is not None:
            model = batman.TransitModel(p, t, supersample_factor=SUPERSAMPLE, exp_time=EXP_TIME)
        return f0 * model.light_curve(p)

    def log_prior(self, theta):
        T0, lnP, b, lnk, u1, u2, log_rho, log_M, f0, ln_sigma = theta
        P, k = np.exp(lnP), np.exp(lnk)
        if not (T0_RANGE[0] < T0 < T0_RANGE[1]):
            return -np.inf
        if not (self.p_min < P < self.p_max):
            return -np.inf
        if not (1e-3 < k < 0.1 and 0 <= b < 1 + k):
            return -np.inf
        # physical limb darkening (Kipping 2013): u1 > 0, u1 + u2 < 1, u1 + 2 u2 > 0
        if not (u1 > 0 and u1 + u2 < 1 and u1 + 2 * u2 > 0):
            return -np.inf
        if not (0.99 < f0 < 1.01 and -15 < ln_sigma < -5):
            return -np.inf

        lp = PERIOD_PRIORS[self.period_prior] * lnP
        lp += -0.5 * ((u1 - U1_PRIOR[0]) / U1_PRIOR[1]) ** 2
        lp += -0.5 * ((u2 - U2_PRIOR[0]) / U2_PRIOR[1]) ** 2
        lp += -0.5 * ((log_rho - LOG_RHO_PRIOR[0]) / LOG_RHO_PRIOR[1]) ** 2
        lp += -0.5 * ((log_M - LOG_M_PRIOR[0]) / LOG_M_PRIOR[1]) ** 2
        return lp

    def log_likelihood(self, theta):
        sigma = np.exp(theta[-1])
        r = self.f - self.light_curve(theta)
        return -0.5 * np.sum((r / sigma) ** 2) - len(r) * np.log(sigma)

    def __call__(self, theta):
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(theta)


class EccentricTransitModel(CircularTransitModel):
    """Physical model with free eccentricity: theta = PARAM_NAMES + ["e", "w"], w in rad
    (batman convention: transit at true anomaly pi/2 - w). Used for the direct dynesty
    cross-check of the reweighted posteriors; the priors live in run_dynesty.py."""

    def light_curve(self, theta, t=None):
        T0, lnP, b, lnk, u1, u2, log_rho, _, f0, _, e, w = theta
        P = np.exp(lnP)
        aR = a_over_rstar(P, log_rho)
        p = self.params
        p.t0, p.per, p.rp, p.a, p.ecc, p.w = T0, P, np.exp(lnk), aR, e, np.degrees(w)
        p.inc = np.degrees(np.arccos(b * (1 + e * np.sin(w)) / (aR * (1 - e**2))))
        p.u = [u1, u2]
        model = self.batman
        if t is not None:
            model = batman.TransitModel(p, t, supersample_factor=SUPERSAMPLE, exp_time=EXP_TIME)
        return f0 * model.light_curve(p)

    def log_likelihood(self, theta):
        sigma = np.exp(theta[9])
        r = self.f - self.light_curve(theta)
        return -0.5 * np.sum((r / sigma) ** 2) - len(r) * np.log(sigma)


# --- shape-space parametrisation -------------------------------------------------
#
# For a/R* >> 1 the transit chord is a straight line, so the light curve depends on
# the orbit only through the chord speed v = a/(R* P) * g(e, w)  [1/d] and b.
# We therefore sample the light-curve "shape" with a flat prior in ln v, and map to
# the physical parameters afterwards:
#     P = G rho* DAY^2 g^3 / (3 pi v^3)          (Kepler's 3rd law, ρ* in cgs)
# Since d lnP / d ln v = -3 is constant, the posterior under any prior p(P, rho*, e, w)
# is obtained by importance weights  w = p(lnP) p(rho*) p(e, w | transit)
# with rho*, e, w drawn from their priors (see `shape_to_physical`).

SHAPE_PARAM_NAMES = ["T0", "ln_v", "b", "lnk", "u1", "u2", "f0", "ln_sigma"]
P_REF = 365.0      # d, nominal period for batman; irrelevant for a/R* >> 1
LN_V_RANGE = (np.log(0.02), np.log(20.0))


def period_from_v(v, log_rho, g=1.0):
    """Orbital period [d] for chord speed v = a/(R* P) * g [1/d]."""
    return G * 10.0**log_rho * DAY**2 * g**3 / (3 * np.pi * v**3)


class ShapeTransitModel(CircularTransitModel):
    """Posterior over light-curve shape parameters (SHAPE_PARAM_NAMES), flat in ln v."""

    def __init__(self, t, f):
        super().__init__(t, f, p_min=0.0)

    def light_curve(self, theta, t=None):
        T0, ln_v, b, lnk, u1, u2, f0, _ = theta
        aR = np.exp(ln_v) * P_REF
        p = self.params
        p.t0, p.per, p.rp, p.a = T0, P_REF, np.exp(lnk), aR
        p.inc = np.degrees(np.arccos(b / aR))
        p.u = [u1, u2]
        model = self.batman
        if t is not None:
            model = batman.TransitModel(p, t, supersample_factor=SUPERSAMPLE, exp_time=EXP_TIME)
        return f0 * model.light_curve(p)

    def log_prior(self, theta):
        T0, ln_v, b, lnk, u1, u2, f0, ln_sigma = theta
        k = np.exp(lnk)
        if not (T0_RANGE[0] < T0 < T0_RANGE[1] and LN_V_RANGE[0] < ln_v < LN_V_RANGE[1]):
            return -np.inf
        if not (1e-3 < k < 0.1 and 0 <= b < 1 + k):
            return -np.inf
        if not (u1 > 0 and u1 + u2 < 1 and u1 + 2 * u2 > 0):
            return -np.inf
        if not (0.99 < f0 < 1.01 and -15 < ln_sigma < -5):
            return -np.inf
        lp = -0.5 * ((u1 - U1_PRIOR[0]) / U1_PRIOR[1]) ** 2
        lp += -0.5 * ((u2 - U2_PRIOR[0]) / U2_PRIOR[1]) ** 2
        return lp


def draw_e_w(e_prior, n, rng):
    """Draw (e, w) from an eccentricity prior (or e = 0) and a uniform w."""
    om = rng.uniform(0, 2 * np.pi, n)
    e = np.zeros(n) if e_prior == "circular" else E_PRIORS[e_prior].rvs(n, random_state=rng)
    return e, om


def log_orbit_weight(P, e, om, log_rho, period_prior, p_min, p_max):
    """Unnormalised log prior of (lnP, e, w) given 'transits once', relative to p(e) p(w):
    P^x (single-transit period prior) x (1 + e sin w)/(1 - e^2) (transit probability).
    -inf outside the allowed region: period range, and periastron a(1 - e) > 2 R*
    (no collision with the star; ~ Roche limit of a rocky planet for this star's density)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        aR = a_over_rstar(P, log_rho)
        ok = (e < 1) & (P > p_min) & (P < p_max) & (aR * (1 - e) > 2)
        logw = PERIOD_PRIORS[period_prior] * np.log(P) + np.log(transit_prob_factor(e, om))
    return np.where(ok, logw, -np.inf)


def log_prior_norm(period_prior, e_prior, p_min, p_max, n=4_000_000, rng=None):
    """ln C, C = normalisation of the orbit prior: integral of exp(log_orbit_weight)
    over lnP in [p_min, p_max] and p(e) p(w) p(rho*) (Monte Carlo)."""
    rng = np.random.default_rng(rng)
    lnP = rng.uniform(np.log(p_min), np.log(p_max), n)
    e, om = draw_e_w(e_prior, n, rng)
    lw = log_orbit_weight(np.exp(lnP), e, om, rng.normal(*LOG_RHO_PRIOR, n), period_prior, p_min, p_max)
    return np.log(np.log(p_max / p_min)) + _log_mean_exp(lw)


def _log_mean_exp(x):
    m = x[np.isfinite(x)].max()
    return m + np.log(np.mean(np.exp(x - m)))


def shape_to_physical(shape, period_prior="kipping", e_prior="circular", p_min=64.2,
                      p_max=P_MAX, n_draws=10, rng=None, return_lnz=False, data=None):
    """Turn shape-space samples into weighted physical samples under a given prior.

    Each shape sample is paired with `n_draws` draws of (log_rho, log_M, e, w) from
    their priors. Returns (theta, e, w, weights) where theta rows follow PARAM_NAMES
    and the weights are normalised importance weights.

    With return_lnz=True also returns ln(Z_model / Z_shape), the evidence of this prior
    relative to the flat-in-ln-v shape model (needs only the prior normalisation):
        Z_model / Z_shape = E_post[pi_model(ln v)] * |LN_V_RANGE|,
        pi_model(ln v)    = 3 / C * E_{rho, e, w}[exp(log_orbit_weight)].

    The shape model is exact only for circular orbits: on an eccentric orbit the planet's
    distance (hence b and speed) drifts during the transit. Passing data=(t, f) multiplies
    the weights by L_exact / L_shape, with L_exact the full eccentric batman likelihood,
    which makes the reweighting exact.
    """
    rng = np.random.default_rng(rng)
    s = np.repeat(shape, n_draws, axis=0)
    n = len(s)
    log_rho = rng.normal(*LOG_RHO_PRIOR, n)
    log_M = rng.normal(*LOG_M_PRIOR, n)
    e, om = draw_e_w(e_prior, n, rng)
    with np.errstate(invalid="ignore", divide="ignore"):
        P = period_from_v(np.exp(s[:, 1]), log_rho, g_factor(e, om))
    logw = log_orbit_weight(P, e, om, log_rho, period_prior, p_min, p_max)
    if data is not None and e_prior != "circular":
        shp, ecc = ShapeTransitModel(*data), EccentricTransitModel(*data, p_min)
        lnl_shape = np.repeat([shp.log_likelihood(x) for x in shape], n_draws)
        failed = []
        for i in np.flatnonzero(np.isfinite(logw)):
            x = np.r_[s[i, 0], np.log(P[i]), s[i, 2:6], log_rho[i], log_M[i], s[i, 6:8], e[i], om[i]]
            try:
                logw[i] += ecc.log_likelihood(x) - lnl_shape[i]
            except RuntimeError:  # batman's Kepler solver fails for e -> 1: keep shape approximation
                failed.append(i)
        if failed:
            lw = np.where(np.isfinite(logw), logw, -np.inf)
            frac = np.exp(_log_mean_exp(lw[failed]) - _log_mean_exp(lw)) * len(failed) / len(lw)
            warnings.warn(f"{len(failed)} draws (max e = {e[failed].max():.4f}, weight {frac:.1e}) "
                          "kept the shape approximation")
    w = np.exp(logw - logw[np.isfinite(logw)].max())
    w /= w.sum()
    theta = np.column_stack([s[:, 0], np.log(P), s[:, 2], s[:, 3], s[:, 4], s[:, 5],
                             log_rho, log_M, s[:, 6], s[:, 7]])
    if not return_lnz:
        return theta, e, om, w
    lnz = (_log_mean_exp(logw) + np.log(3 * (LN_V_RANGE[1] - LN_V_RANGE[0]))
           - log_prior_norm(period_prior, e_prior, p_min, p_max, rng=rng))
    return theta, e, om, w, lnz


def derived(samples, rng=None, e=None):
    """Physical quantities from posterior samples (rows = PARAM_NAMES).

    T14 and v_tr are the circular-orbit values. If e is given, insolation and Teq are
    orbit-averaged: <S> = L / (a^2 sqrt(1 - e^2)).
    """
    rng = np.random.default_rng(rng)
    T0, lnP, b, lnk, u1, u2, log_rho, log_M, f0, ln_sigma = samples.T
    ecc_fac = 1.0 if e is None else 1 / np.sqrt(1 - e**2)
    P, k = np.exp(lnP), np.exp(lnk)
    aR = a_over_rstar(P, log_rho)
    inc = np.arccos(b / aR)
    M = 10**log_M * MSUN
    rho = 10**log_rho
    R = (3 * M / (4 * np.pi * rho)) ** (1 / 3)
    a = aR * R
    arg = np.sqrt(np.clip((1 + k) ** 2 - b**2, 0, None)) / (aR * np.sin(inc))
    T14 = P / np.pi * np.arcsin(np.clip(arg, 0, 1))  # circular-orbit total duration [d]
    # stellar L and Teff are not fitted: draw from their (Gaussian) uncertainties
    L = rng.normal(*LUM, size=len(P))
    teff = rng.normal(*TEFF, size=len(P))
    a_au = a / AU
    # observed depth at mid-transit for quadratic limb darkening (small-planet limit)
    mu = np.sqrt(np.clip(1 - b**2, 0, None))
    depth = k**2 * (1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2) / (1 - u1 / 3 - u2 / 6)
    return {
        "P": P,
        "k": k,
        "b": b,
        "depth_ppm": depth * 1e6,
        "T14_hr": T14 * 24,
        "v_tr": 2 * np.sqrt(np.clip((1 + k) ** 2 - b**2, 0, None)) / T14,  # R*/d, Osborn+16
        "a_over_R": aR,
        "inc_deg": np.degrees(inc),
        "R_star": R / RSUN,
        "Rp_earth": k * R / REARTH,
        "a_au": a_au,
        "insolation": L / a_au**2 * ecc_fac,
        "Teq_0.3": teff * np.sqrt(R / (2 * a)) * (1 - 0.3) ** 0.25 * ecc_fac**0.25,
    }
