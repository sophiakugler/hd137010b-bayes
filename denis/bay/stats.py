import numpy as np

# Sampled vector. Saved chains store parameters by position: never reorder.
NAMES = (
    "mid_transit_time",                            # BKJD
    "ln_orbital_period",                           # ln(d)
    "impact_parameter",                            # R_star
    "ln_planet_to_star_radius_ratio",
    "limb_darkening_linear_coefficient",
    "limb_darkening_quadratic_coefficient",
    "log10_stellar_density",                       # log10(g cm^-3)
    "log10_stellar_mass",                          # log10(M_sun)
    "out_of_transit_flux",
    "ln_white_noise_amplitude",
    "sqrt_eccentricity_cos_periastron_argument",   # h = sqrt(e) cos w
    "sqrt_eccentricity_sin_periastron_argument",   # q = sqrt(e) sin w
)
IDX = {name: i for i, name in enumerate(NAMES)}
LD = [IDX["limb_darkening_linear_coefficient"], IDX["limb_darkening_quadratic_coefficient"]]
ECC = [IDX["sqrt_eccentricity_cos_periastron_argument"],
       IDX["sqrt_eccentricity_sin_periastron_argument"]]


def log_likelihood(flux, model_flux, sigma):
    r = (flux - model_flux) / sigma
    return -0.5 * np.sum(r ** 2) - flux.size * np.log(sigma * np.sqrt(2 * np.pi))


def physical_ld(u1, u2):
    return u1 > 0 and u1 + u2 < 1 and u1 + 2 * u2 > 0


def transit_factor(e, sin_w): # (1 + e sin w) / (1 - e^2): we only see planets that transit
    return (1 + e * sin_w) / (1 - e ** 2)


def e_logpdf(e, p):
    if p["dist"] == "rayleigh":
        return np.log(e) - 2 * np.log(p["scale"]) - 0.5 * (e / p["scale"]) ** 2
    if p["dist"] == "beta":
        return (p["a"] - 1) * np.log(e) + (p["b"] - 1) * np.log1p(-e)
    return 0.0          # uniform, circular: the (h, q) disc is already flat in e


def e_draw(p, rng):
    if p["dist"] == "rayleigh":
        return rng.rayleigh(p["scale"])
    if p["dist"] == "beta":
        return rng.beta(p["a"], p["b"])
    if p["dist"] in ("uniform", "circular"):
        return rng.uniform(0, p["e_max"])
    raise ValueError(f"eccentricity: unknown dist {p['dist']!r}")


def log_prior(theta, priors):
    lp = 0.0
    
    for key, p in priors.items():
        dist = p["dist"]
        
        if key == "limb_darkening":            # normal, cut to the physical triangle
            if not physical_ld(*theta[LD]):
                return -np.inf
            lp -= 0.5 * np.sum(((theta[LD] - p["mu"]) / p["sigma"]) ** 2)
            
        elif key == "eccentricity":            # circular: sampled, ignored by the model
            h, q = theta[ECC]
            e = h * h + q * q
            if not e < p["e_max"]:
                return -np.inf
            if dist != "circular":
                lp += e_logpdf(e, p)
                if p["transit_weighted"]:
                    lp += np.log(transit_factor(e, q / np.sqrt(e) if e > 0 else 0.0))
        else:
            x = theta[IDX[key]]
            if dist in ("uniform", "single_transit"):
                if not p["low"] <= x < p["high"]:
                    return -np.inf
                # single transit: p(ln P) ~ P^exponent
                lp += -np.log(p["high"] - p["low"]) if dist == "uniform" else p["exponent"] * x
                
            elif dist == "normal":
                lp -= 0.5 * ((x - p["mu"]) / p["sigma"]) ** 2
                
            elif dist == "impact":             # b ~ U(0, 1 + k)
                k = np.exp(theta[IDX["ln_planet_to_star_radius_ratio"]])
                if not 0 <= x <= 1 + k:
                    return -np.inf
                lp -= np.log1p(k)

    return lp


def _sample_one(priors, rng): # impact last: its range depends on the radius ratio
    theta = np.full(len(NAMES), np.nan)
    
    for key, p in sorted(priors.items(), key=lambda item: item[1]["dist"] == "impact"):
        dist = p["dist"]
        if key == "limb_darkening":
            u = rng.normal(p["mu"], p["sigma"])
            while not physical_ld(*u):
                u = rng.normal(p["mu"], p["sigma"])
            theta[LD] = u
            
        elif key == "eccentricity":
            weighted = dist != "circular" and p["transit_weighted"]
            while True:
                e, w = e_draw(p, rng), rng.uniform(0, 2 * np.pi)
                if e < p["e_max"] and (not weighted or rng.uniform(0, transit_factor(p["e_max"], 1))
                                       < transit_factor(e, np.sin(w))):
                    break
            theta[ECC] = np.sqrt(e) * np.cos(w), np.sqrt(e) * np.sin(w)
            
        elif dist == "uniform":
            theta[IDX[key]] = rng.uniform(p["low"], p["high"])
            
        elif dist == "normal":
            theta[IDX[key]] = rng.normal(p["mu"], p["sigma"])
            
        elif dist == "single_transit":         # inverse CDF of p(ln P) ~ P^exponent
            x = p["exponent"]
            lo, hi = np.exp(p["low"] * x), np.exp(p["high"] * x)
            theta[IDX[key]] = np.log(lo + rng.uniform() * (hi - lo)) / x
            
        elif dist == "impact":
            k = np.exp(theta[IDX["ln_planet_to_star_radius_ratio"]])
            theta[IDX[key]] = rng.uniform(0, 1 + k)
        
    if np.isnan(theta).any():
        raise ValueError(f"no prior for {[NAMES[i] for i in np.flatnonzero(np.isnan(theta))]}")
    return theta


def sample_prior(priors, n, rng=None):
    rng = np.random.default_rng(rng)
    return np.array([_sample_one(priors, rng) for _ in range(n)])
