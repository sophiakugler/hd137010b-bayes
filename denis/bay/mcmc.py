import batman
import emcee
import numpy as np

from .stats import IDX, NAMES, log_likelihood, log_prior, sample_prior
from .transit import K2_LONG_CADENCE, TransitModel, inc_from_b


class Posterior:
    def __init__(self, spec, data_path):
        data = np.loadtxt(data_path, skiprows=15)
        window = spec["window"]
        mask = np.abs(data[:, 0] - window["center"]) < window["half_width"]
        self.t, self.flux = data[mask, 0], data[mask, 2]      # detrended flux
        self.priors = spec["priors"]
        self.circular = self.priors["eccentricity"]["dist"] == "circular"
        # build batman once: the supersampled time grid is the slow part
        reference = TransitModel(window["center"], 350.0, 0.0136, inc_from_b(0.28, 267),
                                 267.0, 0.0, 90.0, 0.60, 0.12)
        self._batman = batman.TransitModel(reference.to_batman(), self.t,
                                           supersample_factor=spec["supersample"],
                                           exp_time=K2_LONG_CADENCE)

    def model_flux(self, theta):
        model = TransitModel.from_theta(theta, self.circular)
        return theta[IDX["out_of_transit_flux"]] * self._batman.light_curve(model.to_batman())

    def __call__(self, theta):
        lp = log_prior(theta, self.priors)
        if not np.isfinite(lp):
            return -np.inf
        sigma = np.exp(theta[IDX["ln_white_noise_amplitude"]])
        ll = log_likelihood(self.flux, self.model_flux(theta), sigma)
        return lp + ll if np.isfinite(ll) else -np.inf


def _sample(sampler, state, n, report, phase):
    """sampler.run_mcmc, reporting every quarter."""
    for i, _ in enumerate(sampler.sample(state, iterations=n), 1):
        if i % max(n // 4, 1) == 0 or i == n:
            report(f"{phase} {i}/{n}")


def run_chain(spec, data_path, report=print):
    """Burn in from the prior, restart every walker from the top `restart_keep`
    walkers, jittered by `restart_jitter` x their spread per parameter, sample."""
    rng = np.random.default_rng(spec["seed"])
    moves = [(getattr(emcee.moves, name)(), weight) for name, weight in spec["moves"].items()]
    sampler = emcee.EnsembleSampler(spec["nwalkers"], len(NAMES),Posterior(spec, data_path), moves=moves)
    sampler.random_state = np.random.RandomState(spec["seed"]).get_state()   # the moves' rng
    p0 = sample_prior(spec["priors"], spec["nwalkers"], rng)
    _sample(sampler, p0, spec["nburn"], report, "burn-in")
    last, lp = sampler.get_chain()[-1], sampler.get_log_prob()[-1]
    report(f"restart from top {spec['restart_keep']}, best log prob {lp.max():.1f}")
    top = last[np.argsort(lp)[-spec["restart_keep"]:]]
    p1 = top[rng.integers(len(top), size=spec["nwalkers"])]
    p1 += spec["restart_jitter"] * top.std(axis=0) * rng.normal(size=p1.shape)
    sampler.reset()                    # the burn-in does not stay in the chain
    _sample(sampler, p1, spec["nsteps"], report, "sampling")
    return sampler
