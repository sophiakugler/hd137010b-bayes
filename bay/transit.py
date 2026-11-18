import dataclasses

import batman
import numpy as np

from .stats import ECC, IDX, LD

# cgs
G = 6.67430e-8
R_SUN = 6.957e10
R_EARTH = 6.371e8
M_SUN = 1.98892e33
AU = 1.495978707e13
DAY = 86400.0
K2_LONG_CADENCE = 29.4 / 1440.0   # exposure time [d]


def a_from_rho(per, rho_star): # a / R_star from Kepler's third law: (a/R)^3 = G rho P^2 / (3 pi)
    return (G * rho_star * (per * DAY) ** 2 / (3 * np.pi)) ** (1 / 3)


def inc_from_b(b, a, ecc=0.0, w=90.0): # Inclination [deg] for impact parameter b
    r_transit = a * (1 - ecc ** 2) / (1 + ecc * np.sin(np.radians(w)))
    return np.degrees(np.arccos(b / r_transit))


@dataclasses.dataclass
class TransitModel:
    t0: float                   # mid-transit time [d]
    per: float                  # orbital period [d]
    rp: float                   # R_planet / R_star
    inc: float                  # inclination [deg]
    a: float                    # a / R_star
    ecc: float
    w: float                    # argument of periastron [deg]
    u1: float                   # quadratic limb darkening
    u2: float
    R_star: float = 0.71        # [R_sun]
    L_star: float = 0.23        # [L_sun]
    rho_star: float = 2.85      # [g cm^-3]
    
    # habitable zone [AU], Kopparapu et al. 2013 as quoted by Venner et al. 2026
    hz_conservative_au: tuple = (0.48, 0.88)   # 1.00 .. 0.30 S_earth
    hz_optimistic_au: tuple = (0.38, 0.93)     # 1.60 .. 0.27 S_earth

    @classmethod
    def from_theta(cls, theta, circular=False): # Model for one sampler vector; R_star from stellar mass and density.

        per = np.exp(theta[IDX["ln_orbital_period"]])
        rho = 10.0 ** theta[IDX["log10_stellar_density"]]
        mass = 10.0 ** theta[IDX["log10_stellar_mass"]]
        h, q = (0.0, 0.0) if circular else theta[ECC]
        ecc, w = h ** 2 + q ** 2, np.degrees(np.arctan2(q, h))
        a = a_from_rho(per, rho)
        u1, u2 = theta[LD]
        return cls(t0=theta[IDX["mid_transit_time"]], per=per,
                   rp=np.exp(theta[IDX["ln_planet_to_star_radius_ratio"]]),
                   inc=inc_from_b(theta[IDX["impact_parameter"]], a, ecc, w),
                   a=a, ecc=ecc, w=w, u1=u1, u2=u2, rho_star=rho,
                   R_star=(3 * mass * M_SUN / (4 * np.pi * rho)) ** (1 / 3) / R_SUN)

    def to_batman(self):
        p = batman.TransitParams()
        p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w = (
            self.t0, self.per, self.rp, self.a, self.inc, self.ecc, self.w)
        p.limb_dark, p.u = "quadratic", [self.u1, self.u2]
        return p

    def light_curve(self, t, exp_time=K2_LONG_CADENCE, supersample_factor=15): # in F/F0 at time in days
        params = self.to_batman()
        kw = {} if exp_time is None else dict(exp_time=exp_time,
                                              supersample_factor=supersample_factor)
        return batman.TransitModel(params, np.asarray(t, dtype=float), **kw).light_curve(params)

    @property
    def b(self): # Impact parameter in R_star
        r_transit = self.a * (1 - self.ecc ** 2) / (1 + self.ecc * np.sin(np.radians(self.w)))
        return r_transit * np.cos(np.radians(self.inc))

    @property
    def R_planet_earth(self):
        return self.rp * self.R_star * R_SUN / R_EARTH


    ###############################
    #  period
    ###############################
    
    @property
    def period_years(self):
        return self.per / 365.25

    @property
    def rho_star_implied(self): # Stellar density implied by (a, per) in g cm^-3
        return 3 * np.pi * self.a ** 3 / (G * (self.per * DAY) ** 2)

    @property
    def per_from_density(self): # Period implied by (a, rho_star) in days
        return np.sqrt(3 * np.pi * self.a ** 3 / (G * self.rho_star)) / DAY


    ###############################
    # orbit and habitable zone
    ###############################
    
    @property
    def a_au(self):
        return self.a * self.R_star * R_SUN / AU

    @property
    def r_min_au(self): # Periapsis AU
        return self.a_au * (1 - self.ecc)

    @property
    def r_max_au(self): # Apoapsis AU
        return self.a_au * (1 + self.ecc)

    @property
    def r_min(self): #Periapsis R_Star
        return self.a * (1 - self.ecc)

    @property
    def r_max(self): # Apoapsis R_Star
        return self.a * (1 + self.ecc)

    @property
    def hz_conservative(self): # Edge Habitable Zone CONSERVATIVE in R_Star
        return tuple(r * AU / (self.R_star * R_SUN) for r in self.hz_conservative_au)

    @property
    def hz_optimistic(self): # Edge Habitable Zone OPTIMISTIC in R_Star
        return tuple(r * AU / (self.R_star * R_SUN) for r in self.hz_optimistic_au)

    @property
    def insolation(self): # Orbit-averaged stellar flux at the planet in S_earth
        return self.L_star / (self.a_au ** 2 * np.sqrt(1 - self.ecc ** 2))
