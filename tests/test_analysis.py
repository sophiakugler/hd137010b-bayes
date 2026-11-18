"""Numerical invariants for the transit model and prior measure."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np

spec=importlib.util.spec_from_file_location('analysis',Path(__file__).parents[1]/'src/analysis.py')
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)

class AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.x=np.array([np.log(355),.3,np.log(.01368),0,.462,.56,.10])
    def test_kepler_units(self):
        self.assertAlmostEqual(a.scaled_axis(355,2.90),267.7,delta=1.)
    def test_prior_log_period_jacobian(self):
        y=self.x.copy();y[0]+=np.log(2)
        self.assertAlmostEqual(a.log_prior(y,'circular')-a.log_prior(self.x,'circular'),-2/3*np.log(2))
    def test_circular_limit(self):
        t=np.linspace(a.TREF-.4,a.TREF+.4,100)
        model=a.Transit(t)
        np.testing.assert_allclose(model(self.x),model(np.r_[self.x,0.,0.]),atol=1e-12)
    def test_reused_batman_matches_fresh_model(self):
        t=np.linspace(a.TREF-.4,a.TREF+.4,100)
        model=a.Transit(t);model(self.x)
        for e,w in [(.6,-1.2),(.1,2.4)]:
            x=np.r_[self.x,np.sqrt(e)*np.cos(w),np.sqrt(e)*np.sin(w)]
            np.testing.assert_allclose(model(x),a.Transit(t)(x),atol=1e-12)
    def test_physical_flux_limits(self):
        d=a.derived(np.tile(self.x,(100,1)),np.random.default_rng(1))
        np.testing.assert_allclose(d['mean_flux_Earth'],d['flux_at_a_Earth'])
        np.testing.assert_allclose(d['periastron_flux_Earth'],d['apastron_flux_Earth'])
    def test_duration_central_limit(self):
        x=self.x.copy();x[1]=0
        expected=24*355/np.pi*np.arcsin((1+.01368)/a.scaled_axis(355,10**.462))
        self.assertAlmostEqual(a.duration(x),expected)

if __name__=='__main__':unittest.main()
