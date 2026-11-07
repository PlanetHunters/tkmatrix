import os
import tempfile
import unittest
from unittest import mock

import astropy.units as u
import numpy as np
import pandas as pd
import ellc


class TestEllcDependency(unittest.TestCase):
    def test_inject_model_uses_installed_ellc(self):
        """tkmatrix must use the ellc package from PyPI, not a bundled copy."""
        import tkmatrix
        import tkmatrix.inject_model as inject_model

        package_dir = os.path.dirname(os.path.abspath(ellc.__file__))
        tkmatrix_dir = os.path.dirname(os.path.abspath(tkmatrix.__file__))
        self.assertEqual("ellc", os.path.basename(package_dir))
        self.assertFalse(package_dir.startswith(tkmatrix_dir + os.sep),
                         f"ellc loaded from inside tkmatrix: {package_dir}")
        self.assertIs(inject_model.lc, ellc.lc)

    def test_lc_transit_depth(self):
        times = np.linspace(-0.1, 0.1, 201)
        flux = ellc.lc(times, radius_1=0.1, radius_2=0.01, sbratio=0.0, incl=90.0,
                       t_zero=0.0, period=3.1, ld_1="quad", ldc_1=[0.4, 0.2], verbose=0)
        self.assertAlmostEqual(1.0, flux[0], places=6)
        self.assertLess(flux[100], 1.0)

    def test_make_model_injects_transit_with_ellc(self):
        """InjectModel must build the transit with ellc, without falling back to batman."""
        from tkmatrix.inject_model import InjectModel

        time = np.linspace(0.0, 10.0, 7201)
        flux = np.ones_like(time)
        flux_err = np.full_like(time, 1e-4)
        with tempfile.TemporaryDirectory() as inject_dir, \
                mock.patch("tkmatrix.inject_model.batman.TransitModel",
                           side_effect=AssertionError("batman fallback used")):
            inject_model = InjectModel(inject_dir, time, flux, flux_err, 1.0 * u.R_sun, 1.0 * u.M_sun,
                                       t0=2.0, period=3.0, rplanet=2.0 * u.R_earth,
                                       exposure_time=120, ab=(0.4, 0.2))
            InjectModel.make_model(inject_model)
            files = os.listdir(inject_dir)
            self.assertEqual(1, len(files))
            injected = pd.read_csv(os.path.join(inject_dir, files[0]))
        self.assertEqual(len(time), len(injected))
        in_transit = np.abs(((injected["time"] - 2.0 + 1.5) % 3.0) - 1.5) < 0.02
        self.assertLess(injected["flux"][in_transit].mean(), 1.0 - 1e-4)
        self.assertAlmostEqual(1.0, injected["flux"][~in_transit].max(), places=6)


if __name__ == "__main__":
    unittest.main()
