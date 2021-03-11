import multiprocessing
import os
import shutil
import unittest
import astropy.units as u
from transitleastsquares import transitleastsquares

from tkmatrix.tkmatrix_class import TIRMA


class TestsTirma(unittest.TestCase):
    def test_inject_one(self):
        tirma = TIRMA("TIC 220513363", [1], ".")
        try:
            tirma.inject(1, 5, 5, 0.1, 3, 3, 0.1)
            self.assertEquals(1, len(os.listdir(tirma.build_inject_dir())))
            self.assertEquals([1], tirma.sectors)
            self.assertAlmostEqual(0.47, tirma.mass, 2)
            self.assertAlmostEqual(0.44, tirma.massmin, 2)
            self.assertAlmostEqual(0.5, tirma.massmax, 2)
            self.assertAlmostEqual(0.18, tirma.radius, 2)
            self.assertAlmostEqual(0.076, tirma.radiusmin, 3)
            self.assertAlmostEqual(0.284, tirma.radiusmax, 3)
            self.assertEquals("TIC 220513363", tirma.object_info.mission_id())
            self.assertEquals("TIC 220513363", tirma.id)
            self.assertAlmostEqual(0.47, tirma.mstar.value, 2)
            self.assertAlmostEqual(0.44, tirma.mstar_min.value, 2)
            self.assertAlmostEqual(0.5, tirma.mstar_max.value, 2)
            self.assertAlmostEqual(0.18, tirma.rstar.value, 2)
            self.assertAlmostEqual(0.076, tirma.rstar_min.value, 3)
            self.assertAlmostEqual(0.284, tirma.rstar_max.value, 3)
            self.assertEquals(".", tirma.dir)
            tirma.recovery(multiprocessing.cpu_count() - 1, 0)
            self.assertEquals(3, len(os.listdir(tirma.build_inject_dir())))
        finally:
            shutil.rmtree(tirma.build_inject_dir(), ignore_errors=True)

    def test_inject_four(self):
        tirma = TIRMA("TIC 220513363", [1], ".")
        try:
            tirma.inject(1, 5, 5.1, 0.1, 3, 3.1, 0.1)
            self.assertEquals(4, len(os.listdir(tirma.build_inject_dir())))
            tirma.recovery(multiprocessing.cpu_count() - 1, 0, None, 0.5)
            self.assertEquals(6, len(os.listdir(tirma.build_inject_dir())))
            with open(tirma.build_inject_dir() + "/a_tls_report.csv") as f:
                for i, l in enumerate(f):
                    pass
            self.assertEquals(5, i + 1)
        finally:
            shutil.rmtree(tirma.build_inject_dir(), ignore_errors=True)

    def test_inject_eight_multiphase(self):
        tirma = TIRMA("TIC 220513363", [1], ".")
        try:
            tirma.inject(2, 5, 5.1, 0.1, 3, 3.1, 0.1)
            self.assertEquals(8, len(os.listdir(tirma.build_inject_dir())))
            tirma.recovery(multiprocessing.cpu_count() - 1, 0, None, 0.5)
            self.assertEquals(10, len(os.listdir(tirma.build_inject_dir())))
            with open(tirma.build_inject_dir() + "/a_tls_report.csv") as f:
                for i, l in enumerate(f):
                    pass
            self.assertEquals(9, i + 1)
        finally:
            shutil.rmtree(tirma.build_inject_dir(), ignore_errors=True)

    def test_transit_mask(self):
        tirma = TIRMA("TIC 305048087", [2], ".")
        try:
            tirma.retrieve_object_data()
            power_args = {"transit_template": "default", "period_min": 1, "period_max": 10,
                          "n_transits_min": 2, "show_progress_bar": True, "use_threads": multiprocessing.cpu_count() - 1,
                          "R_star": tirma.radius, "R_star_min": tirma.radiusmin, "R_star_max": tirma.radiusmax,
                          "M_star": tirma.mass, "M_star_min": tirma.massmin, "M_star_max": tirma.massmax}
            intransit = tirma.transit_masks([{"P": 5.43, "T0": 2458355.249756 - 2457000.0, "D": 1.172585}],
                                            tirma.lc.time.value)
            model = transitleastsquares(tirma.lc.time.value, tirma.lc.flux.value)
            results = model.power(**power_args)
            self.assertAlmostEqual(5.43, results.period, 2)
            model = transitleastsquares(tirma.lc.time.value[~intransit], tirma.lc.flux.value[~intransit])
            results = model.power(**power_args)
            self.assertNotAlmostEqual(5.43, results.period, 2)
        finally:
            shutil.rmtree(tirma.build_inject_dir(), ignore_errors=True)
            

if __name__ == '__main__':
    unittest.main()
