import multiprocessing
import os
import shutil
import unittest
import astropy.units as u
import pytest
from transitleastsquares import transitleastsquares

from tkmatrix.tkmatrix_class import MATRIX


class TestsMatrix(unittest.TestCase):
    def test_inject_one(self):
        matrix = MATRIX("TIC 220513363", [1], ".")
        inject_dir = None
        try:
            inject_dir = matrix.inject(1, 5, 5, 0.1, 3, 3, 0.1, 120)
            self.assertEquals(1, len(os.listdir(inject_dir)))
            self.assertEquals([1], matrix.sectors)
            self.assertAlmostEqual(0.47, matrix.mass, 2)
            self.assertAlmostEqual(0.44, matrix.massmin, 2)
            self.assertAlmostEqual(0.5, matrix.massmax, 2)
            self.assertAlmostEqual(0.18, matrix.radius, 2)
            self.assertAlmostEqual(0.076, matrix.radiusmin, 3)
            self.assertAlmostEqual(0.284, matrix.radiusmax, 3)
            self.assertEquals("TIC 220513363", matrix.object_info.mission_id())
            self.assertEquals("TIC 220513363", matrix.id)
            self.assertAlmostEqual(0.47, matrix.mstar.value, 2)
            self.assertAlmostEqual(0.44, matrix.mstar_min.value, 2)
            self.assertAlmostEqual(0.5, matrix.mstar_max.value, 2)
            self.assertAlmostEqual(0.18, matrix.rstar.value, 2)
            self.assertAlmostEqual(0.076, matrix.rstar_min.value, 3)
            self.assertAlmostEqual(0.284, matrix.rstar_max.value, 3)
            self.assertEquals(".", matrix.dir)
            matrix.recovery(multiprocessing.cpu_count() - 1, inject_dir, 5, 0)
            self.assertEquals(2, len(os.listdir(inject_dir)))
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)

    def test_inject_one_preserve(self):
        matrix = MATRIX("TIC 220513363", [1], ".", True)
        inject_dir = None
        try:
            inject_dir = matrix.inject(1, 5, 5, 0.1, 3, 3, 0.1, 120)
            matrix.recovery(multiprocessing.cpu_count() - 1, inject_dir, 5, 0)
            self.assertEquals(3, len(os.listdir(inject_dir)))
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)

    def test_inject_four(self):
        matrix = MATRIX("TIC 220513363", [1], ".")
        inject_dir = None
        try:
            inject_dir = matrix.inject(1, 5, 5.1, 0.1, 3, 3.1, 0.1, 120)
            self.assertEquals(4, len(os.listdir(inject_dir)))
            matrix.recovery(multiprocessing.cpu_count() - 1, inject_dir, 5, 0, None, 0.5)
            self.assertEquals(2, len(os.listdir(inject_dir)))
            with open(inject_dir + "/a_tls_report.csv") as f:
                for i, l in enumerate(f):
                    pass
            self.assertEquals(5, i + 1)
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)

    def test_inject_multiphase(self):
        matrix = MATRIX("TIC 220513363", [1], ".")
        inject_dir = None
        try:
            inject_dir = matrix.inject(2, 5, 5, 0.1, 3, 3, 0.1, 120)
            self.assertEquals(2, len(os.listdir(inject_dir)))
            matrix.recovery(multiprocessing.cpu_count() - 1, inject_dir, 5, 0, None, 0.5)
            self.assertEquals(2, len(os.listdir(inject_dir)))
            with open(inject_dir + "/a_tls_report.csv") as f:
                for i, l in enumerate(f):
                    pass
            self.assertEquals(2, i)
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)

    def test_transit_mask(self):
        matrix = MATRIX("TIC 305048087", [2], ".")
        matrix.retrieve_object_data()
        power_args = {"transit_template": "default", "period_min": 1, "period_max": 10,
                      "n_transits_min": 2, "show_progress_bar": True, "use_threads": multiprocessing.cpu_count() - 1,
                      "R_star": matrix.radius, "R_star_min": matrix.radiusmin, "R_star_max": matrix.radiusmax,
                      "M_star": matrix.mass, "M_star_min": matrix.massmin, "M_star_max": matrix.massmax}
        intransit = matrix.transit_masks([{"P": 5.43, "T0": 2458355.249756 - 2457000.0, "D": 1.172585}],
                                        matrix.lc.time.value)
        model = transitleastsquares(matrix.lc.time.value, matrix.lc.flux.value)
        results = model.power(**power_args)
        self.assertAlmostEqual(5.43, results.period, 2)
        model = transitleastsquares(matrix.lc.time.value[~intransit], matrix.lc.flux.value[~intransit])
        results = model.power(**power_args)
        self.assertNotAlmostEqual(5.43, results.period, 2)

    def test_inject_inputs(self):
        matrix = MATRIX("TIC 305048087", [2], ".")
        inject_dir = None
        with(pytest.raises(AssertionError)):
            inject_dir = matrix.inject(2, 5, 5.1, 0.1, 3, 3.1, 0.1, None)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir = matrix.inject(2, 5, 5.1, 0.1, 3, 3.1, 0.1, -2)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir = matrix.inject(2, 5, 5.1, 0.1, 3, 3.1, 0.1, "bla")
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)

    def test_inject_dir(self):
        inject_dir1 = None
        inject_dir2 = None
        matrix = MATRIX("TIC 305048087", [2], ".")
        try:
            inject_dir1 = matrix.inject(2, 5, 5, 0.1, 3, 3, 0.1, 120)
            self.assertTrue(os.path.isdir(inject_dir1))
            inject_dir2 = matrix.inject(2, 5, 5, 0.1, 3, 3, 0.1, 120)
            self.assertTrue(os.path.isdir(inject_dir2))
        finally:
            if inject_dir1 is not None:
                shutil.rmtree(inject_dir1)
            if inject_dir2 is not None:
                shutil.rmtree(inject_dir2)

if __name__ == '__main__':
    unittest.main()
