import multiprocessing
import os
import shutil
import unittest
import pytest
from lcbuilder.star.starinfo import StarInfo
from foldedleastsquares import transitleastsquares

from tkmatrix.tkmatrix_class import MATRIX


class TestsMatrix(unittest.TestCase):
    def test_inject_one(self):
        target = "TIC 220513363"
        matrix = MATRIX(target, [1], ".", exposure_time=120)
        inject_dir = None
        try:
            inject_dir = matrix.inject(1, 5, 5, 1, 3, 3, 1)
            self.assertEquals(5, len(os.listdir(inject_dir)))
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
            matrix.recovery(inject_dir, 5, 0)
            matrix.plot_results(target, inject_dir)
            self.assertEquals(8, len(os.listdir(inject_dir)))
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)

    def test_inject_one_preserve(self):
        target = "TIC 220513363"
        matrix = MATRIX(target, [1], ".", True, exposure_time=120)
        inject_dir = None
        try:
            inject_dir = matrix.inject(1, 5, 5, 1, 3, 3, 1)
            matrix.recovery(inject_dir, 5, 0)
            matrix.plot_results(target, inject_dir)
            self.assertEquals(9, len(os.listdir(inject_dir)))
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)

    def test_inject_four(self):
        target = "TIC 220513363"
        matrix = MATRIX(target, [1], ".", exposure_time=120)
        inject_dir = None
        try:
            inject_dir = matrix.inject(1, 5, 5.1, 2, 3, 3.1, 2)
            self.assertEquals(8, len(os.listdir(inject_dir)))
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)

    def test_inject_multiphase(self):
        target = "TIC 220513363"
        matrix = MATRIX(target, [1], ".", exposure_time=120)
        inject_dir = None
        try:
            inject_dir = matrix.inject(2, 5, 5, 1, 3, 3, 1)
            self.assertEquals(6, len(os.listdir(inject_dir)))
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)

    def test_transit_mask(self):
        target = "TIC 305048087"
        matrix = MATRIX(target, [2], ".", exposure_time=120)
        inject_dir = matrix.retrieve_object_data()
        power_args = {"transit_template": "default", "period_min": 1, "period_max": 10,
                      "n_transits_min": 2, "show_progress_bar": True, "use_threads": multiprocessing.cpu_count() - 1,
                      "R_star": matrix.radius, "R_star_min": matrix.radiusmin, "R_star_max": matrix.radiusmax,
                      "M_star": matrix.mass, "M_star_min": matrix.massmin, "M_star_max": matrix.massmax}
        intransit = matrix.transit_masks([{"P": 5.43, "T0": 2458355.249756 - 2457000.0, "D": 1.172585}],
                                        matrix.lc_build.lc.time.value)
        model = transitleastsquares(matrix.lc_build.lc.time.value, matrix.lc_build.lc.flux.value)
        results = model.power(**power_args)
        self.assertAlmostEqual(5.43, results.period, 2)
        model = transitleastsquares(matrix.lc_build.lc.time.value[~intransit], matrix.lc_build.lc.flux.value[~intransit])
        results = model.power(**power_args)
        self.assertNotAlmostEqual(5.43, results.period, 2)

    def test_inject_inputs(self):
        target = "TIC 305048087"
        matrix = MATRIX(target, [2], ".", exposure_time=120)
        inject_dir = None
        with(pytest.raises(AssertionError)):
            inject_dir = matrix.inject(2, 5, 5.1, 2, 3, 3.1, "ho")
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir = matrix.inject(2, 5, 5.1, 2, 3, 3.1, -0.1)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir = matrix.inject(2, 5, 5.1, 2, 3, -3.1, 2)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir = matrix.inject(2, 5, 5.1, 2, -3, 3.1, 2)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir = matrix.inject(2, 5, 5.1, -2, 3, 3.1, 2)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir = matrix.inject(2, 5, -5.1, 2, 3, 3.1, 2)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir = matrix.inject(2, -5, 5.1, 2, 3, 3.1, 2)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir = matrix.inject(-2, 5, 5.1, 2, 3, 3.1, 2)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)

    def test_inject_dir(self):
        inject_dir1 = None
        inject_dir2 = None
        target = "TIC 305048087"
        matrix = MATRIX(target, [2], ".", exposure_time=120)
        try:
            inject_dir1 = matrix.inject(2, 5, 5, 1, 3, 3, 1)
            self.assertTrue(os.path.isdir(inject_dir1))
            inject_dir2 = matrix.inject(2, 5, 5, 1, 3, 3, 1)
            self.assertTrue(os.path.isdir(inject_dir2))
        finally:
            if inject_dir1 is not None:
                shutil.rmtree(inject_dir1)
            if inject_dir2 is not None:
                shutil.rmtree(inject_dir2)

    def test_star_info(self):
        target = "TIC 220513363"
        star_info = StarInfo(target, (0.2, 0.5), 2000, 1.2, None, None, 0.5, 0.1, 0.2, 0.7, 0.15, 0.05, None, None)
        matrix = MATRIX(target, [1], ".", False, star_info, exposure_time=120)
        inject_dir = None
        try:
            inject_dir = matrix.inject(1, 5, 5, 1, 3, 3, 1)
            self.assertEquals(5, len(os.listdir(inject_dir)))
            self.assertEquals((0.2, 0.5), matrix.star_info.ld_coefficients)
            self.assertEquals(2000, matrix.star_info.teff)
            self.assertAlmostEqual(0.7, matrix.star_info.mass)
            self.assertAlmostEqual(0.55, matrix.star_info.mass_min)
            self.assertAlmostEqual(0.75, matrix.star_info.mass_max)
            self.assertAlmostEqual(0.5, matrix.star_info.radius)
            self.assertAlmostEqual(0.4, matrix.star_info.radius_min)
            self.assertAlmostEqual(0.7, matrix.star_info.radius_max)
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)
        matrix = MATRIX(target, [1], ".", exposure_time=120)
        inject_dir = None
        try:
            inject_dir = matrix.inject(1, 5, 5, 1, 3, 3, 1)
            self.assertEquals(5, len(os.listdir(inject_dir)))
            self.assertEquals((0.1258, 0.235), matrix.star_info.ld_coefficients)
            self.assertEquals(31000.0, matrix.star_info.teff)
            self.assertAlmostEqual(0.47, matrix.star_info.mass)
            self.assertAlmostEqual(0.44, matrix.star_info.mass_min)
            self.assertAlmostEqual(0.5, matrix.star_info.mass_max)
            self.assertAlmostEqual(0.18, matrix.star_info.radius)
            self.assertAlmostEqual(0.076, matrix.star_info.radius_min)
            self.assertAlmostEqual(0.284, matrix.star_info.radius_max)
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
