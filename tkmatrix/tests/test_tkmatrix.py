import os
import shutil
import unittest

import lcbuilder.constants
import pytest
from lcbuilder.star.starinfo import StarInfo

from tkmatrix.tkmatrix_class import MATRIX


class TestsMatrix(unittest.TestCase):
    def test_inject_one(self):
        target = "TIC 220513363"
        matrix = MATRIX(target, [1], lcbuilder.constants.SPOC_AUTHOR, ".", exposure_time=120)
        inject_dir = None
        try:
            inject_dir, period_grid, radius_grid = matrix.inject(1, 5, 5, 1, 3, 3, 1)
            self.assertEquals(10, len(os.listdir(inject_dir)))
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
            matrix.recovery(inject_dir, 5, detrend_ws=0, oversampling=0.1)
            matrix.plot_results(target, inject_dir)
            self.assertEquals(9, len(os.listdir(inject_dir)))
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)

    def test_inject_one_preserve(self):
        target = "TIC 220513363"
        matrix = MATRIX(target, [1], lcbuilder.constants.SPOC_AUTHOR, ".", True, exposure_time=120)
        inject_dir = None
        try:
            inject_dir, period_grid, radius_grid = matrix.inject(1, 5, 5, 1, 3, 3, 1)
            matrix.recovery(inject_dir, 5, detrend_ws=0, oversampling=0.1)
            matrix.plot_results(target, inject_dir)
            self.assertEquals(10, len(os.listdir(inject_dir)))
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)

    def test_inject_four(self):
        target = "TIC 220513363"
        matrix = MATRIX(target, [1], lcbuilder.constants.SPOC_AUTHOR, ".", exposure_time=120)
        inject_dir = None
        try:
            inject_dir, period_grid, radius_grid = matrix.inject(1, 5, 5.1, 2, 3, 3.1, 2)
            self.assertEquals(13, len(os.listdir(inject_dir)))
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)

    def test_inject_multiphase(self):
        target = "TIC 220513363"
        matrix = MATRIX(target, [1], lcbuilder.constants.SPOC_AUTHOR, ".", exposure_time=120)
        inject_dir = None
        try:
            inject_dir, period_grid, radius_grid = matrix.inject(2, 5, 5, 1, 3, 3, 1)
            self.assertEquals(11, len(os.listdir(inject_dir)))
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)

    def test_inject_inputs(self):
        target = "TIC 305048087"
        matrix = MATRIX(target, [2], lcbuilder.constants.SPOC_AUTHOR, ".", exposure_time=120)
        inject_dir = None
        with(pytest.raises(AssertionError)):
            inject_dir, period_grid, radius_grid = matrix.inject(2, 5, 5.1, 2, 3, 3.1, "ho")
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir, period_grid, radius_grid = matrix.inject(2, 5, 5.1, 2, 3, 3.1, -0.1)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir, period_grid, radius_grid = matrix.inject(2, 5, 5.1, 2, 3, -3.1, 2)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir, period_grid, radius_grid = matrix.inject(2, 5, 5.1, 2, -3, 3.1, 2)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir, period_grid, radius_grid = matrix.inject(2, 5, 5.1, -2, 3, 3.1, 2)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir, period_grid, radius_grid = matrix.inject(2, 5, -5.1, 2, 3, 3.1, 2)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir, period_grid, radius_grid = matrix.inject(2, -5, 5.1, 2, 3, 3.1, 2)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)
        with(pytest.raises(AssertionError)):
            inject_dir, period_grid, radius_grid = matrix.inject(-2, 5, 5.1, 2, 3, 3.1, 2)
        if inject_dir is not None:
            shutil.rmtree(inject_dir, ignore_errors=True)

    def test_inject_dir(self):
        inject_dir1 = None
        inject_dir2 = None
        target = "TIC 305048087"
        matrix = MATRIX(target, [2], lcbuilder.constants.SPOC_AUTHOR, ".", exposure_time=120)
        try:
            inject_dir1, period_grid, radius_grid = matrix.inject(2, 5, 5, 1, 3, 3, 1)
            self.assertTrue(os.path.isdir(inject_dir1))
            inject_dir2, period_grid, radius_grid = matrix.inject(2, 5, 5, 1, 3, 3, 1)
            self.assertTrue(os.path.isdir(inject_dir2))
        finally:
            if inject_dir1 is not None:
                shutil.rmtree(inject_dir1, ignore_errors=True)
            if inject_dir2 is not None:
                shutil.rmtree(inject_dir2, ignore_errors=True)

    def test_star_info(self):
        target = "TIC 220513363"
        star_info = StarInfo(target, (0.2, 0.5), 2000, 1.2, None, None, 0.5, 0.1, 0.2, 0.7, 0.15, 0.05, None, None)
        matrix = MATRIX(target, [1], lcbuilder.constants.SPOC_AUTHOR, ".", False, star_info, exposure_time=120)
        inject_dir = None
        try:
            inject_dir, period_grid, radius_grid = matrix.inject(1, 5, 5, 1, 3, 3, 1)
            self.assertEquals(10, len(os.listdir(inject_dir)))
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
        matrix = MATRIX(target, [1], lcbuilder.constants.SPOC_AUTHOR, ".", exposure_time=120)
        inject_dir = None
        try:
            inject_dir, period_grid, radius_grid = matrix.inject(1, 5, 5, 1, 3, 3, 1)
            self.assertEquals(10, len(os.listdir(inject_dir)))
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

    def test_inject_grids(self):
        target = "TIC 220513363"
        matrix = MATRIX(target, [1], lcbuilder.constants.SPOC_AUTHOR, ".", exposure_time=120)
        inject_dir = None
        period_grid_expected = [1, 2, 3, 8, 20]
        radius_grid_expected = [1.1, 1.4, 2, 2.5]
        try:
            inject_dir, period_grid, radius_grid = matrix.inject(1, 5, 5.1, 2, 3, 3.1, 2, period_grid=period_grid_expected,
                                                                 radius_grid=radius_grid_expected)
            self.assertEquals(29, len(os.listdir(inject_dir)))
            self.assertEquals(period_grid_expected, period_grid)
            self.assertEquals(radius_grid_expected, radius_grid)
        finally:
            if inject_dir is not None:
                shutil.rmtree(inject_dir, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
