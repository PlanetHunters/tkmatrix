"""Integration tests for the MATRIX RV stages (issue #34). Run offline: only star loading is stubbed."""
import inspect
import os
import re
from unittest import mock

import astropy.units as u
import numpy as np
import pandas as pd
import pytest
from lcbuilder.star.starinfo import StarInfo

from tkmatrix.tkmatrix_class import MATRIX

REPORT_COLUMNS = ['period', 'mass', 'epoch', 'period_found', 'epoch_found', 'mass_found', 'found', 'snr', 'sde', 'run']
PERIOD = 6.3
K = 8.0


def write_rv_file(path, shuffled=False, points=70):
    rng = np.random.default_rng(5)
    time = np.sort(rng.uniform(2458000.0, 2458250.0, points))
    frame = pd.DataFrame({"bjd": time, "rv": rng.normal(0, 1.0, points), "rv_err": np.full(points, 1.0)})
    if shuffled:
        frame = frame.sample(frac=1, random_state=3)
    frame.to_csv(path, index=False)
    return frame


@pytest.fixture
def matrix(tmp_path):
    star = StarInfo("TIC 1", (0.3, 0.3), 3500, None, 4.7, 0.1, 0.5, None, None, 0.5, None, None, None, None)
    instance = MATRIX("TIC 1", [1], "SPOC", str(tmp_path), star_info=star, exposure_time=120, cores=2)
    inject_dir = str(tmp_path / "TIC1_ir_0") + "/"
    os.mkdir(inject_dir)

    def loaded(search_input, inject_dir_arg=None):
        search_input.mstar = 0.5 * u.M_sun
        search_input.rstar = 0.5 * u.R_sun
        return inject_dir, None, None, search_input

    with mock.patch.object(MATRIX, "retrieve_object_data", side_effect=loaded):
        yield instance


def inject(matrix, rv_file, phases=2, masses=(2.0, 40.0), **kwargs):
    return matrix.inject_rv(None, rv_file, phases, 3.0, 8.0, 2, min(masses), max(masses), len(masses), **kwargs)


def test_recovery_rv_periods_writes_the_diagnostic_plot_and_returns_the_directory(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    inject_dir = matrix.recovery_rv_periods(rv_file, 12, cores=2)
    assert os.path.exists(os.path.join(inject_dir, "rv_thresholds.png"))


def test_injected_files_keep_the_contract(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    inject_dir, period_grid, mass_grid = inject(matrix, rv_file)
    names = [name for name in os.listdir(inject_dir) if name.startswith("RV_P")]
    assert len(names) == 2 * 2 * 2
    assert all(re.fullmatch(r"RV_P\d{3}\.\d{2}_M\d{2}\.\d{2}_T\d+\.\d{2}\.csv", name) for name in names)
    assert list(pd.read_csv(os.path.join(inject_dir, names[0])).columns) == ['bjd', 'rv', 'rv_err']


def test_recovery_writes_one_contract_row_per_scenario_and_finds_the_massive_planet(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    inject_dir, _, _ = inject(matrix, rv_file)
    matrix.recovery_rv(inject_dir, None, 5, 3, 12, 1)
    report = pd.read_csv(os.path.join(inject_dir, "a_rv_report.csv"))
    assert list(report.columns) == REPORT_COLUMNS
    assert len(report) == 8
    assert report[report["mass"] == 40.0]["found"].all()
    assert report["run"].min() >= 0


def test_results_plot_is_written(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    # plot_results uses the grid values as bin edges, so it needs at least three of each
    inject_dir, _, _ = matrix.inject_rv(None, rv_file, 1, 3.0, 8.0, 3, 2.0, 40.0, 3)
    matrix.recovery_rv(inject_dir, None, 5, 3, 12, 1)
    MATRIX.plot_results("TIC 1", inject_dir, is_rv=True)
    assert os.path.exists(os.path.join(inject_dir, "inj-rec-rv.png"))


def test_injected_files_are_removed_unless_preserved(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    inject_dir, _, _ = inject(matrix, rv_file)
    matrix.recovery_rv(inject_dir, None, 5, 3, 12, 1)
    assert not [name for name in os.listdir(inject_dir) if name.startswith("RV_P")]


def test_epochs_are_anchored_to_the_earliest_time_for_unsorted_files(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    frame = write_rv_file(rv_file, shuffled=True)
    inject_dir, _, _ = inject(matrix, rv_file, phases=1, masses=(5.0, 6.0))
    epochs = [float(re.search(r"_T([0-9.]+)\.csv", name)[1]) for name in os.listdir(inject_dir)
              if name.startswith("RV_P")]
    assert min(epochs) >= frame["bjd"].min()
    assert max(epochs) <= frame["bjd"].min() + 8.0


def test_log_mass_grid_uses_the_mass_steps(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    _, _, mass_grid = matrix.inject_rv(None, rv_file, 1, 3.0, 8.0, 2, 1.0, 100.0, 3, mass_grid_geom="log")
    assert len(mass_grid) == 3


def test_rows_with_missing_values_are_dropped_before_injection(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    frame = write_rv_file(rv_file)
    frame.loc[frame.index[3], "rv"] = np.nan
    frame.to_csv(rv_file, index=False)
    inject_dir, _, _ = inject(matrix, rv_file, phases=1, masses=(5.0, 6.0))
    name = next(name for name in os.listdir(inject_dir) if name.startswith("RV_P"))
    assert len(pd.read_csv(os.path.join(inject_dir, name))) == len(frame) - 1


@pytest.mark.parametrize("bad", ["missing_column", "zero_error"])
def test_invalid_rv_files_are_rejected_with_a_clear_error(matrix, tmp_path, bad):
    rv_file = str(tmp_path / "rv.csv")
    frame = write_rv_file(rv_file)
    if bad == "missing_column":
        frame.drop(columns=["rv_err"]).to_csv(rv_file, index=False)
    else:
        frame.loc[frame.index[0], "rv_err"] = 0
        frame.to_csv(rv_file, index=False)
    with pytest.raises(ValueError, match="rv_err"):
        inject(matrix, rv_file, phases=1, masses=(5.0, 6.0))


def test_public_signatures_keep_the_contract():
    assert list(inspect.signature(MATRIX.inject_rv).parameters)[:13] == [
        'self', 'inject_dir', 'rv_file', 'phases', 'min_period', 'max_period', 'steps_period', 'min_mass', 'max_mass',
        'steps_mass', 'period_grid', 'mass_grid', 'period_grid_geom']
    assert list(inspect.signature(MATRIX.recovery_rv).parameters) == [
        'self', 'inject_dir', 'rv_masks', 'snr_threshold', 'run_limit', 'max_period_search', 'oversampling']
    assert list(inspect.signature(MATRIX.recovery_rv_periods).parameters)[:6] == [
        'self', 'filename', 'max_period_search', 'rv_masks', 'oversampling', 'cores']
