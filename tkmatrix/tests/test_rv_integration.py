"""Integration tests for the MATRIX RV stages (issue #34). Run offline: only star loading is stubbed."""
import inspect
import os
import re
import signal
import subprocess
import sys
from concurrent.futures.process import BrokenProcessPool
from unittest import mock

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from lcbuilder.star.starinfo import StarInfo

from tkmatrix.rv import RvFitter
from tkmatrix.tkmatrix_class import MATRIX

REPORT_COLUMNS = ['period', 'mass', 'epoch', 'period_found', 'epoch_found', 'mass_found', 'found', 'snr', 'sde', 'run',
                  'fap']
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
    assert report["fap"].between(0, 1).all()


def test_results_plot_is_written(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    inject_dir, _, _ = inject(matrix, rv_file, phases=1)
    matrix.recovery_rv(inject_dir, None, 5, 3, 12, 1)
    MATRIX.plot_results("TIC 1", inject_dir, is_rv=True)
    assert os.path.exists(os.path.join(inject_dir, "inj-rec-rv.png"))


def write_report(inject_dir, periods, masses, found_cell):
    rows = [{"period": period, "mass": mass, "found": (period, mass) == found_cell, "sde": 1.0}
            for period in periods for mass in masses]
    pd.DataFrame(rows).to_csv(os.path.join(inject_dir, "a_rv_report.csv"), index=False)


def plotted_grid(inject_dir):
    """Runs plot_results and returns the axes, the cell edges and the plotted rates of the heat map."""
    with mock.patch.object(plt, "close"):
        MATRIX.plot_results("TIC 1", inject_dir, is_rv=True)
    axes = plt.gca().figure.axes[0]
    mesh = axes.collections[0]
    corners = mesh.get_coordinates()
    return axes, corners[0, :, 0], corners[:, 0, 1], np.asarray(mesh.get_array()).reshape(corners.shape[0] - 1, -1)


@pytest.mark.parametrize("periods, scale", [(np.round(np.logspace(0, 2, 7), 2), "log"),
                                            (np.linspace(2.0, 14.0, 7), "linear")])
def test_every_grid_value_owns_the_cell_drawn_around_it(tmp_path, periods, scale):
    masses = [1.0, 2.0, 3.0, 4.0]
    write_report(str(tmp_path), periods, masses, (periods[2], masses[1]))
    axes, x_edges, y_edges, rates = plotted_grid(str(tmp_path))
    plt.close("all")
    assert axes.get_xscale() == scale
    assert axes.get_yscale() == "linear"
    assert rates.shape == (len(masses), len(periods))
    assert np.argwhere(rates == 100).tolist() == [[1, 2]]
    assert x_edges[2] < periods[2] < x_edges[3]
    assert y_edges[1] < masses[1] < y_edges[2]


@pytest.mark.parametrize("planet_period, left_limit", [(5.0, None), (0.5, 0.5)])
def test_planet_markers_leave_no_margin_around_the_map(tmp_path, planet_period, left_limit):
    periods = np.round(np.logspace(0, 2, 7), 2)
    masses = [1.0, 2.0, 3.0, 4.0]
    write_report(str(tmp_path), periods, masses, (periods[2], masses[1]))
    planets = pd.DataFrame([{"period": planet_period, "mass": 2.0, "mass_err_bottom": 0.1, "mass_err_up": 0.1}])
    with mock.patch.object(plt, "close"):
        MATRIX.plot_results("TIC 1", str(tmp_path), is_rv=True, planets_df=planets)
    axes = plt.gca().figure.axes[0]
    edges = axes.collections[0].get_coordinates()[0, :, 0]
    limits = axes.get_xlim()
    plt.close("all")
    assert limits == pytest.approx((left_limit or edges[0], edges[-1]))


def test_diff_plot_is_written_for_a_log_grid(tmp_path):
    periods = np.round(np.logspace(0, 2, 7), 2)
    for name, found in (("one", True), ("two", False)):
        os.mkdir(tmp_path / name)
        rows = [{"period": period, "radius": radius, "found": found, "sde": 1.0}
                for period in periods for radius in (1.0, 2.0, 3.0)]
        pd.DataFrame(rows).to_csv(tmp_path / name / "a_tls_report.csv", index=False)
    MATRIX.plot_diff("TIC 1", str(tmp_path / "one"), str(tmp_path / "two"), str(tmp_path))
    assert os.path.exists(tmp_path / "inj-rec-diff.png")


def test_the_last_grid_values_keep_separate_cells(tmp_path):
    periods = [2.0, 4.0, 6.0]
    masses = [1.0, 2.0, 3.0]
    write_report(str(tmp_path), periods, masses, (6.0, 3.0))
    rates = plotted_grid(str(tmp_path))[3]
    plt.close("all")
    assert rates[-1, -1] == 100
    assert rates[-1, -2] == 0
    assert rates[-2, -1] == 0


def test_rv_search_settings_reach_the_recovery(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    inject_dir, _, _ = inject(matrix, rv_file)
    with mock.patch("tkmatrix.rv.RvFitter.recover_strongest_signal", wraps=RvFitter.recover_strongest_signal) as search:
        matrix.search_input.cores = 1
        scenario = (os.path.join(inject_dir, sorted(f for f in os.listdir(inject_dir) if f.startswith("RV_P"))[0]),
                    3.0, 2.0, 1.0, None, 0.5, 100, 0.7, 12, 5, 1, 0.2, 0.5)
        report, failure = MATRIX.recover_rv_scenario(scenario)
    assert failure is None
    assert search.call_args.args[2] == 100
    assert search.call_args.args[3] == 0.7
    assert report["found"] == (abs(report["period_found"] - 3.0) / 3.0 < 0.2 and report["fap"] < 0.5)


@pytest.mark.parametrize("settings", [{"min_period_search": 13}, {"min_period_search": 0},
                                      {"period_match_tolerance": 0}])
def test_invalid_rv_search_settings_are_rejected_before_searching(matrix, tmp_path, settings):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    inject_dir, _, _ = inject(matrix, rv_file)
    with pytest.raises(ValueError, match="period"):
        matrix.recovery_rv(inject_dir, None, 5, 3, 12, 1, **settings)


@pytest.mark.parametrize("threshold", [0, -0.1])
def test_invalid_fap_threshold_is_rejected_before_searching(matrix, tmp_path, threshold):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    inject_dir, _, _ = inject(matrix, rv_file)
    with pytest.raises(ValueError, match="FAP threshold"):
        matrix.recovery_rv(inject_dir, None, 5, 3, 12, 1, fap_threshold=threshold)


def test_recovery_rv_periods_uses_the_minimum_period(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    with mock.patch("tkmatrix.tkmatrix_class.RvFitter.recover_periods", wraps=RvFitter.recover_periods) as search:
        matrix.recovery_rv_periods(rv_file, 12, cores=2, min_period_search=0.8)
    assert search.call_args.kwargs["period_min"] == 0.8
    assert search.call_args.kwargs["steps_period"] == int((12 - 0.8) * 100)


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


@pytest.mark.parametrize("content,message", [
    ("bjd,rv,rv_err\n", "valid rows"),
    ("bjd,rv,rv_err\n1,1,1\n2,2,1\n3,3,1\n", "valid rows"),
    ("bjd,rv,rv_err\n1,inf,1\n2,2,1\n3,3,1\n4,4,1\n", "column rv "),
    ("bjd,rv,rv_err\n1,1,inf\n2,2,1\n3,3,1\n4,4,1\n", "column rv_err"),
    ("bjd,rv,rv_err\n1,text,1\n2,2,1\n3,3,1\n4,4,1\n", "column rv "),
    ("bjd,rv,rv_err\n1,1,1\n1,2,1\n1,3,1\n1,4,1\n", "identical"),
])
def test_unusable_rv_files_name_the_problem(tmp_path, content, message):
    rv_file = tmp_path / "rv.csv"
    rv_file.write_text(content)
    with pytest.raises(ValueError, match=message):
        MATRIX.load_rv_file(str(rv_file))


def test_whole_number_errors_work_with_a_mass_mask(tmp_path):
    rv_file = tmp_path / "rv.csv"
    frame = write_rv_file(str(rv_file))
    frame["rv_err"] = 2
    frame.to_csv(rv_file, index=False)
    loaded = MATRIX.load_rv_file(str(rv_file))
    assert np.all(np.isfinite(RvFitter.mask_signals(loaded, [{"P": 5.0, "M": 3.0}], 0.5)))


def test_invalid_file_is_reported_before_retrieving_the_object(matrix, tmp_path):
    rv_file = tmp_path / "rv.csv"
    rv_file.write_text("bjd,rv\n1,1\n")
    with pytest.raises(ValueError, match="rv_err"):
        matrix.recovery_rv_periods(str(rv_file), 12)
    MATRIX.retrieve_object_data.assert_not_called()


def test_when_all_scenarios_fail_the_cause_is_raised_and_the_files_are_kept(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    inject_dir, _, _ = inject(matrix, rv_file)
    with mock.patch("tkmatrix.rv.RvFitter.recover_signal", side_effect=ValueError("engine failure")):
        with pytest.raises(RuntimeError, match="engine failure"):
            matrix.recovery_rv(inject_dir, None, 5, 3, 12, 1)
    assert len([name for name in os.listdir(inject_dir) if name.startswith("RV_P")]) == 8


def test_partial_failures_are_logged_and_keep_the_other_rows(matrix, tmp_path, caplog):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    inject_dir, _, _ = inject(matrix, rv_file)
    broken = sorted(name for name in os.listdir(inject_dir) if name.startswith("RV_P"))[0]
    with open(os.path.join(inject_dir, broken), "w") as file:
        file.write("not,a,valid,file\n")
    with caplog.at_level("WARNING"):
        matrix.recovery_rv(inject_dir, None, 5, 3, 12, 1)
    assert len(pd.read_csv(os.path.join(inject_dir, "a_rv_report.csv"))) == 7
    assert "1 of 8 RV scenarios failed" in caplog.text
    assert broken in caplog.text


def test_report_is_saved_while_recovery_is_running(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    inject_dir, _, _ = inject(matrix, rv_file)
    saved_rows = []
    with mock.patch.object(MATRIX, "RV_REPORT_SAVE_SECONDS", 0), \
            mock.patch.object(MATRIX, "save_rv_report", side_effect=lambda rows, _: saved_rows.append(len(rows))):
        matrix.recovery_rv(inject_dir, None, 5, 3, 12, 1)
    assert saved_rows == [1, 2, 3, 4, 5, 6, 7, 8, 8]


def test_report_is_not_rewritten_before_the_save_interval_elapses(matrix, tmp_path):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    inject_dir, _, _ = inject(matrix, rv_file)
    saved_rows = []
    with mock.patch.object(MATRIX, "RV_REPORT_SAVE_SECONDS", 3600), \
            mock.patch.object(MATRIX, "save_rv_report", side_effect=lambda rows, _: saved_rows.append(len(rows))):
        matrix.recovery_rv(inject_dir, None, 5, 3, 12, 1)
    assert saved_rows == [8]


def test_a_dead_worker_raises_and_keeps_the_finished_rows(matrix, tmp_path, caplog):
    rv_file = str(tmp_path / "rv.csv")
    write_rv_file(rv_file)
    inject_dir, _, _ = inject(matrix, rv_file)
    matrix.search_input.cores = 1
    recover_signal = RvFitter.recover_signal

    def kill_on_the_long_period(df, period, *args, **kwargs):
        if period > 5:
            os.kill(os.getpid(), signal.SIGKILL)
        return recover_signal(df, period, *args, **kwargs)

    with mock.patch("tkmatrix.rv.RvFitter.recover_signal", side_effect=kill_on_the_long_period):
        with caplog.at_level("ERROR"), pytest.raises(BrokenProcessPool):
            matrix.recovery_rv(inject_dir, None, 5, 3, 12, 1)
    assert len(pd.read_csv(os.path.join(inject_dir, "a_rv_report.csv"))) == 4
    assert "RV recovery stopped after 4 of 8 scenarios" in caplog.text


def test_built_package_contains_ravex_and_its_license(tmp_path):
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    subprocess.run([sys.executable, "setup.py", "-q", "build_py", "--build-lib", str(tmp_path)], cwd=root, check=True,
                   capture_output=True)
    assert (tmp_path / "tkmatrix" / "ravex" / "src" / "ravex" / "core.py").exists()
    assert (tmp_path / "tkmatrix" / "ravex" / "LICENSE").exists()


def test_public_signatures_keep_the_contract():
    assert list(inspect.signature(MATRIX.inject_rv).parameters)[:13] == [
        'self', 'inject_dir', 'rv_file', 'phases', 'min_period', 'max_period', 'steps_period', 'min_mass', 'max_mass',
        'steps_mass', 'period_grid', 'mass_grid', 'period_grid_geom']
    assert list(inspect.signature(MATRIX.recovery_rv).parameters)[:7] == [
        'self', 'inject_dir', 'rv_masks', 'snr_threshold', 'run_limit', 'max_period_search', 'oversampling']
    assert list(inspect.signature(MATRIX.recovery_rv_periods).parameters)[:6] == [
        'self', 'filename', 'max_period_search', 'rv_masks', 'oversampling', 'cores']
