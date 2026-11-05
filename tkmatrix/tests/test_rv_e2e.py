"""End-to-end tests running the real CLI on an RV block (issue #34). They need MAST access, as the other suites do."""
import glob
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def run_cli(tmp_path, photometry):
    rng = np.random.default_rng(11)
    time = np.sort(rng.uniform(2458000.0, 2458200.0, 60))
    rv_file = tmp_path / "rv.csv"
    pd.DataFrame({"bjd": time, "rv": rng.normal(0, 1.0, 60), "rv_err": np.full(60, 1.0)}).to_csv(rv_file, index=False)
    properties = {
        "TARGET": "TIC 220513363", "SECTORS": [1], "AUTHOR": "SPOC", "EXPOSURE_TIME": 120, "CPUS": 2,
        "RV": {"FILE": str(rv_file), "PHASES": 1, "MIN_PERIOD": 3.0, "MAX_PERIOD": 6.0, "STEPS_PERIOD": 3,
               "PERIOD_GRID_GEOM": "lin", "MIN_MASS": 5.0, "MAX_MASS": 60.0, "STEPS_MASS": 3,
               "MASS_GRID_GEOM": "lin", "MIN_PERIOD_SEARCH": 0.6, "MAX_PERIOD_SEARCH": 10,
               "PERIOD_MATCH_TOLERANCE": 0.02, "SNR_THRESHOLD": 5, "RUN_LIMIT": 2,
               "OVERSAMPLING": 1},
    }
    if photometry:
        properties.update({"PHASES": 1, "MIN_PERIOD": 4, "MAX_PERIOD": 6, "STEPS_PERIOD": 3, "MIN_RADIUS": 3,
                           "MAX_RADIUS": 5, "STEPS_RADIUS": 3, "OVERSAMPLING": 0.1, "SNR_THRESHOLD": 5})
    properties_file = tmp_path / "properties.yaml"
    properties_file.write_text(yaml.safe_dump(properties))
    environment = dict(os.environ, PYTHONPATH=ROOT, MPLBACKEND="Agg")
    return subprocess.run([sys.executable, "-m", "tkmatrix", "--properties", str(properties_file),
                           "--dir", str(tmp_path)], cwd=str(tmp_path), env=environment, capture_output=True,
                          text=True, timeout=1800)


def assert_rv_outputs(tmp_path, result):
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
    directories = glob.glob(str(tmp_path / "TIC220513363_ir_*"))
    assert len(directories) == 1
    inject_dir = directories[0]
    report = pd.read_csv(os.path.join(inject_dir, "a_rv_report.csv"))
    assert report["period_found"].min() >= 0.6
    assert report["found"].eq((report["period_found"] - report["period"]).abs() / report["period"] < 0.02).all()
    assert len(report) == 9
    assert report[report["mass"] == 60.0]["found"].all()
    assert os.path.exists(os.path.join(inject_dir, "rv_thresholds.png"))
    assert os.path.exists(os.path.join(inject_dir, "inj-rec-rv.png"))
    return inject_dir


def test_cli_rv_only(tmp_path):
    inject_dir = assert_rv_outputs(tmp_path, run_cli(tmp_path, photometry=False))
    assert not os.path.exists(os.path.join(inject_dir, "a_tls_report.csv"))


def test_cli_rv_and_photometry(tmp_path):
    inject_dir = assert_rv_outputs(tmp_path, run_cli(tmp_path, photometry=True))
    assert os.path.exists(os.path.join(inject_dir, "a_tls_report.csv"))
    assert os.path.exists(os.path.join(inject_dir, "inj-rec.png"))
