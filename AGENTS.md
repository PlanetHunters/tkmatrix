# AGENTS.md — tkmatrix

> Agent-facing guide for the `tkmatrix` (MATRIX ToolKit) project.

## Before Making Changes

Always understand the user's question or command fully before writing any code. If there is any ambiguity or doubt about what to change, ask clarifying questions. Only make edits when the action is clear and unambiguous. Do not assume intent — confirm it.

## Project Intent

**tkmatrix** (MATRIX ToolKit) is an injection-recovery software for exoplanet research. It creates grids of synthetic transiting exoplanet signals (varying periods, radii, and epochs) injected into real light curves, then attempts to recover them using transit search algorithms (TLS, BLS, etc.). A key feature is **multi-phase analysis**, where multiple epochs are tested for each period-radius combination, improving accuracy over traditional 2D approaches.

- **Package:** `tkmatrix`
- **Version:** `0.13.1`
- **Authors:** M. Devora-Pajares & F.J. Pozuelos
- **License:** MIT
- **Repository:** https://github.com/PlanetHunters/tkmatrix
- **Documentation:** https://tkmatrix.readthedocs.io/

## Architecture

```
tkmatrix/
├── __init__.py              # Builds libellc.so via make; patches ellc module
├── __main__.py              # CLI entrypoint
├── tkmatrix_class.py        # MATRIX class — core engine
├── inject_model.py          # Light curve injection (batman + ellc)
├── inject_rv_model.py       # Radial Velocity injection
├── rv.py                    # RV recovery and fitting
├── bayesian.py              # Bayesian analysis utilities
├── custom_algorithms/        # Pluggable search algorithms
│   └── CustomSearchAlgorithm.py
├── ellc/                    # Bundled Fortran/C library (git submodule)
├── properties.yaml          # Default/reference configuration schema
├── programmatic.py           # Example programmatic API usage
└── tests/
    └── test_tkmatrix.py
```

### Key Classes

| Class | Role |
|-------|------|
| `MATRIX` | Core engine. Methods: `inject()`, `recovery()`, `plot_results()`. |
| `SearchInput` | Dataclass for search configuration. |
| `InjectModel` | Light curve injection logic using `batman` and `ellc`. |
| `InjectRvModel` | Radial velocity injection logic. |
| `RvFitter` | RV recovery and fitting. |
| `CustomSearchAlgorithm` | ABC for pluggable search algorithms (e.g., `BlsCustomSearchAlgorithm`). |

### Execution Flow

1. `__main__.py` parses CLI args, loads YAML properties, instantiates `MATRIX`.
2. `MATRIX.retrieve_object_data()` downloads/reads light curve via `lcbuilder`.
3. `MATRIX.inject()` creates a grid of period x radius x phase scenarios, spawns `InjectModel` per scenario, writes CSV files.
4. `MATRIX.recovery()` reads injected CSVs, runs TLS/BLS/custom search per file, writes `a_tls_report.csv`.
5. `MATRIX.plot_results()` reads reports, generates `inj-rec.png` heatmap.

## Entrypoints

No `console_scripts` in `setup.py`. Entry is via `python -m` module execution or library import.

### CLI Entrypoint

```bash
python3 -m tkmatrix --properties user-properties.yaml
```

Arguments:
- `--properties <path>` (required) — path to user YAML config
- `--dir <path>` (optional, default `./`) — working directory
- `--preserve` (optional flag) — keep injected curve CSVs after execution

### Programmatic / Library Entrypoint

```python
from tkmatrix.tkmatrix_class import MATRIX
from lcbuilder.star.starinfo import StarInfo
import lcbuilder.constants

target = "TIC 220513363"
matrix = MATRIX(
    target=target,
    sectors=[1],
    author=lcbuilder.constants.SPOC_AUTHOR,
    dir=".",
    exposure_time=120
)

inject_dir, period_grid, radius_grid = matrix.inject(
    phases=2,
    min_period=5, max_period=5, steps_period=1,
    min_radius=3, max_radius=3, steps_radius=1
)

matrix.recovery(
    inject_dir,
    snr_threshold=5,
    detrend_ws=0,
    detrend_method='biweight',
    fit_method='bls-periodogram',
    run_limit=2,
    custom_search_algorithm=None,
    min_period_search=0.5,
    max_period_search=25,
    oversampling=1,
    signal_selection_mode='period-epoch'
)

matrix.plot_results(target, inject_dir)
```

## Environment Setup

### Prerequisites

- Python `>= 3.11`
- `make`, `gfortran`, `gcc` — required to build the bundled `ellc` library (`libellc.so`)
- `conda` (tox-conda is used for tests)

### Conda Environment Setup

```bash
cd /path/to/tkmatrix
conda create -n tkmatrix311 python=3.11
conda activate tkmatrix311

# Install the package (triggers make for libellc.so)
pip install -e .
# or install pinned deps first
pip install -r requirements.txt
```

### Docker (alternative)

```bash
docker build -t tkmatrix .
docker run -it -v $(pwd):/work tkmatrix
```

### Key Dependencies

| Package | Purpose |
|---------|---------|
| `lcbuilder==0.26.0` | Light curve building, object info, star catalogs |
| `foldedleastsquares==1.1.11` | TLS search engine |
| `lightkurve==2.5.0` | Mission data download |
| `batman-package==2.5.3` | Transit modeling |
| `numpy==2.1.1` | Core numerics |
| `astropy==7.0.1`, `astroquery==0.4.10` | Astronomy stack |
| `matplotlib==3.10.1`, `seaborn==0.13.2`, `corner==2.2.3` | Plotting |
| `torch==2.7.0` | PyTorch |
| `Cython==3.0.6`, `numba==0.61.2`, `pybind11==2.11.1` | Compiled extensions |

## Build & Test

### Run Tests (tox — recommended)

```bash
tox
```

- Uses `tox-conda` to provision a Python 3.11 environment.
- Installs `numpy==2.1.1`, `pytest`, `setuptools`, `wheel`, `Cython`.
- Command: `pytest -v tkmatrix/tests/`

### Run Tests (pytest directly)

```bash
pytest -v tkmatrix/tests/
```

### Test Suite

Located in `tkmatrix/tests/test_tkmatrix.py`:
- `test_inject_one` — basic injection + recovery
- `test_inject_9_preserve` — injection with preserve and plotting
- `test_inject_four`, `test_inject_multiphase` — multi-phase injection
- `test_inject_inputs` — bad input validation
- `test_inject_dir` — directory collision handling
- `test_star_info` — custom vs catalog star info
- `test_inject_grids` — custom period/radius grids

Tests require internet access (downloads TESS/Kepler data) and take non-trivial time.

## Execution Recipes

### Use Case A: CLI Execution (Recommended)

1. Create a properties YAML file based on `user-properties.yaml` or files in `examples/`.
2. Run:
   ```bash
   python3 -m tkmatrix --properties my-properties.yaml
   ```
3. Optional flags:
   ```bash
   python3 -m tkmatrix --properties my-properties.yaml --preserve --dir /path/to/output
   ```

### Use Case B: Programmatic / Notebook Execution

```python
from tkmatrix.tkmatrix_class import MATRIX
from lcbuilder.star.starinfo import StarInfo
import lcbuilder.constants

target = "TIC 220513363"
matrix = MATRIX(
    target=target,
    sectors=[1],
    author=lcbuilder.constants.SPOC_AUTHOR,
    dir=".",
    exposure_time=120
)

inject_dir, period_grid, radius_grid = matrix.inject(
    phases=2,
    min_period=5, max_period=5, steps_period=1,
    min_radius=3, max_radius=3, steps_radius=1
)

matrix.recovery(
    inject_dir,
    snr_threshold=5,
    detrend_ws=0,
    detrend_method='biweight',
    fit_method='bls-periodogram',
    run_limit=2,
    custom_search_algorithm=None,
    min_period_search=0.5,
    max_period_search=25,
    oversampling=1,
    signal_selection_mode='period-epoch'
)

matrix.plot_results(target, inject_dir)
```

### Use Case C: Docker Execution

```bash
docker build -t tkmatrix .
docker run -it -v $(pwd):/work tkmatrix
# inside container:
python3 -m tkmatrix --properties /work/user-properties.yaml
```

### Use Case D: RV Injection & Recovery

The properties YAML supports an `RV:` block. When `RV.FILE` is provided, the tool will:
1. Run `recovery_rv_periods()` on the original RV file.
2. Inject synthetic RV signals (`inject_rv`).
3. Run RV recovery (`recovery_rv`).
4. Plot RV results.
