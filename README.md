
# Installation

## Requirements

- Linux or WSL (this project is built and run under WSL on Windows).
- Python 3.12.
- CMake >= 3.22, a C++17 compiler (e.g. GCC 11).
- Eigen3 (`libeigen3-dev`).
- Boost 1.74+ (`libboost-all-dev`).
- **GTSAM 4.3a2**, built from source with Python bindings disabled (this
  project uses its own pybind11 bindings in `src/bindings.cpp`, not
  GTSAM's). Install it to a prefix and point `GTSAM_INCLUDE_DIR` /
  `GTSAM_LIB_DIR` at it when building (see below); it defaults to
  `/usr/local/include` / `/usr/local/lib` if built with the default prefix.
- `pybind11` (installed automatically as a build dependency, see below).

## Build

```bash
# Build and install the tendon_robot extension 
# GTSAM_INCLUDE_DIR / GTSAM_LIB_DIR only needed if GTSAM was installed to a
# non-default prefix.
pip install -e .
```

Verify the build:

```bash
python3 -c "import tendon_robot; print(tendon_robot.GaussianBasisForceEstimator)"
```



# Running experiments

Run from `scripts/experiments/`:

- `camera_transform.py` — camera-to-model transform helpers (imported, not run).
- `experiment_parameters.py` — shared constants (imported, not run).
- `calibrate_experiment.py` — fits Young's modulus + friction from `experimental_data/calibration_data/`.
- `estimate_experiment.py` — runs the estimator on one trial (library; demo under `__main__`).
- `run_single_force_batch.py` / `run_two_force_batch.py` — run all benchmark cases and write results to this folder.

# Running simulations

`simulation_data/` holds the SoRoSim-generated ground-truth datasets
(`generateDataset.m` + `solveSample.m`, run in MATLAB with SoRoSim V6.3 on the
path) and `manipulator_new.mat`, the rod model used to generate them. Run the
Python side from `scripts/simulation/`:

- `simulation_parameters.py` — shared constants (imported, not run).
- `estimate_simulation.py` — runs the estimator on one SoRoSim sample (library; demo under `__main__`).
- `run_single_force_batch.py` / `run_two_force_batch.py` — run the estimator over every sample in `simulation_data/dataset_1forces.h5` / `dataset_2forces.h5` and write results to this folder.
