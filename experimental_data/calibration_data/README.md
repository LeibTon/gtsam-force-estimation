# Cable Drive Length Dataset

This folder contains multi-camera reconstruction results for a cable-driven soft robot / continuum structure at different **cable drive lengths**.

## Robot Configuration

Key geometry and material parameters used by the Cosserat multi-segment model (single cable actuation, analytic Jacobians, no finite differences).

### Geometry

| Parameter | Symbol / variable | Value | Notes |
|-----------|-------------------|-------|-------|
| Number of segments | `n` | 10 | Constant-curvature Cosserat segments |
| Segment length | `Llist` | 0.027 m (27 mm) each | Uniform; total backbone length ≈ 0.27 m |
| Rod radius | `r` | 0.0005 m (0.5 mm) | 1 mm diameter backbone rod |
| Cable routing radius | `Rrad` | 0.008 m (8 mm) | Radial offset of the actuation cable |
| Actuation length | `Lact` | 0.27 m (full rod length) | The tendon is routed from the base disc to the tip disc (start to end), so its actuation length equals the total backbone length -- **not** 0.02 m as an earlier version of this table stated |
| Disk thickness | `disk_thickness` | 0.001 m (1 mm) | Disk thickness used in visualization |
| Disk twist angles | `alpha` | all zeros (`n+1` entries) | No pre-rotation / no intentional precurvature |
| Plot samples per segment | `N` | 11 | Used by `plot_cosserat_multiseg` |

### Material / stiffness

| Parameter | Symbol / variable | Value | Notes |
|-----------|-------------------|-------|-------|
| Young's modulus | `E` | 50e9 Pa (50 GPa) | Backbone rod |
| Poisson's ratio | `nu` | 0.30 | |
| Shear modulus | `G` | `E / (2*(1+nu))` | Derived |
| Second moment of area | `I` | `pi*r^4/4` | Bending |
| Polar moment | `J` | `pi*r^4/2` | Torsion |
| Bending stiffness | `K11`, `K22` | `E*I` | About local 1- and 2-axes |
| Torsional stiffness | `K33` | `G*J` | About local 3-axis |
| Coupling stiffness | `K12`, `K13`, `K23` | 0 | Uncoupled stiffness matrix |

Stiffness units: N·m².

### Loading / gravity

| Parameter | Symbol / variable | Value | Notes |
|-----------|-------------------|-------|-------|
| Segment mass | `mL` | 0.002 kg per segment | Uniform |
| Gravity | `g` | -9.81 m/s² | Along the model gravity axis |

### Base pose (for reconstruction / plotting)

| Parameter | Variable | Value | Notes |
|-----------|----------|-------|-------|
| Base position | `p0` | `[0, 0, 0]` | |
| Base orientation | `R0` | `eye(3)` with `R0(3,3) = -1` | Flips the local z-axis |

### Solver setup (reference)

| Item | Value |
|------|-------|
| Unknowns | Curvature/twist per segment `k` (`3 x n`) + scalar cable tension |
| Initial guess | `k0 = 0.1 * ones(3,n)`, `tension0 = 0.1` |
| Solver | `fsolve` with Levenberg–Marquardt |

## Folder Naming

Each experiment folder is named by drive length:

| Folder | Drive length |
|--------|--------------|
| `lact_0` | 0 |
| `lact_0.5` | 0.5 |
| `lact_1` | 1 |
| `lact_1.5` | 1.5 |
| `lact_2` | 2 |
| `lact_2.5` | 2.5 |
| `lact_3` | 3 |
| `lact_3.5` | 3.5 |
| `lact_4` | 4 |

Inside each folder:

- `lact_<L>.txt` — original reconstruction log (multiple repeated entries; content is essentially the same within one length)
- `images/` — captured images
- `annotation_records/` — per-frame annotation records

## MATLAB Data File

All useful fields from the first log entry of each drive length are collected in:

```text
cable_drive_data.m
```

Camera positions are **not** included. For force, only the force magnitude (`force_n`) is kept.

### Load data

```matlab
data = cable_drive_data();
```

### Struct fields

| Field | Type | Description |
|-------|------|-------------|
| `length` | scalar | Cable drive length |
| `points_3d` | `11 x 3` | Reconstructed 3D points `[x, y, z]` |
| `force_n` | scalar | Force magnitude in Newtons (N) |
| `timestamp` | char | Capture timestamp |

### Access examples

```matlab
% Load all entries
data = cable_drive_data();

% Index access (1 = length 0, 2 = length 0.5, ...)
data(1).length
data(1).points_3d
data(1).force_n

% Lookup by drive length
d = data([data.length] == 1.5);
pts = d.points_3d;
F = d.force_n;

% Collect all lengths and forces
L = [data.length];
F = [data.force_n];
plot(L, F, '-o');
xlabel('Cable drive length');
ylabel('Force (N)');
```

### Index mapping

| Index `i` | `data(i).length` |
|-----------|------------------|
| 1 | 0 |
| 2 | 0.5 |
| 3 | 1 |
| 4 | 1.5 |
| 5 | 2 |
| 6 | 2.5 |
| 7 | 3 |
| 8 | 3.5 |
| 9 | 4 |

## Notes

- Within each `lact_<L>.txt`, repeated Data Count entries share the same reconstructed points; force values may vary slightly over time. The MATLAB file uses the **first** entry for each length.
- Units for `points_3d` follow the original reconstruction log (same coordinate frame as Cam1 origin in the raw logs).

## Model Calibration

The nominal parameters above (`E = 50e9 Pa`, tendon at `angle_offset = 0`) do not
match this rig exactly. Two calibrations were fit against all 9 `lact_*`
trials jointly (pooled RMSE over all 99 backbone points), each correcting the
**tendon's angular position** and one of two candidate stiffness explanations:

### 1. Young's modulus + angle offset

Assumes the load-cell readings (`force_n`) are the true tendon tension, and the
rod is stiffer than the nominal value.

| Parameter | Calibrated value |
|-----------|-------------------|
| Young's modulus `E` | `8.413808e10 Pa` (≈ 84.1 GPa, vs. nominal 50 GPa) |
| Shear modulus `G` | `E / (2*(1+nu))` ≈ `3.236080e10 Pa` |
| Tendon `angle_offset` | `-21.0119°` (vs. nominal 0°) |
| Pooled RMSE | **2.716 mm** |

Files: `scripts/calibrate_youngs_modulus_and_angle_all.py` →
`scripts/calibrated_youngs_modulus_and_angle_all.npz`. This is the default used
in `scripts/config.py`.

### 2. Tendon tension scale + angle offset

Keeps `E` fixed at the nominal 50 GPa and instead assumes the load-cell reading
is proportionally off (e.g. pulley friction, sensor mounting) — true tension
`= k * force_n`.

| Parameter | Calibrated value |
|-----------|-------------------|
| Young's modulus `E` | `5.0e10 Pa` (fixed, not calibrated) |
| Tension scale factor `k` | `0.6947` |
| Tendon `angle_offset` | `-21.0188°` |
| Pooled RMSE | **4.651 mm** |

Files: `scripts/calibrate_tension_scale_and_angle_all.py` →
`scripts/calibrated_tension_scale_and_angle_all.npz`.

The Young's-modulus calibration fits noticeably better (2.72 mm vs. 4.65 mm
pooled RMSE) and is the one used by default. Both calibrations found nearly
identical `angle_offset` values (≈ -21°), confirming that the tendon's true
angular position on the disc — not the stiffness/tension split — is the
dominant, well-determined correction; a joint 3-parameter fit (E, tension
scale, and angle together) confirmed E and the tension scale are degenerate
with each other (only their product is identifiable from shape data), so nothing
is gained by calibrating both at once.

### Per-trial fit error

RMSE between the calibrated model's predicted 11-point backbone shape and the
camera-reconstructed shape (mapped into the model frame via the
camera→model rigid transform, `scripts/calibrate_camera_frame.py` →
`scripts/camera_to_model_transform.npz`), evaluated with
`scripts/evaluate_all_drive_lengths.py <calibration.npz>`:

| Drive length (mm) | Force (N) | RMSE, E-calibrated (mm) | RMSE, tension-scale-calibrated (mm) |
|---|---|---|---|
| 0.0 | 0.006 | 1.122 | 1.121 |
| 0.5 | 0.599 | 1.745 | 1.806 |
| 1.0 | 1.271 | 1.555 | 2.255 |
| 1.5 | 1.938 | 2.840 | 3.544 |
| 2.0 | 2.624 | 2.191 | 3.763 |
| 2.5 | 3.312 | 2.202 | 4.293 |
| 3.0 | 3.972 | 2.467 | 4.686 |
| 3.5 | 4.593 | 3.098 | 5.897 |
| 4.0 | 5.252 | 5.127 | 9.121 |
| **Pooled (all 9 trials)** | | **2.716** | **4.651** |
