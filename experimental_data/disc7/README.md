# Cable Drive + Hanging Weight Dataset

Multi-camera reconstruction data for a cable-driven continuum manipulator under combined cable actuation and an external hanging-weight force applied at **disk 7**.

All 3D coordinates are expressed in the **Cam1-origin** frame.

## Naming Convention

Each trial folder is named:

```text
lact{L}-weight{W}-disc{D}
```

| Token | Meaning |
|-------|---------|
| `lact{L}` | Cable drive length |
| `weight{W}` | Hanging weight mass in grams (e.g. `weight10` = 10 g) |
| `disc{D}` | Disk index where the external force is applied |

### Disk indexing

Disks are counted from the **base**:

- **Disk 1** = base disk (fixed / proximal end)
- Disk numbers increase toward the tip
- **`disc7`** = external force attached at the **7th disk from the base**

### Trials included

| Folder | Drive length (`lact`) | Weight (g) | Force disk |
|--------|----------------------|------------|------------|
| `lact0-weight10-disc7` | 0 | 10 | 7 |
| `lact0-weight30-disc7` | 0 | 30 | 7 |
| `lact0-weight50-disc7` | 0 | 50 | 7 |
| `lact1-weight10-disc7` | 1 | 10 | 7 |
| `lact1-weight30-disc7` | 1 | 30 | 7 |
| `lact1-weight50-disc7` | 1 | 50 | 7 |
| `lact2-weight10-disc7` | 2 | 10 | 7 |
| `lact2-weight30-disc7` | 2 | 30 | 7 |
| `lact2-weight50-disc7` | 2 | 50 | 7 |
| `lact3-weight10-disc7` | 3 | 10 | 7 |
| `lact3-weight30-disc7` | 3 | 30 | 7 |
| `lact3-weight50-disc7` | 3 | 50 | 7 |
| `lact4-weight10-disc7` | 4 | 10 | 7 |
| `lact4-weight30-disc7` | 4 | 30 | 7 |

Each folder contains `reconstruction_log.txt`, `images/`, and `annotation_records/`.

## MATLAB Data File

```matlab
data = cable_drive_weight_data();
```

Uses the **first** log entry of each trial. Camera positions are not included.

| Field | Description |
|-------|-------------|
| `folder` | Trial folder name |
| `lact` | Cable drive length |
| `weight_g` | Hanging weight (g) |
| `disc` | Force-application disk index (base = 1) |
| `points_3d` | `11 x 3` reconstructed backbone points `[x y z]` (Cam1 origin) |
| `force_n` | Load-cell tendon/cable tension (N) |
| `external_force_n` | Hanging-weight force magnitude (N) = `weight_g * 9.81 / 1000` |
| `force_direction_unit_3d` | Unit direction of the force-application rope (`ROPE1`, Cam1 origin) |
| `force_direction_3d` | Non-unit rope direction (`ROPE1`) |
| `force_vector_n` | External force vector (N) = `external_force_n * force_direction_unit_3d` |
| `rope_points_3d` | Reconstructed rope points |
| `rope_length_mm` | Reconstructed rope length (mm) |
| `timestamp` | Capture timestamp |

```matlab
% Lookup by drive length and hanging weight
d = data([data.lact] == 2 & [data.weight_g] == 30);
Fdir = d.force_direction_unit_3d;
Fvec = d.force_vector_n;
```

## Robot Configuration (summary)

| Parameter | Value |
|-----------|-------|
| Segments `n` | 10 |
| Segment length | 0.027 m each |
| Rod radius | 0.0005 m |
| Cable routing radius `Rrad` | 0.008 m |
| Model actuation length `Lact` | 0.27 m (full rod length -- tendon is routed base-to-tip) |
| Young's modulus `E` | 50 GPa |
| Poisson's ratio `nu` | 0.30 |
| Segment mass `mL` | 0.002 kg |
| Gravity `g` | -9.81 m/s² |
| Base pose | `p0 = [0,0,0]`, `R0 = eye(3)` with `R0(3,3) = -1` |

Bending stiffness `K11 = K22 = E*I`, torsion `K33 = G*J`, couplings `K12 = K13 = K23 = 0`.
