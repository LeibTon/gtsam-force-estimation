# Cable Drive + Hanging Weight Dataset

Single-force experiments with the external force applied at **disk 10**. Disk 1 is the fixed base disk, and disk numbers increase toward the tip.

## Naming

`lact{L}-weight{W}-disc{D}` means cable drive length `L`, hanging mass `W` in grams, and force applied at disk `D`.

## Trials

- `lact3-weight10-disc10`: lact = 3, weight = 10 g, disk = 10
- `lact3-weight20-disc10`: lact = 3, weight = 20 g, disk = 10
- `lact3-weight30-disc10`: lact = 3, weight = 30 g, disk = 10
- `lact4-weight10-disc10`: lact = 4, weight = 10 g, disk = 10
- `lact4-weight20-disc10`: lact = 4, weight = 20 g, disk = 10
- `lact4-weight30-disc10`: lact = 4, weight = 30 g, disk = 10

Each trial retains the key names `reconstruction_log.txt`, `images/`, and `annotation_records/`.

## MATLAB

```matlab
data = cable_drive_weight_data();
d = data([data.lact] == 3 & [data.weight_g] == 20);
```

The MATLAB file uses the first log record from each trial. It includes reconstructed backbone points, measured tendon tension (`force_n`), hanging-weight force magnitude, force direction, force vector, rope points, and timestamp. Camera positions are omitted.

Camera 2 and Camera 3 images are used for reconstruction, while all reconstructed 3D coordinates and force directions are expressed in the **Camera 1 frame**.
