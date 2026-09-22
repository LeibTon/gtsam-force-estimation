# Camera Calibration

File: `calibration_information.py`

## Reconstruction cameras

- **Image pair used for 3D reconstruction:** Camera 2 and Camera 3  
  (Camera 1 images are not used for triangulation)
- **3D coordinate frame / origin:** Camera 1  
  (all reconstructed points and rope directions are in the Camera 1 frame)
- Coordinate convention: OpenCV right-handed; extrinsic translations in mm

## What to use for image overlays

| Item | Source in file |
|------|----------------|
| Camera 2 intrinsics | `Camera2.K`, `Camera2.dist_coeffs` |
| Camera 3 intrinsics | `Camera3.K`, `Camera3.dist_coeffs` |
| Extrinsics Cam1 → Cam2 | `Stereo_1_2.R`, `Stereo_1_2.T` (`R12`, `T12`) |
| Extrinsics Cam2 → Cam3 | `Stereo_2_3.R`, `Stereo_2_3.T` (`R23`, `T23`) |

Transform convention:

```text
P_cam2 = R12 * P_cam1 + T12
P_cam3 = R23 * P_cam2 + T23
```

Because reconstructed points are in the Camera 1 frame, project into Camera 2 / Camera 3 using the extrinsics above when overlaying on annotated images.
