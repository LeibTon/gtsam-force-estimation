"""
Complete Camera Calibration Parameters
All Three Cameras + Stereo Parameters
Generated: 05-Dec-2025 20:07:16

=============================================================================
RECONSTRUCTION SETUP (IMPORTANT)
=============================================================================
- Image pair used for 3D reconstruction: Camera 2 and Camera 3
  (Camera 1 images are NOT used for triangulation).
- Reconstructed 3D points / rope directions are expressed in the Camera 1 frame
  (Camera 1 origin).

- For overlaying estimates on annotated images, use Camera2 / Camera3 intrinsics
  and the extrinsics that relate them to the Camera 1 world frame
  (see Stereo_1_2 and Stereo_2_3).

=============================================================================
COORDINATE / TRANSFORM CONVENTIONS
=============================================================================
1. COORDINATE SYSTEM:
   OpenCV Standard (Right-Handed): X-Right, Y-Down, Z-Forward.
   Translation units: millimeters (mm).

2. TRANSFORM NOTATION (R_AB, T_AB):
   - "Stereo_A_B" transforms a point FROM Camera A TO Camera B.
   - Mathematical relation: P_B = R_AB * P_A + T_AB

   - T_AB: Translation from Camera A origin to Camera B origin.
   - R_AB: Rotation of Camera A relative to Camera B.

3. INVERSE TRANSFORM (B -> A):
   - P_A = R_AB.T * (P_B - T_AB)
=============================================================================
"""

import numpy as np

# ========================================
# CAMERA 1 PARAMETERS (Left Camera)
# ========================================
class Camera1:
    """Intrinsics for Camera (Right/New)."""
    # Focal length [px]
    fx = 1403.6137
    fy = 1407.1332
    focal_length = np.array([fx, fy])
    
    # Principal point [px]
    cx = 951.2505
    cy = 536.1726
    principal_point = np.array([cx, cy])
    
    # Arm start point in image [px] 
    # (Note: This value was not present in the image, set to 0,0 as placeholder)
    # Arm start point in image [px]
    arm_start_point = np.array([949, 130])
    
    # Intrinsic Matrix K
    K = np.array([
        [1403.6137, 0.0, 951.2505],
        [0.0, 1407.1332, 536.1726],
        [0.0, 0.0, 1.0]
    ])
    
    # Distortion Coefficients (k1, k2, p1, p2, k3)
    # Image lists k3 as "unused", so it is set to 0.0
    dist_coeffs = np.array([-0.053873, 0.117899, 0.0, 0.0, 0.0])
    
    # Image size (H, W)
    image_size = (1080, 1920)
    
    # Calibration Error
    mean_reprojection_error = 0.8074



# ========================================
# CAMERA 2 PARAMETERS (Middle Camera)
# ========================================
class Camera2:
    """Intrinsics for the next camera (from image_eaec01)."""
    # Focal length [px]
    fx = 1445.4101
    fy = 1444.0307
    focal_length = np.array([fx, fy])
    
    # Principal point [px]
    cx = 967.2628
    cy = 523.0121
    principal_point = np.array([cx, cy])
    
    # Arm start point in image [px]
    arm_start_point = np.array([1025, 158])
    
    # Intrinsic Matrix K
    # Standard CV format: [[fx, 0, cx], [0, fy, cy], [0, 0, 1]]
    K = np.array([
        [1445.4101, 0.0, 967.2628],
        [0.0, 1444.0307, 523.0121],
        [0.0, 0.0, 1.0]
    ])
    
    # Distortion Coefficients (k1, k2, p1, p2, k3)
    # Based on image: k1, k2 present; p1, p2 are 0.0; k3 is unused
    dist_coeffs = np.array([-0.056712, 0.447845, 0.0, 0.0, 0.0])
    
    # Image size (H, W)
    # Assuming standard 1080p based on previous context
    image_size = (1080, 1920)
    
    # Calibration Error
    mean_reprojection_error = 0.6950
# ========================================
# CAMERA 3 PARAMETERS (Right Camera)
# ========================================
class Camera3:
    """Intrinsics for the third camera (from image_eaf6cb)."""
    # Focal length [px]
    fx = 1394.0849
    fy = 1393.7046
    focal_length = np.array([fx, fy])
    
    # Principal point [px]
    cx = 925.4751
    cy = 540.2804
    principal_point = np.array([cx, cy])
    
    # Arm start point in image [px]
    arm_start_point = np.array([1028, 185])
    
    # Intrinsic Matrix K
    # Standard CV format: [[fx, 0, cx], [0, fy, cy], [0, 0, 1]]
    K = np.array([
        [1394.0849, 0.0, 925.4751],
        [0.0, 1393.7046, 540.2804],
        [0.0, 0.0, 1.0]
    ])
    
    # Distortion Coefficients (k1, k2, p1, p2, k3)
    # k3 is unused; p1 and p2 are 0.0
    dist_coeffs = np.array([-0.176872, 0.048886, 0.0, 0.0, 0.0])
    
    # Image size (H, W)
    image_size = (1080, 1920)
    
    # Calibration Error
    mean_reprojection_error = 0.6327


# ========================================
# STEREO PARAMETERS: CAM 0 -> CAM 1 (Likely Unused)
# ========================================
class Stereo_0_1:
    """Parameters for transform from Cam0 to Cam1 (R01, T01)."""
    baseline = 460.31
    T = np.array([-460.2657, 29.5538, -9.5489])
    R = np.array([
        [0.7829 , -0.0039 , -0.6221],
        [-0.0040 , 0.9999 , -0.0113],
        [0.6221 , 0.0113 , 0.7828]
    ])
    rotation_x = 0.83
    rotation_y = -38.47
    rotation_z = 0.29
    mean_reprojection_error = 0.9532


# ========================================
# STEREO PARAMETERS: CAM 1 -> CAM 2 (CRITICAL)
# ========================================

class Stereo_1_2:
    """
    TRANSFORMATION: Camera 1 (Left) -> Camera 2 (Middle)
    
    Mathematics:
    P_cam2 = R12 * P_cam1 + T12
    """
    # Baseline [mm]
    baseline = 385.65
    
    # ---------------------------------------------------------
    # T12: Translation from Cam1 to Cam2
    # Vector points from Cam1 origin to Cam2 origin.
    # ---------------------------------------------------------
    T = np.array( [-274.4366, -26.2725, 269.6734])  # Unit: mm
    
    # ---------------------------------------------------------
    # R12: Rotation from Cam1 to Cam2
    # ---------------------------------------------------------
    R = np.array([
        [0.7289 , 0.0136 , -0.6845],
        [-0.0013 , 0.9998 , 0.0185],
        [0.6846 , -0.0126 , 0.7288]
    ])
    
    # Euler angles (for reference)
    rotation_x = -1.46
    rotation_y = -43.19  # Significant Y-rotation (approx 45 deg)
    rotation_z = -1.07
    
    mean_reprojection_error = 0.8769  


# ========================================
# STEREO PARAMETERS: CAM 2 -> CAM 3 (CRITICAL)
# ========================================

class Stereo_2_3:
    """
    TRANSFORMATION: Camera 2 (Middle) -> Camera 3 (Right)
    
    Mathematics:
    P_cam3 = R23 * P_cam2 + T23
    """
    # Baseline [mm]
    baseline = 443.05
    
    # ---------------------------------------------------------
    # T23: Translation from Cam2 to Cam3
    # Vector points from Cam2 origin to Cam3 origin.
    # ---------------------------------------------------------
    T = np.array( [-436.3091, 44.5657, 62.8006])  # Unit: mm
    
    # ---------------------------------------------------------
    # R23: Rotation from Cam2 to Cam3
    # ---------------------------------------------------------
    R = np.array([
        [0.6949  ,-0.0152  ,-0.7189],
        [0.0157  ,0.9999  ,-0.0060],
        [0.7189  ,-0.0071  ,0.6951]
    ])
    
    # Euler angles (for reference)
    rotation_x = 0.49
    rotation_y =  -45.96 # Significant Y-rotation (approx 49 deg)
    rotation_z = 1.25
    
    mean_reprojection_error = 1.1993


# ========================================
# STEREO PARAMETERS: CAM 3 -> CAM 4 (Likely Unused)
# ========================================
class Stereo_3_4:
    baseline = 477.45
    T = np.array( [-346.7616, -10.9668, 328.0207])
    R = np.array([
        [0.7295,  0.0164 , -0.6838],
        [-0.0039 , 0.9998 , 0.0198],
        [0.6839 , -0.0117 , 0.7294]
    ])
    rotation_x = -1.55
    rotation_y =  -43.14
    rotation_z = -1.29
    mean_reprojection_error = 0.8353


# ========================================
# ALIASES FOR CODE READABILITY
# ========================================

# Camera Matrices (Intrinsics)
K1 = Camera1.K
K2 = Camera2.K
K3 = Camera3.K

# Distortion Coefficients
dist1 = Camera1.dist_coeffs
dist2 = Camera2.dist_coeffs
dist3 = Camera3.dist_coeffs

# Stereo Extrinsics (Explicit Naming for Logic)
# --------------------------------------------------
# R12: Rotation matrix FROM Camera 1 TO Camera 2
# T12: Translation vector FROM Camera 1 TO Camera 2
# Usage: P2 = R12 @ P1 + T12
# --------------------------------------------------
R12 = Stereo_1_2.R
T12 = Stereo_1_2.T

# --------------------------------------------------
# R23: Rotation matrix FROM Camera 2 TO Camera 3
# T23: Translation vector FROM Camera 2 TO Camera 3
# Usage: P3 = R23 @ P2 + T23
# --------------------------------------------------
R23 = Stereo_2_3.R
T23 = Stereo_2_3.T


# ========================================
# HELPER FUNCTIONS
# ========================================

def get_camera_params(camera_id):
    """Get (K, dist, size) by ID."""
    cameras = {
        1: (Camera1.K, Camera1.dist_coeffs, Camera1.image_size),
        2: (Camera2.K, Camera2.dist_coeffs, Camera2.image_size),
        3: (Camera3.K, Camera3.dist_coeffs, Camera3.image_size)
    }
    return cameras.get(camera_id, None)

def get_stereo_params(pair):
    """Get (R, T, baseline) by pair string (e.g. '1-2')."""
    pairs = {
        "1-2": (Stereo_1_2.R, Stereo_1_2.T, Stereo_1_2.baseline),
        "2-3": (Stereo_2_3.R, Stereo_2_3.T, Stereo_2_3.baseline)
    }
    return pairs.get(pair, None)

def print_summary():
    """Debug print for verification."""
    print("=" * 60)
    print("CALIBRATION PARAMETERS SUMMARY (Standard: OpenCV Right-Handed)")
    print("=" * 60)
    print(f"R12 (Cam1->Cam2) Loaded. T12 Z-shift: {T12[2]:.2f} mm")
    print(f"R23 (Cam2->Cam3) Loaded. T23 Z-shift: {T23[2]:.2f} mm")
    print("=" * 60)

if __name__ == "__main__":
    print_summary()