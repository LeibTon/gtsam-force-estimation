"""
Calibration between the UIUC-GT experiment's camera frame (Cam1-origin, millimetres)
and the Cosserat-rod estimator's spatial/model frame (robot base, metres).

The transform is fit once from the zero-tendon-tension trial, where the measured
backbone shape should correspond to the model's unloaded straight-rod shape.
Given matched, ordered point sets (base -> tip, 11 points each), we solve for the
rigid rotation + translation (Kabsch/Procrustes) that best maps camera-frame points
(after mm -> m unit conversion) onto the model-frame points. Scale is fixed at
1/1000 (mm -> m); only rotation and translation are fit.
"""
import json

import numpy as np

# Camera-1 reconstruction frame C1 -> robot base/model frame B0 = W.
# Fitted from the five full-slack centreline observations: Disc 0 is the B0
# origin and the best-fit full-slack axis is +B0-z.  The otherwise
# unobservable roll is retained nearest to the supplied C1-to-B0 convention.
CAMERA_TO_MODEL_R = np.array([[0.9997654447560451, -0.0102739807593613, -0.0190656967117271],
                              [-0.0194970026292784, -0.0436577215860512, -0.9988562810706998],
                              [0.0094298653342028, 0.9989937180309350, -0.0438477932683925]])
CAMERA_TO_MODEL_T_M = np.array([0.0277030067913273, 0.4548766378382083, 0.1465360181040651])

# Original supplied C1 -> B0 axis convention.  Use this for annotated force
# directions; translations never apply to a direction vector.
ORIGINAL_CAMERA_TO_MODEL_R = np.array([[1.0, 0.0, 0.0],
                                       [0.0, 0.0, -1.0],
                                       [0.0, 1.0, 0.0]])


def load_arm_points_mm(json_path):
    """Load the 11x3 points_3d array (millimetres) from an arm_points.json file."""
    with open(json_path, 'r') as f:
        data = json.load(f)
    return np.array(data['points_3d'], dtype=float)


def camera_to_model_transform():
    return CAMERA_TO_MODEL_R.copy(), CAMERA_TO_MODEL_T_M.copy()


def force_direction_camera_to_model_rotation():
    # Rope points and backbone points are reconstructed by the same camera
    # pair.  A direction is transformed by the rotational part of the exact
    # same C1 -> B0 calibration used for positions (translation is omitted).
    return CAMERA_TO_MODEL_R.copy()


def apply_transform(camera_points_mm, R, t):
    """Map camera-frame points (mm) into the model frame (m) using a fitted (R, t)."""
    camera_points_m = camera_points_mm / 1000.0
    return (R @ camera_points_m.T).T + t


def apply_transform_direction(direction_camera, R):
    """
    Map a direction vector (not a point) from the camera frame into the model
    frame: rotation only, no translation, no mm->m scale (directions are unitless).
    """
    direction_camera = np.asarray(direction_camera, dtype=float)
    return R @ direction_camera
