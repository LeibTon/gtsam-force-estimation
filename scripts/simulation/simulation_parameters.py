"""Canonical experimental parameters used by the estimation workflow.

Camera-frame transform parameters and helpers live in camera_transform.py.
"""
import numpy as np

# The physical robot hangs downward, but the experimental data are
# reoriented so its backbone points upward in the model convention.
# Consequently physical gravity maps to +B0-z.
GRAVITY_M_S2 = 0.0
GRAVITY_DIRECTION = np.array([0.0, 0.0, 1.0])

ROUTING_ANGLE_OFFSET_RAD = np.deg2rad(85.40654)
ROD_MASS_PER_LENGTH_KG_M = 0.002 / 0.027
YOUNGS_MODULUS_PA = 71.72316664415858e9
SHEAR_MODULUS_PA = YOUNGS_MODULUS_PA / (2.0 * (1.0 + 0.3))
TENDON_FRICTION_COEFFICIENT = 0.0

TENSION_STD = 1e-3
POS_MES_STD = 1e-5 # In presence of noise this should be 1e-3.
POSITION_NOISE_STD_M = 1e-3