"""Canonical experimental parameters used by the estimation workflow.

Camera-frame transform parameters and helpers live in camera_transform.py.
"""
import numpy as np

# The physical robot hangs downward, but the experimental data are
# reoriented so its backbone points upward in the model convention.
# Consequently physical gravity maps to +B0-z.
GRAVITY_M_S2 = 9.81
GRAVITY_DIRECTION = np.array([0.0, 0.0, 1.0])
EXTERNAL_FORCE_TRANSMISSION_SCALE = 0.7365246258492919

ROUTING_ANGLE_OFFSET_RAD = np.deg2rad(85.40654)
ROD_MASS_PER_LENGTH_KG_M = 0.002 / 0.027
YOUNGS_MODULUS_PA = 71.72316664415858e9
TENDON_FRICTION_COEFFICIENT = 0.225614214748237
