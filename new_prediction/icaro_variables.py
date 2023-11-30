#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Variables useful for the bulk of ICARO project- new prediction
"""

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "ICARO: IC50 binding Affinity Regression Optimized"

"""
Folder locations
"""
SYSTEM_SEP = "/"
CSV_SEP = ","
CSV_TERMINATION = ".csv"
PARAGRAPH_SEP = "\n"
TXT_TERMINATION = ".txt"

DEFAULT_LOCATION = "/home/user/ICARO/new_prediction" #"Change/to/user/location"
DATA_FOLDER = DEFAULT_LOCATION + SYSTEM_SEP
FEATURES_FOLDER = DEFAULT_LOCATION + SYSTEM_SEP + "features"
NORMALIZE_FOLDER = FEATURES_FOLDER + SYSTEM_SEP + "normalization"
FEATURES_MORDRED_FOLDER = FEATURES_FOLDER + SYSTEM_SEP + "results_mordred"
MODEL_FOLDER = DEFAULT_LOCATION + SYSTEM_SEP + "model"
MORDRED_FOLDER = NORMALIZE_FOLDER + SYSTEM_SEP + "mordred_norm/"

RANDOM_SEED = 42
