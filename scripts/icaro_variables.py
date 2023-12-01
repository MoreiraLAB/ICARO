#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Variables useful for the bulk of ICARO project
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

DEFAULT_LOCATION = "Change/to/user/location/ICARO"
DATA_FOLDER = DEFAULT_LOCATION + SYSTEM_SEP + "data"
FEATURES_FOLDER = DEFAULT_LOCATION + SYSTEM_SEP + "features"
FEATURES_MORDRED_FOLDER = FEATURES_FOLDER + SYSTEM_SEP + "mordred"
FEATURES_DAVIS_MORDRED_FOLDER = FEATURES_FOLDER + SYSTEM_SEP + "davis_mordred"
SPLITS_FOLDER = DEFAULT_LOCATION + SYSTEM_SEP + "splits"
PLOTS_FOLDER = DEFAULT_LOCATION + SYSTEM_SEP + "plots"
PERFORMANCE_FOLDER = DEFAULT_LOCATION + SYSTEM_SEP + "performance"
RESULTS_FOLDER = DEFAULT_LOCATION + SYSTEM_SEP + "results"
SCRIPTS_FOLDER = DEFAULT_LOCATION + SYSTEM_SEP + "scripts"
SUMMARY_RESULTS_FOLDER = DEFAULT_LOCATION + SYSTEM_SEP + "summary_results"
SUPPORT_FOLDER = DATA_FOLDER + SYSTEM_SEP + "support"
CLUSTERFIND_FOLDER = DEFAULT_LOCATION + SYSTEM_SEP + "clusterfind"
DISTANCES_FOLDER = DEFAULT_LOCATION + SYSTEM_SEP + "distances"
CLUSTERING_PICKLE_FILE = SUPPORT_FOLDER + SYSTEM_SEP + "clustering_dictionary.pkl"
DISTANCES_FILE = RESULTS_FOLDER + SYSTEM_SEP + "distances.csv"

"""
Useful files
"""
ICARO_BASE_FILE = SUPPORT_FOLDER + SYSTEM_SEP + "IC50_values&statistics.csv"
MORDRED_FEATURES_FILE = FEATURES_FOLDER  + SYSTEM_SEP + "icaro_mordred_features.csv"
MORDRED_FEATURES_FILE_H5 = FEATURES_FOLDER  + SYSTEM_SEP + "icaro_mordred_features.h5"
MORDRED_DAVIS_FEATURES_FILE_H5 = FEATURES_FOLDER  + SYSTEM_SEP + "icaro_davis_mordred_features.h5"
MORDRED_NORMALIZED_FEATURES_FILE_H5 = FEATURES_FOLDER  + SYSTEM_SEP + "icaro_mordred_normalized_features.h5"
MORDRED_DAVIS_NORMALIZED_FEATURES_FILE_H5 = FEATURES_FOLDER  + SYSTEM_SEP + "icaro_mordred_normalized_features_davis.h5"
TRAIN_SPLITS_FILE = RESULTS_FOLDER + SYSTEM_SEP + "train_ligands.txt"
CLUSTERFIND_LOG_FILE = CLUSTERFIND_FOLDER + SYSTEM_SEP + "clusterfind_log.csv"
MORDRED_NORMALIZATION_FILE = DEFAULT_LOCATION + SYSTEM_SEP + "new_norm.csv"
MOLECULE_INFO_FILE = FEATURES_FOLDER + SYSTEM_SEP + "molecule_target_info.csv"
DRUGBANK_LABELLED_FILE = FEATURES_FOLDER + SYSTEM_SEP + "molecule_name_drugbank.csv"

"""
Numeric Variables
"""
RANDOM_SEED = 42
