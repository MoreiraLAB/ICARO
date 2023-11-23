#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DESCIPTION: Train ERT model based on ifeature, mordred, fp-admet and protein family information.

INPUT FILES: (path: "./features/") "icaro_normalized_ifeature.h5", "icaro_mordred_normalized_features_davis.h5", "icaro_normalized_prt_family.h5", "icaro_normalized_prt_family.h5" 
OUTPUT FILES: 
- path: "./results/ERT_model/" where models of Grid Search are saved.
- path: "./results/ERT_model/metrics/" where training, test and validation metrics are saved.
"""

__author__ = "A.T. Gaspar"
__email__ = "ateresa.work@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "ICARO: IC50 binding Affinity Regression Optimized"

import h5py
import pandas as pd
import icaro_variables
import sys
import os
import numpy as np
import random
import pickle
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error
import time

start_time = time.time()
random.seed(750)
np.random.seed(750)

def XTreesBoost(x_train, y_train, x_test, x_val_prot, x_val_lig, x_val_prot_lig, parameters, file, mode=""):
    from sklearn.ensemble import ExtraTreesRegressor
    if mode=="GridSearch":
        model = ExtraTreesRegressor()
        clf = GridSearchCV(model, parameters, scoring="neg_mean_squared_error", cv=3, verbose=1)
        model = clf.fit(x_train, y_train)
        pickle.dump(model.best_params_, open(PICKLES_FOLDER+"/GS_parameters/XTrees_GSparameters.sav", 'wb'))
        file.write("XTrees:"+ str(model.best_params_)+"\n")
        y_pred_train = model.predict(x_train)
        y_pred_test = model.predict(x_test)
        return y_pred_train, y_pred_test
    else:
        parameters = pickle.load(open(PICKLES_FOLDER+"/GS_parameters/XTrees_GSparameters.sav", "rb"))
        model = ExtraTreesRegressor(**parameters).fit(x_train, y_train)
        #model = ExtraTreesRegressor(n_jobs=-1).fit(x_train, y_train) #default model
        pickle.dump(model, open(PICKLES_FOLDER+"models/XTrees.sav", 'wb'))
        y_pred_train = model.predict(x_train)
        y_pred_test = model.predict(x_test)
        y_pred_val_prot = model.predict(x_val_prot)
        y_pred_val_lig = model.predict(x_val_lig)
        y_pred_val_prot_lig = model.predict(x_val_prot_lig)
        return y_pred_train, y_pred_test, y_pred_val_prot, y_pred_val_lig, y_pred_val_prot_lig


def metrics(y_true, y_pred, file_name):
    from sklearn.metrics import r2_score
    from sklearn.metrics import mean_absolute_error
    from sklearn.metrics import median_absolute_error #chosen because is robust to outliers
    from sklearn.metrics import explained_variance_score
    from lifelines.utils import concordance_index
    from scipy import stats
    from csv import writer

    """Evaluation metrics' calculation"""
    r2 = r2_score(y_true, y_pred)
    mean_abs_error = mean_absolute_error(y_true, y_pred)
    mean_sqr_error = mean_squared_error(y_true, y_pred, squared=False)
    median_abs_error = median_absolute_error(y_true, y_pred)
    explained_var = explained_variance_score(y_true, y_pred)
    spearman = stats.spearmanr(y_true, y_pred)[0]
    pearson_coef = stats.pearsonr(y_true, y_pred)[0]
    concordance_index = concordance_index(y_true, y_pred)

    """Write results in file"""
    file = open(METRIC_FOLDER+file_name+".csv", "w")
    writer = writer(file, delimiter=";")
    writer.writerow(["r2_score", "mean_absolute_error","RMSE","median_absolute_error",
                    "explained_variance_score","spearman", "pearson", "CI"])
    writer.writerow([r2, mean_abs_error, mean_sqr_error, median_abs_error,
                    explained_var, spearman, pearson_coef, concordance_index])

def retrieve_features(pd_class):
    output = []
    for index, row in pd_class.iterrows():
        features_row = np.array([])
        lig = row["Molecule_ChEMBLID"]
        uniprot = row[0] #0 is the name of the column
        with h5py.File(FEATURES_FOLDER+"icaro_mordred_normalized_features_davis.h5") as mordred:
            features_row = np.concatenate((features_row, mordred[lig][:]))
        with h5py.File(FEATURES_FOLDER+"icaro_normalized_admet.h5") as admet:
            features_row = np.concatenate((features_row, admet[lig][:]))
        with h5py.File(FEATURES_FOLDER+"icaro_normalized_ifeature.h5") as ifeature:
            features_row = np.concatenate((features_row, ifeature[uniprot][:]))
        with h5py.File(FEATURES_FOLDER+"icaro_normalized_prt_family.h5") as family:
            features_row = np.concatenate((features_row, family[uniprot][:]))
        output.append(features_row)
    return output

FEATURES_FOLDER = icaro_variables.FEATURES_FOLDER +"/"
PICKLES_FOLDER = icaro_variables.RESULTS_FOLDER + "/ERT_model/"
METRIC_FOLDER = icaro_variables.RESULTS_FOLDER + "/ERT_model/metrics/"
RESULTS_FOLDER = icaro_variables.RESULTS_FOLDER +"/"
DATA_FOLDER = icaro_variables.DATA_FOLDER +"/"

if not os.path.exists(PICKLES_FOLDER):
    os.makedirs(PICKLES_FOLDER)
if not os.path.exists(METRIC_FOLDER):
    os.makedirs(METRIC_FOLDER)

if __name__== "__main__":
    print("#### ICARO- ERT model Training")
    mapping = pd.read_csv(DATA_FOLDER+"chembl_uniprot_mapping.txt", sep="\t", header=None, skiprows=1, usecols=[0,1])
    mapping = mapping.drop(labels=[9118,9120], axis=0, inplace=False, errors='raise')

    """Open file with train, test and validations interactions"""
    inter_train = pd.read_csv(RESULTS_FOLDER+"train_interactions.txt", sep=";", header=0)
    inter_test = pd.read_csv(RESULTS_FOLDER+"test_interactions.txt", sep=";", header=0)
    inter_validation_prot = pd.read_csv(RESULTS_FOLDER+"validation_inter_prots.txt", sep=";", header=0)
    inter_validation_lig = pd.read_csv(RESULTS_FOLDER+"validation_inter_ligs.txt", sep=";", header=0)
    inter_validation_prot_lig = pd.read_csv(RESULTS_FOLDER+"validation_inter_prots_ligs.txt", sep=";", header=0)
    out_file = open(RESULTS_FOLDER+"XTrees.txt", "w")

    inter_train = inter_train[["Molecule_ChEMBLID", "Target_ChEMBLID", "pChEMBL_median"]]\
                 .merge(mapping, how="left", left_on="Target_ChEMBLID", right_on=1).drop(labels=1, axis=1)
    inter_test = inter_test[["Molecule_ChEMBLID", "Target_ChEMBLID", "pChEMBL_median"]]\
                 .merge(mapping, how="left", left_on="Target_ChEMBLID", right_on=1).drop(labels=1, axis=1)
    inter_validation_prot = inter_validation_prot[["Molecule_ChEMBLID", "Target_ChEMBLID", "pChEMBL_median"]]\
                 .merge(mapping, how="left", left_on="Target_ChEMBLID", right_on=1).drop(labels=1, axis=1)
    inter_validation_lig = inter_validation_lig[["Molecule_ChEMBLID", "Target_ChEMBLID", "pChEMBL_median"]]\
                 .merge(mapping, how="left", left_on="Target_ChEMBLID", right_on=1).drop(labels=1, axis=1)
    inter_validation_prot_lig = inter_validation_prot_lig[["Molecule_ChEMBLID", "Target_ChEMBLID", "pChEMBL_median"]]\
                 .merge(mapping, how="left", left_on="Target_ChEMBLID", right_on=1).drop(labels=1, axis=1)

    y_train = inter_train["pChEMBL_median"].to_numpy()
    y_test = inter_test["pChEMBL_median"].to_numpy()
    y_val_prot = inter_validation_prot["pChEMBL_median"].to_numpy()
    y_val_lig = inter_validation_lig["pChEMBL_median"].to_numpy()
    y_val_prot_lig = inter_validation_prot_lig["pChEMBL_median"].to_numpy()

    """Retrieve features"""
    x_train = retrieve_features(inter_train)
    x_test = retrieve_features(inter_test)
    x_val_prot = retrieve_features(inter_validation_prot)
    x_val_lig = retrieve_features(inter_validation_lig)
    x_val_prot_lig = retrieve_features(inter_validation_prot_lig)
    features_time = time.time()
    out_file.write("XTrees Time to retrieve features "+str(features_time - start_time)+"\n")

    """Select 10% of the training interactions to perform grid search. The remaining 90% are used as tests"""
    grid_train = inter_train.sample(frac=0.1, random_state=750)
    grid_test = inter_train.drop(labels=grid_train.index, axis=0)

    grid_x_train = np.take(x_train, grid_train.index, axis=0)
    grid_x_test = np.take(x_train, grid_test.index, axis=0)
    grid_y_train = grid_train["pChEMBL_median"].to_numpy()
    grid_y_test = grid_test["pChEMBL_median"].to_numpy()

    """Perform Grid Search"""
    params = {
    "n_estimators": [50, 100, 500], #default= 100 #best= 500
    'max_features' : [1, 50, 500], #default= 1 #best= 500
    'max_depth' : [5, 7, None], #default= None #best= None
    'random_state': [750], #default= None
    'min_samples_split': [2, 3, 5], #default= 2 #best= 3
    'bootstrap': [True, False], #default= False #best= False
    'n_jobs':[-1]} #default= None

    ypred_train, ypred_test = XTreesBoost(grid_x_train, grid_y_train, grid_x_test, [], [], [], params, out_file, mode="GridSearch")
    metrics(grid_y_train, ypred_train, "GS_XTrees_train")
    metrics(grid_y_test, ypred_test, "GS_XTrees_test")
    gridsearch_time = time.time()
    out_file.write("GridSearch XTrees Time of Execution:"+str(gridsearch_time - features_time)+"\n")

    """Train, test and validate the model"""
    ypred_train, ypred_test, ypred_val_prot, ypred_val_lig, ypred_val_prot_lig = XTreesBoost(x_train, y_train, x_test, x_val_prot, x_val_lig, x_val_prot_lig, [], out_file)
    #print(ypred_test)
    metrics(y_train, ypred_train, "Train_XTrees_final")
    metrics(y_test, ypred_test, "Test_XTrees_final")
    metrics(y_val_prot, ypred_val_prot, "Validation_Proteins_XTrees_final")
    metrics(y_val_lig, ypred_val_lig, "Validation_Ligands_XTrees_final")
    metrics(y_val_prot_lig, ypred_val_prot_lig, "Validation_Proteins_Ligands_XTrees_final")
    out_file.write("Final Model XTrees Time of Execution:" + str(time.time() - gridsearch_time)+"\n")
