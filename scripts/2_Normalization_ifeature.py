#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DESCIPTION: Normalize ifeature descriptors, calculating mean and standard deviation on training data and normalize in all datasets.

NOTE: When running this script check if there is no file named iFeature_all_variances.txt in your features directory

INPUT FILES: (path: "./results/) "train_targets.txt", "test_targets.txt" and "validation_unique_targets.txt"
OUTPUT FILES: (path: "./features/) "train_targets.txt", "iFeature_header.txt" and "icaro_normalized_ifeature.h5"
"""

__author__ = "A.T. Gaspar"
__email__ = "ateresa.work@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "ICARO: IC50 binding Affinity Regression Optimized"

import matplotlib
matplotlib.use('Agg')
import h5py
import icaro_variables
import pandas as pd
import csv as txt
import os
import sys
import numpy as np
import pickle

def select_lines(file, train_set, test_set, validation_set): #train_table, test_table, validation_table
    name = file.split(".")[0]
    file_pd = pd.read_csv(file, header=0, sep="\t")
    file_pd.columns = ["id"]+[name+"_"+str(i) for i in range(1,len(file_pd.columns))]
    """Update train_table"""
    train_table = file_pd[file_pd.iloc[:,0].isin(train_set)==True].reset_index(drop=True)
    """Update test_table"""
    test_table = file_pd[file_pd.iloc[:,0].isin(test_set)==True].reset_index(drop=True)
    """Update validation_table"""
    validation_table = file_pd[file_pd.iloc[:,0].isin(validation_set)==True].reset_index(drop=True)
    return train_table, test_table, validation_table

def process_variance(train, test, validation, counter=0):
    """This function checks drops columns without variance, given a pandas table.
    It returns the train, test and validation with the feature set that have variance"""
    if "iFeature_all_variances.txt" in os.listdir(FEATURES_FOLDER):
        file = open(FEATURES_FOLDER+"iFeature_all_variances.txt","a")
    else:
        file = open(FEATURES_FOLDER+"iFeature_all_variances.txt","w")
    name = train.columns[1][:-2]
    size = train.shape[1]-1
    for column in train.columns[1:]:
        train[column] = pd.to_numeric(train[column], errors="coerce")
        file.write(column+","+str(train[column].var())+"\n")
        if train[column].var() <= 0.01:
            train = train.drop(columns = column)
            test = test.drop(columns = column)
            validation = validation.drop(columns = column)
            counter+=1
    print(name+"-Number of columns dropped: "+str(counter)+" in "+str(size)+" features.")
    return train, test, validation

def standardization(train, test, validation, file, path):
    from sklearn.preprocessing import StandardScaler
    """Construct scaler and save normalization parameters in pickle"""
    scaler = StandardScaler().fit(train)
    #print(scaler.mean_, scaler.var_**(1/2))
    pickle.dump(scaler, open(path+file+'_parameters_training_set.sav', 'wb'))
    """Apply normalization to train, test and validation data"""
    train_transformed = scaler.transform(train)
    test_transformed = scaler.transform(test)
    validation_transformed = scaler.transform(validation)
    return train_transformed, test_transformed, validation_transformed


def retrieve_features(files, pca_files, train_ids, test_ids, validation_ids, reduction=""):
    for file in files:
        name = file.split(".")[0]
        train_lines, test_lines, validation_lines = select_lines(file, train_ids, test_ids, validation_ids)
        """Remove variance and normalize the train, test and validation sets"""
        if file in pca_files:
            train_lines, test_lines, validation_lines = process_variance(train_lines, test_lines, validation_lines)
        if file == files[0]:
            train_data, test_data, validation_data = train_lines.iloc[:,0],test_lines.iloc[:,0],validation_lines.iloc[:,0]
        if train_lines.shape[1] > 1:
            train_standard, test_standard, validation_standard = standardization(train_lines.iloc[:,1:].to_numpy(),
                                    test_lines.iloc[:,1:].to_numpy(), validation_lines.iloc[:,1:].to_numpy(), name, PICKLES_FOLDER)
        else:
            continue

        train_lines = pd.DataFrame(train_standard,
                      columns=[name+"_"+str(i) for i in range(1, np.shape(train_standard)[1]+1)])
        test_lines = pd.DataFrame(test_standard,
                     columns=[name+"_"+str(i) for i in range(1, np.shape(test_standard)[1]+1)])
        validation_lines = pd.DataFrame(validation_standard,
                           columns=[name+"_"+str(i) for i in range(1, np.shape(validation_standard)[1]+1)])
        train_data = pd.concat((train_data,train_lines), axis=1)
        test_data = pd.concat((test_data, test_lines), axis=1)
        validation_data = pd.concat((validation_data, validation_lines), axis=1)
    return train_data, test_data, validation_data

def write_h5(h5_file, dataframe):#, mode="train"):
    """Given a data frame and a mode creates a group in a h5 file, where a key is one line of the data frame"""
    #group = h5_file.create_group(mode)
    for index, row in dataframe.iterrows():
        if row.iloc[0] not in h5_file.keys():
            h5_file.create_dataset(name= row.iloc[0],
                                 data= row.iloc[1:].to_numpy(dtype="float32"))
    return dataframe.columns[1:]

PICKLES_FOLDER = icaro_variables.FEATURES_FOLDER + "/ifeature_normalization/"
RESULTS_FOLDER = icaro_variables.RESULTS_FOLDER+"/"
FEATURES_FOLDER = icaro_variables.FEATURES_FOLDER+"/"

if not os.path.exists(PICKLES_FOLDER):
    os.makedirs(PICKLES_FOLDER)

if not os.path.exists(FEATURES_FOLDER+"results_iFeature"):
    os.makedirs(FEATURES_FOLDER+"results_iFeature")

if __name__== "__main__":
    """Extract IDs of each set"""
    print("#### ICARO- iFeature Normalization")
    train_ids = pd.read_csv(RESULTS_FOLDER+"/train_targets.txt", sep=";", header=None)
    test_ids = pd.read_csv(RESULTS_FOLDER+"/test_targets.txt", sep=";", header=None)
    validation_ids = pd.read_csv(RESULTS_FOLDER+"/validation_unique_targets.txt", sep=";", header=None)
    train_ids, test_ids, validation_ids = train_ids.iloc[:,1].sort_values(), test_ids.iloc[:,1].sort_values(),\
                                          validation_ids.iloc[:,1].sort_values()

    """Join iFeature features by train, test and validation sets, iterating through the files in the
    iFeature results directory"""
    sep = "/"
    directory = FEATURES_FOLDER+"results_iFeature"+sep
    os.chdir(directory)
    dir_files = os.listdir(directory)
    dir_files.sort()
    for_pca = []
    to_remove = []
    for f in dir_files:
        if f.endswith("_reduced.tsv"):
            for_pca.append(f.replace("_reduced",""))
            to_remove.append(f)
        elif f.endswith(".tsv.png"):
            to_remove.append(f)
    dir_files = [i for i in dir_files if i not in to_remove]
    print("Nº of features type:", len(dir_files))
    print("Nº of features type for PCA",len(for_pca), for_pca)
    
    train_dataframe, test_dataframe, validation_dataframe = retrieve_features(dir_files, for_pca, train_ids, test_ids, validation_ids)

    # """Write features in a csv - optional"""
    # train_dataframe.to_csv(path_file.results+"normalized_train_ifeature.csv", mode="w", sep=",", header=True, index=False)
    # test_dataframe.to_csv(path_file.results+"normalized_test_ifeature.csv", mode="w", sep=",", header=True, index=False)
    # validation_dataframe.to_csv(path_file.results+"normalized_validation_ifeature.csv", mode="w", sep=",", header=True, index=False)
    print("Arrays Shape (including ID column):", train_dataframe.shape, test_dataframe.shape, validation_dataframe.shape)

    """Write features, by ChEMBL ID, in a h5 file, organised by train (80%), test (20%)
    and validation (10%) sets"""
    h5_ifeature = h5py.File(FEATURES_FOLDER+"icaro_normalized_ifeature.h5", "w")
    write_h5(h5_ifeature, train_dataframe)
    write_h5(h5_ifeature, test_dataframe)#, mode="test")
    header = write_h5(h5_ifeature, validation_dataframe)#, mode="validation")

    """Write iFeature features name"""
    with open(FEATURES_FOLDER+"/iFeature_header.txt", "w") as output_file:
        writer = txt.writer(output_file, delimiter = ",")
        output_file.write("Number of Features:"+str(len(header))+"\n")
        writer.writerow(header)
