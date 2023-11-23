#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DESCIPTION: Normalize fp-admet descriptors, calculating mean and standard deviation on training data and normalize in all datasets.

INPUT FILES: 
- files in path: "./features/fpadmet/RESULTS/"
- (path: "./results/") "train_ligands.txt", "test_ligands.txt", "validation_unique_ligands.txt"
OUTPUT FILES: (path: "./features/) "icaro_normalized_admet.h5"
"""

__author__ = "A.T. Gaspar"
__email__ = "ateresa.work@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "ICARO: IC50 binding Affinity Regression Optimized"

import icaro_variables
import pandas as pd
import h5py
import sys
import os

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

def write_h5(h5_file, dataframe):#, mode="train"):
    """Given a data frame and a mode creates a group in a h5 file, where a key is one line of the data frame"""
    #group = h5_file.create_group(mode)
    for index, row in dataframe.iterrows():
        if row.iloc[0] not in h5_file.keys():
            h5_file.create_dataset(name= row.iloc[0],
                                 data= row.iloc[1:].to_numpy(dtype="float32"))
    return dataframe.columns[1:]

def process_admet(file_name, dict, train_set, test_set, validation_set):
    file = pd.read_csv(file_name, sep=" ", header=0)
    print(file)
    nfeature = file_name.replace(".txt","").split("_")[1]
    if nfeature in dict.keys():
        file["Predicted_"+nfeature] = file.iloc[:,0].map(dict[nfeature])
    else:
        import numpy as np
        train, test, validation = file.loc[train_set, "Predicted"].to_frame(),file.loc[test_set,"Predicted"].to_frame(), \
                                  file.loc[validation_set, "Predicted"].to_frame()
        if train.iloc[:,0].var() == 0.01:
            print(file_name, train.iloc[:,0].var())
            return None
        train_norm, test_norm, validation_norm = standardization(train.to_numpy(), test.to_numpy(), validation.to_numpy(),
                                                "admet"+nfeature, path=PICKLES_FOLDER)
        train_values = pd.DataFrame(train_norm, index = train.index.tolist(), columns=["Predicted_"+nfeature])
        test_values = pd.DataFrame(test_norm, index=test.index.tolist(), columns=["Predicted_"+nfeature])
        validation_values = pd.DataFrame(validation_norm, index=validation.index.tolist(), columns=["Predicted_"+nfeature])

        total = pd.concat((train_values, test_values), axis=0)
        unique_total = total[~total.index.duplicated(keep='first')]
        file = pd.concat((unique_total, validation_values))
        #print(file.iloc[:10,:])
    return file["Predicted_"+nfeature]

PICKLES_FOLDER = icaro_variables.FEATURES_FOLDER + "/fpadmet_normalization/"
RESULTS_FOLDER = icaro_variables.RESULTS_FOLDER+"/"
ADMET_FOLDER = icaro_variables.FEATURES_FOLDER+"/fpadmet/RESULTS/"

if __name__ == "__main__":
    """Extract IDs of each set"""
    print("#### ICARO- FP-Admet Normalization")
    train_ids = pd.read_csv(RESULTS_FOLDER+"train_ligands.txt", sep=";", header=None)
    test_ids = pd.read_csv(RESULTS_FOLDER+"test_ligands.txt", sep=";", header=None)
    validation_ids = pd.read_csv(RESULTS_FOLDER+"validation_unique_ligands.txt", sep=";", header=None)
    train_ids, test_ids, validation_ids = train_ids.iloc[:,0].sort_values(), test_ids.iloc[:,0].sort_values(),\
                                          validation_ids.iloc[:,0].sort_values()
    """Change to ADMET directory and open files. All features are binary or a class except
    from the features 50 to 56, inclusive, that are a numeric value"""
    os.chdir(ADMET_FOLDER)
    files = os.listdir()
    files.sort()
    files = [i for i in files if i.startswith("predicted")]
    binary_classes = {"1":{"Positive":1,"Negative":0}, "2":{"Yes":1,"No":0}, "3":{"High":1,"Low":2},
                      "4":{"active":1,"inactive":0}, "5":{"Stable":3,"Moderate":2,"Unstable":1},
                      "6":{"EPA1":1,"EPA2":2,"EPA3":3,"EPA4":4}, "7":{"P":1,"N":0}, "8":{"P":1,"N":0},
                      "9":{"Positive":1,"Negative":0}, "10":{"Positive":1,"Negative":0},
                      "11":{"Positive":1,"Negative":0}, "12":{"P":1,"N":0},
                      "13":{"Positive":1,"Negative":0}, "14":{"Inhibitor":1,"NonInhibitor":0},
                      "15":{"Positive":1,"Negative":0}, "16":{"Inhibitor":1,"NonInhibitor":0},
                      "17":{"Positive":1,"Negative":0}, "18":{"Positive":1,"Negative":0},
                      "19":{"Active":1,"Inactive":0}, "20":{"Yes":1,"No":0}, "21":{"Yes":1,"No":0},
                      "22":{"Yes":1,"No":0}, "23":{"yes":1,"no":0}, "24":{"positive":1,"negative":0}, "25":{"Yes":1,"No":0},
                      "26":{"P":1,"N":0}, "27":{"Positive":1,"Negative":0}, "28":{"Active":1,"Inactive":0},
                      "29":{"Carcinogen":1,"NonCarcinogen":0}, "30":{"Soluble":1,"Insoluble":0},
                      "31":{"Positive":1,"Negative":0}, "32":{"LOW":1,"MODERATE":2,"HIGH":3},
                      "33":{"Blocker":1,"NonBlocker":0}, "34":{"Inhibitor":1,"Noninhibitor":0}, "35":{"Yes":1,"No":0},
                      "36":{"Yes":1,"No":0}, "37":{"LOW":1,"MEDIUM":2,"HIGH":3}, "38":{"LOW":1,"MEDIUM":2,"HIGH":3},
                      "39":{"LOW":1,"MEDIUM":2,"HIGH":3}, "40":{"P":1,"N":0}, "41":{"Active":1,"Inactive":0},
                      "42":{"Active":1,"Inactive":0}, "43":{"Active":1,"Inactive":0},
                      "44":{"Active":1,"Inactive":0}, "45":{"Active":1,"Inactive":0}, "46":{"Active":1,"Inactive":0},
                      "47":{"Active":1,"Inactive":0}, "48":{"Active":1,"Inactive":0}, "49":{"Active":1,"Inactive":0},
                      "57":{"inhibitor":1,"noninhibitor":0}, "58":{"Active":1,"Inactive":0}
                      }
    h5_admet = h5py.File(icaro_variables.FEATURES_FOLDER+"/icaro_normalized_admet.h5", "w")
    
    final_table = process_admet(files[0], binary_classes, train_ids, test_ids, validation_ids).to_frame()
    for f in files[1:]:
        output_feat = process_admet(f, binary_classes, train_ids, test_ids, validation_ids)
        if output_feat.isnull().values.any() == True:
            print(f, output_feat.isnull().values.any())
        if output_feat is None:
            continue
        if output_feat.name in final_table.columns:
            final_table[output_feat.name].update(output_feat)
        else:
            final_table = pd.concat([final_table, output_feat], axis=1, join='outer')
        #print(final_table.iloc[:10,:])
        if final_table.isnull().values.any() == True:
            print(f, final_table)
    write_h5(h5_admet, final_table)
