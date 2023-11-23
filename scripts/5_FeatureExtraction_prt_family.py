#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DESCIPTION: Get protein family information on dataset proteins.

INPUT FILES: (path: "./data/support/") "unique_prt_ids.csv, "targets_family.csv"
OUTPUT FILES: (path: "./features/) "icaro_normalized_prt_family.h5"
"""

__author__ = "A.T. Gaspar"
__email__ = "ateresa.work@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "ICARO: IC50 binding Affinity Regression Optimized"

import pandas as pd
import icaro_variables
import h5py
from csv import writer
import sys

def write_h5(h5_file, dataframe):#, mode="train"):
    """Given a data frame and a mode creates a group in a h5 file, where a key is one line of the data frame"""
    #group = h5_file.create_group(mode)
    for index, row in dataframe.iterrows():
        if row.iloc[0] not in h5_file.keys():
            h5_file.create_dataset(name= row.iloc[0],
                                 data= row.iloc[1:].to_numpy(dtype="float32"))
    return dataframe.columns[1:]

def onehot_fam_feature(targets_data, fam_list):
    targets = targets_data["Target_ChEMBLID"].unique().tolist()
    df = pd.DataFrame(0, index=targets, columns = fam_list)
    for index, row in df.iterrows():
        one_prot_data = targets_data.loc[targets_data["Target_ChEMBLID"]==index]
        one_prot_fam = one_prot_data["Class_l1"].unique().tolist()
        row[one_prot_fam] = 1
        df.update(row)
    return df

SUPPORT_FOLDER = icaro_variables.SUPPORT_FOLDER+"/"
FEATURES_FOLDER = icaro_variables.FEATURES_FOLDER+"/"

if __name__ == "__main__":
    print("#### ICARO- Protein Family Information Features Extraction")
    mapping = pd.read_csv(SUPPORT_FOLDER+"unique_prt_ids.csv", sep=",", header=0, usecols=["Uniprot ID","ChEMBL ID"])
    mapping.index = [mapping["Uniprot ID"].values]
    mapping = mapping.drop(labels=["A0A2P2HJL8","G4WW85"], axis=0)
    fam_sheet = pd.read_csv(SUPPORT_FOLDER+"targets_family.csv", sep=",", header=0, usecols=["Target_ChEMBLID", "Class_l1"])
    fam_data = fam_sheet.drop_duplicates()
    many_class_targets = fam_data.loc[fam_data.duplicated(subset="Target_ChEMBLID",keep=False)==True] #40 targets have more than 1 family
    families = fam_data["Class_l1"].unique().tolist() #There are 15 unique families
    families.sort()

    """Write dataframe with one-hot-encoding representation of proteins' families"""
    final_data = onehot_fam_feature(fam_data, families)
    final_data = final_data.loc[:, final_data.var(axis=0, numeric_only=True) != 0]
    final_data.insert(0, "index", final_data.index)
    cols = final_data.columns[1:].tolist()
    final_data = final_data.merge(mapping, how="inner", right_on="ChEMBL ID", left_on="index")
    neworder = ["Uniprot ID"]+cols
    final_data=final_data.reindex(columns=neworder)
    header = open(FEATURES_FOLDER+"Features_family_header.txt", "w")
    header.write("Number of Features:"+str(len(neworder[1:]))+"\n")
    header_writer = writer(header, delimiter=',')
    header_writer.writerow(neworder[1:])
    output_file = h5py.File(FEATURES_FOLDER+"icaro_normalized_prt_family.h5", "w")
    write_h5(output_file, final_data)
