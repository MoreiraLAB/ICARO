#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Split Data into training, test and PPI, PID and PIT validation sets
"""

__author__ = "A.T. Gaspar"
__email__ = "ateresa.work@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "ICARO: IC50 binding Affinity Regression Optimized"

#INPUT FILES: (path: "./data/") "chembl_uniprot_mapping.txt" and (path: "./support/") "IC50_values&statistics_final_dataset.csv"
#OUTPUT FILES: (path: "./results/") 
#"train_targets.txt", "train_ligands.txt", train_interactions.txt"- files with train protein ids, ligand ids and corresponding interactions
#"test_targets.txt", "test_ligands.txt", "test_interactions.txt"- files with test protein ids, ligand ids and corresponding interactions
#"validation_unique_targets.txt", "validation_unique_ligands", "validation_all_interactions.txt"- files with validation protein ids, ligand ids and corresponding interactions
#"validation_inter_prots.txt", "validation_inter_ligs.txt", "validation_inter_prots_ligs"- files with interactions on the validation protein dataset, ligand dataset and protein-ligand dataset

import pandas as pd
import numpy as np
from csv import writer
import icaro_variables
import random
import os

random.seed(750)
np.random.seed(750)

def data_selection(samples, percentage = 0.1):
    """This function selects 10% of samples from a list.
    """
    random.shuffle(samples)
    samples_to_extract = samples[:int(percentage*len(samples))]
    return samples_to_extract

def train_test_split(input_list, train_size = 0.8):
    """This function splits interactions for train and test.
       Default is 80% for train and 20% for test."""
    random.shuffle(input_list)
    print("Samples for split: ", input_list[:20])
    train = input_list[:int(train_size*len(input_list))]
    print("Train samples", len(train))
    test = input_list[int(train_size*len(input_list)):]
    print("Test samples", len(test))
    return train, test

RESULTS_FOLDER = icaro_variables.RESULTS_FOLDER + "/"
FEATURES_FOLDER = icaro_variables.FEATURES_FOLDER + "/"

if not os.path.exists(RESULTS_FOLDER):
    os.makedirs(RESULTS_FOLDER)
if not os.path.exists(FEATURES_FOLDER):
    os.makedirs(FEATURES_FOLDER)

if __name__== "__main__":
    """Split proteins and ligands IDs into train, test and validation sets"""
    print("#### ICARO- Dataset split")
    map_file = pd.read_csv(icaro_variables.DATA_FOLDER+"/chembl_uniprot_mapping.txt", sep="\t", header=None, skiprows=1, usecols=[0,1])
    map_file = map_file.drop(labels=[9118,9120], axis=0, inplace=False, errors='raise')
    interactions_file = pd.read_csv(icaro_variables.ICARO_BASE_FILE, sep=",", header=0).iloc[:-2,:]
    print("Total Interactions: ", str(len(interactions_file)))
    interactions_file = interactions_file[interactions_file['Canonical_Smiles'].notna()]

    targets = interactions_file[["Target_ChEMBLID"]].merge(map_file, how="inner", left_on="Target_ChEMBLID", right_on=1)
    ligands = interactions_file[["Molecule_ChEMBLID","Canonical_Smiles"]].dropna(axis=0)

    """Put aside 10% of the proteins and 10% of the ligands and the respectives interactions"""
    targets_validation = data_selection(targets["Target_ChEMBLID"].unique().tolist(), percentage=0.1)
    print("Number of Proteins set apart for Validation: "+str(len(targets_validation)))
    ligands_validation = data_selection(ligands["Molecule_ChEMBLID"].unique().tolist(), percentage=0.1)
    print("Number of Ligands set apart for Validation: "+str(len(ligands_validation)))

    val_inter_by_prot = interactions_file.loc[interactions_file["Target_ChEMBLID"].isin(targets_validation)]
    val_inter_by_lig = interactions_file.loc[interactions_file["Molecule_ChEMBLID"].isin(ligands_validation)]
    val_inter = pd.concat((val_inter_by_prot, val_inter_by_lig), axis=0)
    val_prot_lig = val_inter.loc[val_inter.duplicated()==True]
    print("Predict Pairwise Interactions (PPI) dataset:", val_prot_lig.shape)
    val_inter = val_inter.drop_duplicates()
    val_inter_by_prot = val_inter_by_prot.drop(labels=val_prot_lig.index).sort_index()
    val_inter_by_lig = val_inter_by_lig.drop(labels=val_prot_lig.index).sort_index()
    print("Predict Interaction of Drugs (PID) dataset:", val_inter_by_lig.shape)
    print("Predict Interaction of Targets (PIT) dataset:", val_inter_by_prot.shape)

    """Split the remaining interactions into train and test"""
    remaining_interactions = interactions_file.drop(labels=val_inter.index).sort_index()
    inter_train = remaining_interactions.sample(frac=0.8, random_state=750).sort_index()
    inter_test = remaining_interactions.drop(labels= inter_train.index).sort_index()
    print("Train dataset:", inter_train.shape)
    print("Test dataset:",inter_test.shape)

    """Update the list of unique proteins and ligands in validations sets"""
    ligands_reunion = val_inter.merge(remaining_interactions, how="inner", on="Molecule_ChEMBLID")
    ligands_reunion = ligands_reunion["Molecule_ChEMBLID"].unique().tolist()
    print("Ligands in both validation and training/testing sets:",len(ligands_reunion))
    ligands_validation = list(set(val_inter["Molecule_ChEMBLID"].unique().tolist()) - set(ligands_reunion))
    ligands_validation.sort()

    targets_reunion = val_inter.merge(remaining_interactions, how="inner", on="Target_ChEMBLID")
    targets_reunion = targets_reunion["Target_ChEMBLID"].unique().tolist()
    print("Targets in both validation and training/testing sets:",len(targets_reunion))
    targets_validation = list(set(val_inter["Target_ChEMBLID"].unique().tolist()) - set(targets_reunion))
    targets_validation.sort()

    """Write training sets"""
    prot_train = inter_train[["Target_ChEMBLID"]].merge(map_file, how="inner", left_on="Target_ChEMBLID", right_on=1)
    prot_train = prot_train.drop_duplicates(subset=1).drop(labels=1, axis=1)
    prot_train.to_csv(RESULTS_FOLDER+"/train_targets.txt", sep=";", header=False, index=False, mode="w")
    inter_train.to_csv(RESULTS_FOLDER+"/train_interactions.txt", sep=";", header=True, index=False, mode="w")

    with open(RESULTS_FOLDER+"/train_ligands.txt", "w") as output_file:
        file_writer = writer(output_file, delimiter="\n")
        file_writer.writerow(inter_train["Molecule_ChEMBLID"].unique().tolist())

    """Write testing sets"""
    prot_test = inter_test[["Target_ChEMBLID"]].merge(map_file, how="inner", left_on="Target_ChEMBLID", right_on=1)
    prot_test = prot_test.drop_duplicates(subset=1).drop(labels=1, axis=1)
    prot_test.to_csv(RESULTS_FOLDER+"/test_targets.txt", sep=";", header=False, index=False, mode="w")
    inter_test.to_csv(RESULTS_FOLDER+"/test_interactions.txt", sep=";", header=True, index=False, mode="w")

    with open(RESULTS_FOLDER+"/test_ligands"+".txt", "w") as output_file:
        file_writer = writer(output_file, delimiter="\n")
        file_writer.writerow(inter_test["Molecule_ChEMBLID"].unique().tolist())

    """Write validation sets"""
    prot_validation = pd.DataFrame(targets_validation, columns=["Target_ChEMBLID"])
    prot_validation = prot_validation[["Target_ChEMBLID"]].merge(map_file, how="inner", left_on="Target_ChEMBLID", right_on=1)
    prot_validation = prot_validation.drop_duplicates(subset=1).drop(labels=1, axis=1)
    prot_validation.to_csv(RESULTS_FOLDER+"/validation_unique_targets.txt", sep=";", header=False, index=False, mode="w")
    val_inter.to_csv(RESULTS_FOLDER+"/validation_all_interactions.txt", sep=";", header=True, index=False, mode="w")
    val_inter_by_prot.to_csv(RESULTS_FOLDER+"/validation_inter_prots.txt", sep=";", header=True, index=False, mode="w")
    val_inter_by_lig.to_csv(RESULTS_FOLDER+"/validation_inter_ligs.txt", sep=";", header=True, index=False, mode="w")
    val_prot_lig.to_csv(RESULTS_FOLDER+"/validation_inter_prots_ligs.txt", sep=";", header=True, index=False, mode="w")

    with open(RESULTS_FOLDER+"/validation_unique_ligands.txt", "w") as output_file:
        file_writer = writer(output_file, delimiter="\n")
        file_writer.writerow(ligands_validation)
