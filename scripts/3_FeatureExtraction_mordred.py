#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DESCRIPTION: Retrieve mordred features for all the single files.

NOTE: Initiate icaro_env environment. This environment needs to be python version 3.7.x or below
- pip install rdkit-pypi
- pip install mordred 
- pip install pandas

INPUT FILES: (path: "./data/support/") "Davis_comp_smiles.csv"
OUTPUT FILES: 
- path: "./features/davis_mordred/" with all mordred descriptors
- (path: "./features/") "icaro_davis_mordred_features.h5" an h5 file with mordred features
"""

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "MENSA"

import os
from os.path import exists
import sys
import pandas as pd
import numpy as np							
import h5py as h5
from icaro_variables import ICARO_BASE_FILE, SYSTEM_SEP, CSV_SEP, \
							MORDRED_FEATURES_FILE, FEATURES_MORDRED_FOLDER, \
							CSV_TERMINATION, FEATURES_DAVIS_MORDRED_FOLDER, SUPPORT_FOLDER, \
							MORDRED_FEATURES_FILE_H5, MORDRED_DAVIS_FEATURES_FILE_H5

def retrieve_unique_smiles(input_file_name, smile_column = "Canonical_Smiles", id_column = "Molecule_ChEMBLID"):

	"""
	Open the file containing the protein-drug pairs and respective IC50 values
	"""
	opened_file = pd.read_csv(input_file_name, sep = CSV_SEP, header = 0, usecols = [smile_column, id_column]).dropna()
	return opened_file

def fetch_mordred_features(input_table, smile_column = "Canonical_Smiles", id_column = "Molecule_ChEMBLID", \
							write_mode = True, verbose = True, output_folder = FEATURES_MORDRED_FOLDER):

	"""
	Iterate over SMILEs list and retrieve mordred features, save to output_file
	"""
	from rdkit import Chem
	from mordred import Calculator, descriptors
	input_table = input_table.drop_duplicates(id_column, keep = "first")
	for index, row in input_table.iterrows():
		output_name = output_folder + SYSTEM_SEP + row[id_column].replace("/", "-") + CSV_TERMINATION
		if exists(output_name):
			continue
		current_molecule = Chem.MolFromSmiles(row[smile_column])
		calculator = Calculator(descriptors, ignore_3D = False)
		features_dictionary = calculator(current_molecule)
		current_dataframe = pd.DataFrame.from_dict(features_dictionary, orient = "index")
		clean_dataframe = pd.to_numeric(current_dataframe[0], errors ='coerce').fillna(0).astype(float)
		processed_molecule_dataframe = current_dataframe.transpose()

		processed_molecule_dataframe["ID"] = row[id_column]
		processed_molecule_dataframe["SMILE"] = row[smile_column]
		if write_mode == True:
			processed_molecule_dataframe.to_csv(output_name, sep = CSV_SEP, index = False)
		if verbose == True:
			print("Currently extracting mordred features for entry", index + 1, SYSTEM_SEP, input_table.shape[0])

def merge_files(input_folder = FEATURES_MORDRED_FOLDER, \
				output_file = MORDRED_FEATURES_FILE_H5, verbose = False, \
				droppable_columns = ["ID", "SMILE"]):

	"""
	Fetch the pre-calculated files from the results folder and merge them into a single file
	"""
	with h5.File(output_file, "w") as h5_file:
		for index, files in enumerate(os.listdir(input_folder)):
			if files.endswith(CSV_TERMINATION):
				if (index + 1) % 100 == 0:
					print("Currently merging", index + 1, "/", len(os.listdir(input_folder)))
				opened_file = pd.read_csv(input_folder + SYSTEM_SEP + files)
				smile_value = opened_file["ID"][0]
				processed_table = np.nan_to_num(pd.to_numeric(list(opened_file.drop(droppable_columns, axis = 1).values[0]), errors = "coerce"), nan = 0.0)
				current_dataset = h5_file.create_dataset(smile_value.replace("/",""), dtype = "float", data = processed_table) 

if not os.path.exists(FEATURES_DAVIS_MORDRED_FOLDER):
    os.makedirs(FEATURES_DAVIS_MORDRED_FOLDER)
print("#### ICARO- Feature Extraction MORDRED")
#fetch_mordred_features(retrieve_unique_smiles(ICARO_BASE_FILE))
fetch_mordred_features(retrieve_unique_smiles(SUPPORT_FOLDER+"/Davis_comp_smiles.csv", smile_column = "SMILES", id_column = "Comp_Name"), \
	output_folder = FEATURES_DAVIS_MORDRED_FOLDER, smile_column = "SMILES", id_column = "Comp_Name")

merge_files(input_folder = FEATURES_DAVIS_MORDRED_FOLDER, output_file = MORDRED_DAVIS_FEATURES_FILE_H5, verbose = True)
