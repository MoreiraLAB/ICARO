#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DESCRIPTION: Normalize a dataset via a batch processing for large datasets.
Stores the standard deviation and average for each column. 

INPUT FILES: (path: "./features/") "icaro_mordred_features.h5"
OUTPUT FILES: 
- (path: "./features/") "mordred_norm.csv" with standard deviation and average for each column.
- (path: "./features/") "icaro_mordred_normalized_features.h5"
"""

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "Tools"

import os
import pandas as pd
import numpy as np
import sys
import h5py as h5

def retrieve_ids_list(input_file, file_type = "txt"):

	"""
	I mean, the name of the function is pretty self-explanatory
	"""
	if file_type == "txt":
		opened_file = pd.read_csv(input_file)
		return list(opened_file.iloc[:,0])

def generate_placeholder_dictionary(input_list, mode = "array"):

	"""
	Output a zero-ed dictionary that can later be updated with values
	"""
	if mode == "dictionary":
		return {x:0.0 for x in input_list}
	if mode == "array":
		return np.zeros(len(input_list))

def convert_to_dictionary(input_values_list, input_header):
	
	"""
	Reorganize a list into a named dictionary
	"""
	output_dictionary = {}
	for current_name, current_value in zip(input_header, input_values_list):
		output_dictionary[current_name] = current_value
	return output_dictionary

class normalizer:

	"""
	Deploy the normalizer depending on:
	- The input format:
		- "h5py": an input h5py with the key as entries, can be fed an associated header file for proper descriptor identification.
		The "source" argument must be fed the input file.
		- "csv": will search for several csv files, open them and register each value. 
		The "source" argument must be fed a folder with the input csv files.
	- Logging variables:
		- True: Logs a file, requires a filename to log the results into. Default logging file "norm.csv"
		- False: Does not log a file, not recommended
	- Logging file type:
		- "csv": Logs the file onto a table
		- "dictionary": Saves a pickle dictionary that can later be loaded for direct deployment
	- Normalizing type:
		- "mean_and_scale": standardized_value = (value - average) / standard_deviation
	- Header names:
		- "input_header": True will search for an input file, False will assume the order of the columns as their names
		- "header_file": pass a csv file, its header will be considered the header for the whole dataset.
		Pass a list of columns to ignore if applicable.
	- Usable entries: If not all entries are supposed to be used, pass a list on this argument
	"""

	def __init__(self, input_format = "h5py", verbose = True, \
						log_file_name = "norm.csv", log_file_type = "csv", \
						source = "", normalization_type = "mean_and_scale", \
						input_header = True, header_file = "", \
						droppable_columns = [], usable_entries = "all"):

		self.format = input_format	
		self.log_file_name = log_file_name
		self.log_file_type = log_file_type
		self.source = source
		self.normalization_type = normalization_type
		self.input_header = input_header
		self.header_file = header_file
		self.droppable_columns = droppable_columns
		self.usable_entries = usable_entries
		self.verbose = verbose
		
	def fetch_number_of_columns(self):

		if self.format == "h5py":
			with h5.File(self.source, "r") as input_file:
				if self.usable_entries == "all":
					target_list = list(input_file)
				if self.usable_entries != "all":
					target_list = self.usable_entries 
				for entry in target_list:
					entry_values = list(input_file[entry])
					self.number_of_columns = len(entry_values)
					break

	def attach_header(self):

		if self.input_header == True:
			self.header = list(pd.read_csv(self.header_file, header = 0, nrows = 2).drop(self.droppable_columns, axis = 1))

		elif self.input_header == False:
			self.fetch_number_of_columns()
			self.header = list(range(0, self.number_of_columns))

	def retrieve_length(self):

		if self.format == "h5py":
			with h5.File(self.source, "r") as input_file:
				if self.usable_entries == "all":
					self.number_of_entries = len(list(input_file))
				elif self.usable_entries != "all":
					self.number_of_entries = len(self.usable_entries)

	def retrieve_sums(self):

		if self.format == "h5py":
			sums = generate_placeholder_dictionary(self.header, mode = "array")
			with h5.File(self.source, "r") as input_file:
				if self.usable_entries == "all":
					target_list = list(input_file)
				if self.usable_entries != "all":
					target_list = self.usable_entries 
				for index, entry in enumerate(target_list):
					if ((index + 1) % 100) == 0 and self.verbose == True:
						print("Sums for averages calculation progress:",index + 1)
					entry_values = input_file[entry][:]
					sums = np.sum([sums, entry_values], axis = 0)
			self.sums = convert_to_dictionary(sums, self.header)

	def retrieve_averages(self):
		self.retrieve_sums()
		self.retrieve_length()
		self.averages = generate_placeholder_dictionary(self.header, mode = "dictionary")
		if self.format == "h5py":
			for current_column in self.header:
				self.averages[current_column] = float(self.sums[current_column]) / float(self.number_of_entries)

	def retrieve_standard_deviations(self):
		
		self.retrieve_averages()
		standard_deviation_terms = generate_placeholder_dictionary(self.header, mode = "array")
		if self.format == "h5py":
			with h5.File(self.source, "r") as input_file:
				if self.usable_entries == "all":
					target_list = list(input_file)
				if self.usable_entries != "all":
					target_list = self.usable_entries 
				for index, entry in enumerate(target_list):
					if ((index + 1) % 100) == 0 and self.verbose == True:
						print("Standard deviation calculation progress:",index + 1)
					entry_values = input_file[entry][:]
					standard_sum_terms = np.square(np.subtract(entry_values, np.array(list(self.averages.values()))))
					standard_deviation_terms = np.add(standard_deviation_terms, standard_sum_terms)

		standard_deviation_terms = convert_to_dictionary(standard_deviation_terms, self.header)
		self.standard_deviations = generate_placeholder_dictionary(self.header, mode = "dictionary")
		for current_column in self.header:
			self.standard_deviations[current_column] = (float(standard_deviation_terms[current_column]) / float(self.number_of_entries))**float(1/2)
	
	def report(self):

		if self.normalization_type == "mean_and_scale":
			averages_table = pd.DataFrame.from_dict(self.averages.items())
			averages_table.columns = ["column_name","average_value"]

			standard_deviations_table = pd.DataFrame.from_dict(self.standard_deviations.items())
			standard_deviations_table.columns = ["column_name", "standard_deviation"]

			output_table = averages_table.merge(standard_deviations_table, on = "column_name")
			output_table.to_csv(self.log_file_name, index = False)

	def filter_normalization(self, exclude_useless = True, \
				save_new = False, new_norm_csv = ""):

		"""
		On the normalization table, apply filters, such as exclude some values and drop features with standard deviation 0.0
		Generates the self.average_values and self.standard_deviation_values attributes, which can later be applied to normalize new samples
		"""
		opened_table = pd.read_csv(self.log_file_name, sep = ",")
		self.useless_rows = []
		if exclude_useless == True:
			self.useless_rows = opened_table[opened_table["standard_deviation"].isin([0.0, np.inf])].index.tolist()
			opened_table = opened_table.transpose()
			opened_table = opened_table.drop(opened_table.columns[self.useless_rows], axis = 1)
		
		opened_table.columns = opened_table.iloc[0,:]
		opened_table = opened_table.iloc[1:,:]
		self.average_values = opened_table.loc["average_value"].values
		self.standard_deviation_values = opened_table.loc["standard_deviation"].values
		if save_new == True:
			opened_table.to_csv(new_norm_csv, index = False)

	def apply_normalization(self, output_file_name = "", output_format = "h5py"):

		if self.format == "h5py" and output_format == "h5py":
			with h5.File(self.source, "r") as input_file, h5.File(output_file_name, "w") as output_file:
				if self.usable_entries == "all":
					self.usable_entries = list(input_file)
				self.number_of_entries = len(self.usable_entries)
				for index, current_entry in enumerate(self.usable_entries):
					print(current_entry)
					retrieved_values = np.delete(input_file[current_entry][:], self.useless_rows)
					if self.normalization_type == "mean_and_scale":
						normalized_values = np.divide(np.subtract(retrieved_values, self.average_values), self.standard_deviation_values)
						current_dataset = output_file.create_dataset(current_entry, dtype = "float", data = normalized_values.astype(float)) 
					if (self.verbose == True) and ((index + 1) % 100 == 0):
						print("Normalized values at", index + 1, "/", self.number_of_entries)

from icaro_variables import TRAIN_SPLITS_FILE, MORDRED_FEATURES_FILE_H5, MORDRED_NORMALIZED_FEATURES_FILE_H5, \
							MORDRED_NORMALIZED_FEATURES_FILE_H5, MORDRED_FEATURES_FILE_H5, SUPPORT_FOLDER

"""
Generate the file with the normalization metrics
"""

input_header_file = SUPPORT_FOLDER+"/CHEMBL2.csv" # Example file containing the header

"""
Run first block to calculate normalization vectors, second block to actually normalize
"""
train_ids_list = retrieve_ids_list(TRAIN_SPLITS_FILE)
print("#### ICARO- MORDRED Normalization")

normalizer_object = normalizer(source = MORDRED_FEATURES_FILE_H5, \
						header_file = input_header_file, \
						droppable_columns = ["ID", "SMILE"], usable_entries = train_ids_list)
normalizer_object.attach_header()
normalizer_object.retrieve_standard_deviations()
normalizer_object.report()

to_normalize_object.filter_normalization(save_new = True, new_norm_csv = "mordred_norm.csv")
to_normalize_object.apply_normalization(output_file_name = MORDRED_NORMALIZED_FEATURES_FILE_H5)
