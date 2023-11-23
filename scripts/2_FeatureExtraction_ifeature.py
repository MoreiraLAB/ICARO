#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DESCRIPTION: Call the iFeature commands.

NOTE: Scikit-learn and matplotlib packages are needed!
make sure that you have iFeature folder in "./feature/" folder, you can download this folder with:
> git clone https://github.com/Superzchen/iFeature

INPUT FILES: (path: "./data/") "unique_uniprotid.fa"
OUTPUT FILES: (path: "./features/results_ifeature/") 
for each ifeature feature, a .tsv file is saved on "./features/results_ifeature/" folder
"""

__author__ = "A.J. Preto"
__email__ = "ateresa.work@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "ICARO: IC50 binding Affinity Regression Optimized"



import os
import icaro_variables

def run_ifeature(input_file, output_prefix = ""):

	"""
	Run iFeature for the different feature blocks
	"""
	output_loc = RESULTS_FOLDER + SYSTEM_SEP
	#excluded: EAAC, BINARY, EGAAC, AAINDEX, ZSCALE, BLOSUM62

	for feature in IFEATURE_FEATURES:
		if feature == 'SOCNumber' or feature == 'QSOrder' or feature == 'APAAC' or feature == 'PAAC':
			run_command = "python codes/"+feature+".py " + input_file + \
							" 8 " + output_loc + \
							output_prefix + feature + IFEATURE_TERMINATION
		else:
			run_command = "python iFeature.py --file " + input_file + \
							" --type " + feature + " --out " + output_loc + \
							output_prefix + feature + IFEATURE_TERMINATION
		os.system(run_command)
		if IFEATURES_PCA_DICT[feature] != None:
			cluster_command = "python scripts/pcaAnalysis.py --file " + output_loc + \
							output_prefix + feature + IFEATURE_TERMINATION + \
							" --ncomponents " + str(IFEATURES_PCA_DICT[feature]) + " --out " + output_loc + \
							output_prefix + feature + "_reduced" + IFEATURE_TERMINATION
			os.system(cluster_command)

	#excluded: type5, type13
	for pseudo_feature in PSEUDO_IFEATURES:
		run_second_command = "python iFeaturePseKRAAC.py --file " + input_file + \
						" --type " + pseudo_feature + " --gap_lambda 5 --raactype 5 --out " + output_loc + \
						output_prefix + pseudo_feature + IFEATURE_TERMINATION
		os.system(run_second_command)
		if IFEATURES_PCA_DICT[pseudo_feature] != None:
			second_cluster_command = "python scripts/pcaAnalysis.py --file " + output_loc + \
							output_prefix + pseudo_feature + IFEATURE_TERMINATION + \
							" --ncomponents " + str(IFEATURES_PCA_DICT[pseudo_feature]) + " --out " + output_loc + \
							output_prefix + pseudo_feature + "_reduced" + IFEATURE_TERMINATION
			os.system(second_cluster_command)



SYSTEM_SEP = "/"
FILE_LOCATION = icaro_variables.DATA_FOLDER + SYSTEM_SEP
DEFAULT_LOCATION = icaro_variables.FEATURES_FOLDER + SYSTEM_SEP
IFEATURE_FOLDER = "iFeature" # iFeature folder must be installed in features/ folder!
RESULTS_FOLDER = icaro_variables.FEATURES_FOLDER + SYSTEM_SEP + "results_iFeature"

TXT_TERMINATION = ".txt"
IFEATURE_FEATURES = ['AAC', 'CKSAAP', 'DPC', 'DDE', 'TPC', 'GAAC', 'CKSAAGP', 'GDPC', 'GTPC',
				 'CTDC', 'CTDT', 'CTDD', 'CTriad', 'KSCTriad',
				 'SOCNumber', 'QSOrder', 'PAAC', 'APAAC'] #Moran, NMBroto and Geary features were not chosen
				 # because they depend on the AAIndex.csv file, based on the AAindex descriptor, which is calculated
				 # for sequences of the same length. SSEC was not selected too, because some protein's secondary structures
				 # were predicted wrong. EAAC has to be applied to sequences with equal length

PSEUDO_IFEATURES = ['type1', 'type2', 'type3A', 'type3B', 'type4', 'type6A', 'type6B', 'type6C',
								 'type7', 'type8', 'type9', 'type10', 'type11', 'type12', 'type14',
								 'type15', 'type16']
IFEATURES_PCA_DICT = {'AAC': None, 'EAAC': None, 'CKSAAP': 50, 'DPC': 40, 'DDE': 60, 'TPC': 100,
				'GAAC': None, 'CKSAAGP' : 40, 'GDPC': 10, 'GTPC': 40,
				'NMBroto': 50, 'Moran': 50, 'Geary': 40,
				'CTDC': 20, 'CTDT': 20, 'CTDD': 40, 'CTriad': 40, 'KSCTriad': 40, 'SOCNumber': None,
				'QSOrder': None, 'PAAC': None, 'APAAC': None, 'SSEC': None, 'type1': 10, 'type2': 10, 'type3A': 10,
				'type3B': 20, 'type4': 10, 'type5':10, 'type6A': 10, 'type6B': 10, 'type6C': 10, 'type7': 10,
				'type8': 10, 'type9': 10, 'type10': 10, 'type11': 10, 'type12': 10, 'type13': 10, 'type14': 10,
				'type15': 10, 'type16': 10}
IFEATURE_TERMINATION = ".tsv"

if not os.path.exists(RESULTS_FOLDER + SYSTEM_SEP):
    os.makedirs(RESULTS_FOLDER + SYSTEM_SEP)

print("#### ICARO- Feature Extraction iFeature")
ifeature_path = DEFAULT_LOCATION + IFEATURE_FOLDER
fasta_file = FILE_LOCATION + "unique_uniprotid.fa" # file with proteins in fasta format

os.chdir(ifeature_path)
run_ifeature(fasta_file)
