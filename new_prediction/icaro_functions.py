import os
import icaro_variables
import pandas as pd
import pickle 
import numpy as np
import h5py
from os.path import exists

####### iFeature Feature Extraction Functions
def run_ifeature(input_file, RESULTS_FOLDER, output_prefix = ""):

	"""
	Run iFeature for the different feature blocks
	"""
	output_loc = RESULTS_FOLDER
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


####### iFeature Feature Normalization Functions

def open_txt(txt_file):
    file = open(txt_file, "r")
    lines = file.readlines()
    list = []
    for line in lines:
        line = line.replace("\n", "")
        list.append(str(line))
    return list


def select_lines(file, predict_set):
    name = file.split(".")[0]
    file_pd = pd.read_csv(file, header=0, sep="\t")
    file_pd.columns = ["id"]+[name+"_"+str(i) for i in range(1,len(file_pd.columns))]
    """Update validation_table"""
    predict_table = file_pd[file_pd.iloc[:,0].isin(predict_set)==True].reset_index(drop=True)
    return predict_table

def load_model(filename):
    filehandler = open(filename, 'rb') 
    obj = pickle.load(filehandler)
    filehandler.close()
    return obj

def standardization_ifeature(predict, file, path):
	norm =[i.split("_")[0] for i in os.listdir(path)]
	if file in norm:
		print("Normalizing "+file)
		scaler = load_model(path+file+'_parameters_training_set.sav') 
		predict_transformed = scaler.transform(predict)
		return predict_transformed
	else: # feature blocks without standardization
		return predict

def process_variance(predict, name, PICKLES_FOLDER, counter=0):
    """This function checks drops columns without variance, given a pandas table.
    It returns the train, test and validation with the feature set that have variance"""
    file = pd.read_csv(PICKLES_FOLDER+"iFeature_all_variances.txt",sep=",", header = None)
    file = file[file[0].str.startswith(name)]
    remove = [int(i.split("_")[1])-1 for i in list(file[0][file[1]<= 0.01])]
    predict= np.delete(predict, remove, axis=1)
    return predict

def retrieve_features(files, pca_files, predict_ids, PICKLES_FOLDER, reduction=""):
    for file in files:
        name = file.split(".")[0]
        predict_lines = select_lines(file, predict_ids)
        """Remove variance and normalize"""
        if file == files[0]:
            predict_data = predict_lines.iloc[:,0]
        if name in pca_files:
        	predict = process_variance(predict_lines.iloc[:,1:].to_numpy(), name, PICKLES_FOLDER)
        	predict_standard = standardization_ifeature(predict, name, PICKLES_FOLDER)
        elif predict_lines.shape[1] > 1:
            predict_standard = standardization_ifeature(predict_lines.iloc[:,1:].to_numpy(), name, PICKLES_FOLDER)
        else:
            continue

        predict_lines = pd.DataFrame(predict_standard,
                           columns=[name+"_"+str(i) for i in range(1, np.shape(predict_standard)[1]+1)])
        predict_data = pd.concat((predict_data, predict_lines), axis=1)
    return predict_data

def write_h5(h5_file, dataframe):#, mode="train"):
    """Given a data frame and a mode creates a group in a h5 file, where a key is one line of the data frame"""
    #group = h5_file.create_group(mode)
    for index, row in dataframe.iterrows():
        if row.iloc[0] not in h5_file.keys():
            h5_file.create_dataset(name= row.iloc[0],
                                 data= row.iloc[1:].to_numpy(dtype="float32"))
    return dataframe.columns[1:]


####### Family Information Feature Extraction Functions

def onehot_fam_feature(targets_data, fam_list):
    targets = targets_data["Target_ChEMBLID"].unique().tolist()
    df = pd.DataFrame(0, index=targets, columns = fam_list)
    for index, row in df.iterrows():
        one_prot_data = targets_data.loc[targets_data["Target_ChEMBLID"]==index]
        one_prot_fam = one_prot_data["Class_l1"].unique().tolist()
        row[one_prot_fam] = 1
        df.update(row)
    return df

####### fp-ADMET Feature Normalization Functions

def standardization_admet(predict, file, path):
    from sklearn.preprocessing import StandardScaler
    scaler = load_model(path+file+'_parameters_training_set.sav')
    """Apply normalization"""
    prediction_transformed = scaler.transform(predict)
    return prediction_transformed

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

def process_admet(file_name, lig_predict_ids, PICKLES_FOLDER):
    file = pd.read_csv(file_name, sep=" ", header=0)
    print(file_name)
    nfeature = file_name.replace(".txt","").split("_")[1]
    if nfeature in binary_classes.keys():
        file["Predicted_"+nfeature] = file.iloc[:,0].map(binary_classes[nfeature])
    else:
        import numpy as np
        predict = file.loc[lig_predict_ids, "Predicted"].to_frame()
        
        predict_norm = standardization_admet(predict.to_numpy(), "admet"+nfeature, path=PICKLES_FOLDER)
        predict_values = pd.DataFrame(predict_norm, index=predict.index.tolist(), columns=["Predicted_"+nfeature])

        file["Predicted_"+nfeature]= predict_values
        #print(file.iloc[:10,:])
    return file["Predicted_"+nfeature]


####### Mordred Feature Extraction Functions

def retrieve_unique_smiles(input_file_name, smile_column = "Canonical_Smiles", id_column = "Molecule_ChEMBLID"):

	"""
	Open the file containing the protein-drug pairs and respective IC50 values
	"""
	opened_file = pd.read_csv(input_file_name, sep = icaro_variables.CSV_SEP, header = 0, usecols = [smile_column, id_column]).dropna()
	return opened_file

def fetch_mordred_features(input_table, smile_column = "Canonical_Smiles", id_column = "Molecule_ChEMBLID", \
							write_mode = True, verbose = True, output_folder = icaro_variables.FEATURES_MORDRED_FOLDER):

	"""
	Iterate over SMILEs list and retrieve mordred features, save to output_file
	"""
	from rdkit import Chem
	from mordred import Calculator, descriptors
	input_table = input_table.drop_duplicates(id_column, keep = "first")
	for index, row in input_table.iterrows():
		output_name = output_folder + icaro_variables.SYSTEM_SEP + row[id_column].replace("/", "-") + icaro_variables.CSV_TERMINATION
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
			processed_molecule_dataframe.to_csv(output_name, sep = icaro_variables.CSV_SEP, index = False)
		if verbose == True:
			print("Currently extracting mordred features for entry", index + 1, icaro_variables.SYSTEM_SEP, input_table.shape[0])



####### Mordred Feature Normalization Functions

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
			with h5py.File(self.source, "r") as input_file:
				if self.usable_entries == "all":
					self.number_of_entries = len(list(input_file))
				elif self.usable_entries != "all":
					self.number_of_entries = len(self.usable_entries)

	def retrieve_sums(self):

		if self.format == "h5py":
			sums = generate_placeholder_dictionary(self.header, mode = "array")
			with h5py.File(self.source, "r") as input_file:
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
			with h5py.File(self.source, "r") as input_file:
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
		print(self)
		with h5py.File(self.source, "r") as input_file, h5py.File(output_file_name, "w") as output_file:
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
			print(len(normalized_values))
