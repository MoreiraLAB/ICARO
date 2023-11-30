import os
from icaro_variables import DATA_FOLDER, FEATURES_FOLDER, NORMALIZE_FOLDER, FEATURES_MORDRED_FOLDER, CSV_TERMINATION, MODEL_FOLDER, MORDRED_FOLDER
import icaro_functions
import pandas as pd
import h5py
import shutil
import subprocess
import numpy as np

####### iFeature Feature Extraction
#INPUT FILES: (path: "./new_prediction/") "protein_sequence_predict.fa"
#OUTPUT FILES: (path: "./features/results_ifeature/")
#for each ifeature feature, a .tsv file is saved on "./results/results_ifeature/" folder
#IMPORTANT: make sure that you have iFeature folder in "./feature/" folder, you can download this folder with:
#git clone https://github.com/Superzchen/iFeature

SYSTEM_SEP = "/"
IFEATURE_FOLDER = "iFeature" # iFeature folder must be installed in features/ folder!
RESULTS_FOLDER = FEATURES_FOLDER + SYSTEM_SEP + "results_iFeature/"

if not os.path.exists(RESULTS_FOLDER + SYSTEM_SEP):
    os.makedirs(RESULTS_FOLDER + SYSTEM_SEP)


print("#### ICARO- 1: Feature Extraction iFeature")
ifeature_path = FEATURES_FOLDER + SYSTEM_SEP + IFEATURE_FOLDER
PICKLES_FOLDER= NORMALIZE_FOLDER + SYSTEM_SEP + "iFeature_norm/"
fasta_file = DATA_FOLDER + "protein_sequence_predict.fa"
feature_list = pd.read_csv(PICKLES_FOLDER +"/iFeature_header.txt", sep=",")
os.chdir(ifeature_path)
icaro_functions.run_ifeature(fasta_file, RESULTS_FOLDER)

####### iFeature Feature Normalization
#INPUT FILES: (path: "./features/results_ifeature/")
#OUTPUT FILES: (path: "./features/") "icaro_prediction_normalized_ifeature.h5" with ifeature features

print("#### ICARO- 1: Feature Normalization iFeature")
fasta_file = icaro_functions.open_txt(fasta_file)
prt_predict_ids = [i.replace(">","") for i in fasta_file if ">" in i]

os.chdir(RESULTS_FOLDER)
dir_files = os.listdir(RESULTS_FOLDER)
dir_files.sort()

dir_files = [i for i in dir_files]

for_pca = ["CTriad","KSCTriad"]
predict_data = icaro_functions.retrieve_features(dir_files, for_pca, prt_predict_ids, PICKLES_FOLDER)
predict_data = predict_data[feature_list.columns]
print("#### iFeature prediction shape:",predict_data.shape)
h5_ifeature = h5py.File(FEATURES_FOLDER+"/icaro_prediction_normalized_ifeature.h5", "w")
icaro_functions.write_h5(h5_ifeature, predict_data)

####### Family Information Feature Extraction
#INPUT FILES: (path: "./new_prediction/") "protein_sequence_predict.fa"

print("#### ICARO- 2: Protein Family Information Features Extraction")
FAM_PRT_FOLDER= FEATURES_FOLDER + SYSTEM_SEP + "family_prt/"
mapping = pd.read_csv(FAM_PRT_FOLDER+"unique_prt_ids.csv", sep=",", header=0, usecols=["Uniprot ID","ChEMBL ID"])
mapping.index = [mapping["Uniprot ID"].values]
fam_sheet = pd.read_csv(FAM_PRT_FOLDER+"targets_family.csv", sep=",", header=0, usecols=["Target_ChEMBLID", "Class_l1"])
fam_data = fam_sheet.drop_duplicates()
many_class_targets = fam_data.loc[fam_data.duplicated(subset="Target_ChEMBLID",keep=False)==True] #40 targets have more than 1 family
families = fam_data["Class_l1"].unique().tolist() #There are 15 unique families
families.sort()

#Write dataframe with one-hot-encoding representation of proteins' families
final_data = icaro_functions.onehot_fam_feature(fam_data, families)
final_data = final_data.loc[:, final_data.var(axis=0, numeric_only=True) != 0]
final_data.insert(0, "index", final_data.index)
cols = final_data.columns[1:].tolist()
final_data = final_data.merge(mapping, how="inner", right_on="ChEMBL ID", left_on="index")
neworder = ["Uniprot ID"]+cols
final_data=final_data.reindex(columns=neworder)
final_data = final_data.loc[final_data['Uniprot ID'].isin(prt_predict_ids)]
print("#### Family Information prediction shape:",final_data.shape)
output_file = h5py.File(FEATURES_FOLDER+"/icaro_prediction_normalized_prt_family.h5", "w")
icaro_functions.write_h5(output_file, final_data)


####### fp-ADMET Feature Extraction
#INPUT FILES: (path: "./new_prediction/") "ligand_smiles_predict.smi"
#OUTPUT FILES: (path: "./features/fpadmet/RESULTS/")
#NOTE: fp-admet must be installed in "./features/" folder with:
#> git clone https://github.com/jcheminform/fpadmet.git

print("#### ICARO- 3: fp-ADMET Features Extraction")
print("#### Coping fp_admet.sh file to ./features/fpadmet/ folder")
bash_file = DATA_FOLDER +"fp_admet.sh"
ADMET_FOLDER = FEATURES_FOLDER + "/fpadmet/"
PICKLES_FOLDER= NORMALIZE_FOLDER + SYSTEM_SEP + "fpadmet_norm/"
shutil.copy2(bash_file, ADMET_FOLDER)

os.chdir(ADMET_FOLDER)
print("#### Running fp-ADMET features")
subprocess.run("bash fp_admet.sh", shell=True)
os.chdir(DATA_FOLDER)


####### fp-ADMET Feature Normalization
#INPUT FILES: (path: "./features/fpadmet/RESULTS/")
#OUTPUT FILES: "icaro_prediction_normalized_admet.h5"

print("#### ICARO- 3: fp-Admet Normalization")
smile_file = DATA_FOLDER + "ligand_smiles_predict.smi"
smile_file =  icaro_functions.open_txt(smile_file)
lig_predict_ids = [i.split("\t")[1] for i in smile_file]
lig_predict_smiles = [i.split("\t")[0] for i in smile_file]

#Change to ADMET directory and open files. All features are binary or a class except
#from the features 50 to 56, inclusive, that are a numeric value
ADMET_RESULTS_FOLDER = ADMET_FOLDER +"RESULTS/"
os.chdir(ADMET_RESULTS_FOLDER)
files = os.listdir()
files.sort()
files = [i for i in files if i.startswith("predicted")]

final_table = icaro_functions.process_admet(files[0], lig_predict_ids, PICKLES_FOLDER).to_frame()

for f in files[1:]:
    output_feat = icaro_functions.process_admet(f, lig_predict_ids, PICKLES_FOLDER)

    #if output_feat.isnull().values.any() == True:
    #    print(f, output_feat.isnull().values.any())
    if output_feat is None:
        continue
    if output_feat.name in final_table.columns:
        final_table[output_feat.name].update(output_feat)
    else:
        final_table = pd.concat([final_table.reset_index(drop=True), output_feat.reset_index(drop=True)], axis=1, join='outer')
    #print(final_table.iloc[:10,:])
    #if final_table.isnull().values.any() == True:
    #    print(f, final_table)
final_table.insert(0, 'id', lig_predict_ids)
print("#### fp-ADMET prediction shape:",final_table.shape)
h5_admet = h5py.File(FEATURES_FOLDER+"/icaro_prediction_normalized_admet.h5", "w")
icaro_functions.write_h5(h5_admet, final_table)


####### Mordred Feature Extraction
#INPUT FILES: (path: "./new_prediction/") "ligand_smiles_predict.smi"
#OUTPUT FILES: "icaro_prediction_mordred.h5"

if not os.path.exists(FEATURES_MORDRED_FOLDER+"/"):
    os.makedirs(FEATURES_MORDRED_FOLDER+"/")

print("#### ICARO- 4: Feature Extraction MORDRED")

#fetch_mordred_features(retrieve_unique_smiles(ICARO_BASE_FILE))

icaro_functions.fetch_mordred_features(pd.DataFrame({'Canonical_Smiles': lig_predict_smiles, 'Molecule_ChEMBLID': lig_predict_ids}), \
	output_folder = FEATURES_MORDRED_FOLDER)

droppable_columns = ["ID", "SMILE"]
with h5py.File(FEATURES_FOLDER+"/icaro_prediction_mordred.h5", "w") as h5_file:
	for index, files in enumerate(os.listdir(FEATURES_MORDRED_FOLDER)):
		if files.endswith(CSV_TERMINATION):
			if (index + 1) % 100 == 0:
				print("Currently merging", index + 1, "/", len(os.listdir(input_folder)))
			opened_file = pd.read_csv(FEATURES_MORDRED_FOLDER + SYSTEM_SEP + files)
			smile_value = opened_file["ID"][0]
			processed_table = np.nan_to_num(pd.to_numeric(list(opened_file.drop(droppable_columns, axis = 1).values[0]), errors = "coerce"), nan = 0.0)
			current_dataset = h5_file.create_dataset(smile_value.replace("/",""), dtype = "float", data = processed_table) 

print("#### Mordred prediction shape:",processed_table.shape)

####### Mordred Feature Normalization
#INPUT FILES: "icaro_prediction_mordred.h5"
#OUTPUT FILES: "icaro_prediction_normalized_mordred.h5"


#Generate the file with the normalization metrics


input_header_file = MORDRED_FOLDER+"/CHEMBL2.csv" # Example file containing the header


#Run first block to calculate normalization vectors, second block to actually normalize

print("#### ICARO- 4: MORDRED Normalization")

to_normalize_object = icaro_functions.normalizer(source = "icaro_prediction_mordred.h5", \
						header_file = input_header_file, \
						droppable_columns = ["ID", "SMILE"], usable_entries = lig_predict_ids)
to_normalize_object.attach_header()
to_normalize_object.filter_normalization(save_new = False) #,new_norm_csv = MORDRED_FOLDER+"/new_norm.csv")
to_normalize_object.apply_normalization(output_file_name = FEATURES_FOLDER+"/icaro_prediction_normalized_mordred.h5")

####### New prediction
#INPUT FILES: "icaro_prediction_mordred.h5"
#OUTPUT FILES: "icaro_prediction_normalized_mordred.h5"

print("#### ICARO- 5: ICARO prediction")

def retrieve_features(pd_class):
    output = []
    for index, row in pd_class.iterrows():
        features_row = np.array([])
        lig = row["Molecule_ChEMBLID"]
        uniprot = row["Target_ChEMBLID"] #0 is the name of the column
        with h5py.File(FEATURES_FOLDER+"/icaro_prediction_normalized_mordred.h5") as mordred:
            features_row = np.concatenate((features_row, mordred[lig][:]))
        with h5py.File(FEATURES_FOLDER+"/icaro_prediction_normalized_admet.h5") as admet:
            features_row = np.concatenate((features_row, admet[lig][:]))
        with h5py.File(FEATURES_FOLDER+"/icaro_prediction_normalized_ifeature.h5") as ifeature:
            features_row = np.concatenate((features_row, ifeature[uniprot][:]))
        with h5py.File(FEATURES_FOLDER+"/icaro_prediction_normalized_prt_family.h5") as family:
            features_row = np.concatenate((features_row, family[uniprot][:]))
        output.append(features_row)
    return output


predict_features = retrieve_features(pd.DataFrame({'Target_ChEMBLID': prt_predict_ids, 'Molecule_ChEMBLID': lig_predict_ids}))
print("Total number of features:",len(predict_features[0]))

model = icaro_functions.load_model(MODEL_FOLDER+'/XTrees.sav')
ic50 = model.predict(predict_features)
predictions = pd.DataFrame({'Target_ChEMBLID': prt_predict_ids, 'Molecule_ChEMBLID': lig_predict_ids, 'pIC50_prediction': ic50})
predictions.to_csv(DATA_FOLDER+'/new_predictions.csv', index=False)

print("#### ICARO- 6: ICARO output file with prediction saved in new_predictions.csv")
print(predictions)
