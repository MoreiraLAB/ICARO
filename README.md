# ICARO
GitHub Repository from ICARO

In this Repository you can find the information required to replicate ICARO: an ensemble of Extreme Randomized Trees for Predicition of Protein-Ligand IC50.

To replicate this study you have:
A) 2 files regarding the dataset:
  1) IC50_values.csv - the main dataset where Protein-Ligand Interactions (PLI) and their binding affinities in terms of pIC50 are described, as well as some useful information   on the protein and ligand involved in each interaction.
  2) unique_uniprotid.fa - a fasta file representing all unique protein sequences included in the main dataset (IC50_values.csv).
     
B) Scripts for Feature Extraction (herein organized step-by-step):
  1) protein_ligand_sets.py - to construct the train, test and validation sets. For each subsetset it writes a .txt file with the unique protein, unique ligands and unique interactions included in that set.
  2) FeatureExtraction_ifeature.py - to export protein features from iFeature. This script receives a fasta file with all protein sequences from a dataset (in this case, unique_uniprot.fa)
  3) h5_ifeature_final.py - for iFeature features' normalization. It receives .txt of proteins in the train, test and validation sets and writes a .h5 file with iFeature descriptors normalized by train.
  4) FeatureExtraction-mordred.py - to export ligand features from MORDRED. This script access ligand's smiles through the "Canonical Smiles" column from the main dataset and exports a .h5 file ("")
  5) Folder fp_admet master - uses the fp_admet.sh to extract fp-admet descriptors.
  6) process_admet_features.py - to normalize fp-admet descriptors. It receives .txt of ligands in the train, test and validation sets and writes a .h5 file with ADMET descriptors normalized by train.
  7) FeatureExtraction-family.py - to construct a one-hot encoding feature for protein's family, based on ChEMBL's classification. This script accepts two .csv files: one with protein's ID ("single_proteins_id_noduplicate.csv") and other with the ChEMBL's family classification ("targets_family.csv"). It writes a .txt of each family name and a .h5 file saving the corresponding values.

C) Scripts for model development (the 3 best models are included here):
  1) model_random_forest.py - Construction of a Random Forest (RF). It accesses the .h5 files to retrieve interactions' features.
  2) model_xgboost.py - Construction of an Extreme Gradient Boosting algorithm. It accesses the .h5 files to retrieve interactions' features.
  3) model_xtrees.py - Construction of an Ensemble of Extreme Randomized Trees. It accesses the .h5 files to retrieve interactions' features.
