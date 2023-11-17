## ICARO (IC50 binding Affinity Regression Optimized)

## Abstract: 
The search for new pharmaceutical products has become increasingly difficult and expensive, leading the pharmaceutical industry to explore alternative methods of drug discovery. One approach involves using computational methods to reduce the number of compounds during the early stages of drug development. Protein-ligand interactions play a crucial role in understanding the mechanisms of drug action. Therefore, there is great interest in developing computational methods for determining drug-target binding affinity in the form of inhibitory concentrations (IC50), inhibition constants (Ki), or dissociation constants (Kd) for hit recognition. Although several predictors have been developed, they mainly focus on Ki and Kd, and are restricted to interactions involving kinases. Therefore, a general algorithm that can predict the IC50 values for different targets is required. Here, we propose ICAROERT (IC50 binding Affinity Regression Optimized) for drug-target binding affinity prediction, which is expressed as pIC50. ICAROERT achieved a coefficient of determination (R2) value of 0.726, root mean squared error (RMSE) of 0.701, Pearson correlation coefficient (PCC) of 0.855, and concordance index (CI) of 0.837. ICAROERT is an efficient and cost-effective method for drug discovery.

![Graphical Abstract](Graphical_Abstract.png)⁩

### Requirements:
iFeature - GitHub available at https://github.com/Superzchen/iFeature. You should download the iFeature folder and paste it into this git home directory.

In this Repository one can find the information required to replicate ICARO: an ensemble of Extreme Randomized Trees for Prediction of Protein-Ligand IC50.

### To replicate this study, you need:
A) Two files regarding the dataset
 1) IC50_values.csv - The main dataset where Protein-Ligand Interactions (PLI) and their binding affinities in terms of pIC50 are described, as well as some useful information on the protein and ligand involved in each interaction.
 2) unique_uniprotid.fa - A FASTA file representing all the unique protein sequences included in the main dataset (IC50_values.csv).

B) Scripts for Feature Extraction ( organized step-by-step):
 1) protein_ligand_sets.py - To construct the training, test, and validation sets. For each subset, we write as. txt file with unique proteins, ligands, and interactions included in that set.
 2) FeatureExtraction_ifeature.py - to export protein features from iFeature. This script receives a fasta file with all protein sequences from a dataset (in this case, unique_uniprot.fa)
 3) h5_ifeature_final.py - For the normalization of the iFeature features it receives. txt of proteins in the training, test, and validation sets and writes a. h5 file with iFeature descriptors normalized by training.
 4) FeatureExtraction-mordred.py - To export ligand features from MORDRED. This script access ligand's smiles through the "Canonical Smiles" column from the main dataset and exports a .h5 file ("")
 5) fp_admet.sh - to extract fp-admet descriptors.
 6) process_admet_features.py - To normalize fp-admet descriptors. It receives .txt of ligands in the training, testing, and validation sets, and write a. h5 files with ADMET descriptors normalized by training.
 7) FeatureExtraction-family.py - To construct a one-hot encoding feature for the protein family, based on ChEMBL's classification. The script accepts two words. csv files: one with the protein ID ("single_proteins_id_noduplicate.csv") and other with the ChEMBL's family classification ("targets_family.csv"). It writes a .txt for each family name and a. h5 file to save the corresponding values.

C) Scripts for model development (the 3 best models are included here):
 1) model_random_forest.py - To construct a Random Forest (RF). It accesses the .h5 files to retrieve interaction features.
 2) model_xgboost.py - To construct an Extreme Gradient Boosting algorithm. It accesses the .h5 files to retrieve interaction features.
 3) model_xtrees.py - To construct an Ensemble of Extreme Randomized Trees. It accesses the .h5 files to retrieve interaction features.
