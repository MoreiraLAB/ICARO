## ICARO<sub>ERT</sub> (IC<sub>50</sub> binding Affinity Regression Optimized)

## Abstract: 
The search for new pharmaceutical products has become increasingly difficult and expensive, leading the pharmaceutical industry to explore alternative methods of drug discovery. One approach involves using computational methods to reduce the number of compounds during the early stages of drug development. Protein-ligand interactions play a crucial role in understanding the mechanisms of drug action. Therefore, there is great interest in developing computational methods for determining drug-target binding affinity in the form of inhibitory concentrations (IC<sub>50</sub>), inhibition constants (Ki), or dissociation constants (K<sub>d</sub>) for hit recognition. Although several predictors have been developed, they mainly focus on K<sub>i</sub>and K<sub>d</sub>, and are restricted to interactions involving kinases. Therefore, a general algorithm that can predict the IC<sub>50</sub> values for different targets is required. Here, we propose ICARO<sub>ERT</sub> (IC<sub>50</sub> binding Affinity Regression Optimized) for drug-target binding affinity prediction, which is expressed as pIC<sub>50</sub>. ICARO<sub>ERT</sub> achieved a coefficient of determination (R<sup>2</sup>) value of 0.726, root mean squared error (RMSE) of 0.701, Pearson correlation coefficient (PCC) of 0.855, and concordance index (CI) of 0.837. ICARO<sub>ERT</sub> is an efficient and cost-effective method for drug discovery.
⁩

### Prerequisites:
ICARO<sub>ERT</sub> was developed and tested as follows:
> Python 3.6.9 (default, Mar 10 2023, 16:46:00)

We recommend creating an isolated Conda environment to run our pipeline, which can be performed using the following code:
```bash
conda create --name icaro python=3.6.9
conda activate icaro
```
Note: The environment name, defined after the "--name" argument in the first step, can be whatever the user desires.

### Requirements:
These requirements can be met using a pip.

* iFeature - GitHub available at https://github.com/Superzchen/iFeature. The iFeature folder should be downloaded into. /feature/ folder.

* fp-admet - Github available at https://github.com/jcheminform/fpadmet. YThe fp-admet folder should be downloaded into. /feature/ folder.
  
* MORDRED - version 1.2.0 .
* RDKit - version 2023.9.1 .
* numpy - version 1.26.0 .
* pandas - version 2.1.1 .
* scikit-learn - version 1.3.1 .
* lifelines - version 0.27.8 .
* scipy - version 1.11.3 .
* h5py - version 3.9.0 .
* xgboost - version 2.0.2 .


Required information to replicate ICARO<sub>ERT</sub>: an ensemble of Extreme Randomized Trees for Prediction of Protein-Ligand IC<sub>50</sub> is described in this Repository.

### New pIC50 prediction:
A) Feature Normalization Files are within ./ICARO/new_prediction/features/normalization/ folder
B) You can download the trained XTree model in ... 
   Please save this file in ./ICARO/new_prediction/model/
C) Input files:
   For a new prediction, ICARO needs the protein id and fasta and ligand id and smile.
   Protein input example: ./ICARO/new_prediction/protein_sequence_predict.fa
   Ligand input example: ./ICARO/new_prediction/ligand_smiles_predict.smi
D) Script files:
   ./ICARO/new_prediction/icaro_variables.py - has important variables to run ICARO new predictions.
   Please change your working directory in DEFAULT_LOCATION
   ./ICARO/new_prediction/icaro_functions.py - has important functions to run all feature extraction and normalization.
   ./ICARO/new_prediction/fp_admet.sh - script needed to run fp-ADMET descriptors
   ./ICARO/new_prediction/icaro_predict.py - script that retrieves protein and ligand features for a new prediction and retrieves pIC50 predictions
   This script uses as input files: protein_sequence_predict.fa and ligand_smiles_predict.smi

   NOTE: After ifeature and fp-ADMET installation within the ./ICARO/new_prediction/features/ folder, model download in the ./ICARO/new_prediction/model/ folder, you can run a new prediction with your proteins and ligands, by running the command:
   ```bash
   python3 icaro_predict.py
   ```
E) Output file:
   The new predictions will be stored at ./ICARO/new_prediction/new_predictions.csv
   
