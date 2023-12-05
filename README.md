## ICARO<sub>ERT</sub> (IC<sub>50</sub> binding Affinity Regression Optimized)

## Abstract: 
The search for new pharmaceutical products has become increasingly difficult and expensive, leading the pharmaceutical industry to explore alternative methods of drug discovery. One approach involves using computational methods to reduce the number of compounds during the early stages of drug development. Protein-ligand interactions play a crucial role in understanding the mechanisms of drug action. Therefore, there is great interest in developing computational methods for determining drug-target binding affinity in the form of inhibitory concentrations (IC<sub>50</sub>), inhibition constants (Ki), or dissociation constants (K<sub>d</sub>) for hit recognition. Although several predictors have been developed, they mainly focus on K<sub>i</sub>and K<sub>d</sub>, and are restricted to interactions involving kinases. Therefore, a general algorithm that can predict the IC<sub>50</sub> values for different targets is required. Here, we propose ICARO<sub>ERT</sub> (IC<sub>50</sub> binding Affinity Regression Optimized) for drug-target binding affinity prediction, which is expressed as pIC<sub>50</sub>. ICARO<sub>ERT</sub> achieved a coefficient of determination (R<sup>2</sup>) value of 0.726, root mean squared error (RMSE) of 0.701, Pearson correlation coefficient (PCC) of 0.855, and concordance index (CI) of 0.837. ICARO<sub>ERT</sub> is an efficient and cost-effective method for drug discovery.

![Graphical Abstract](Graphical_Abstract.png)⁩

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

* fp-admet - Github available at https://github.com/jcheminform/fpadmet. The fp-admet folder should be downloaded into. /feature/ folder.
  
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

### Study Replication:
A) Dataset files: All data and support files are available in this [link](https://www.dropbox.com/scl/fo/zfb4syff5vpmbdwvfjzzg/h?rlkey=fxk58b0fvcrqw8ktfpyv7ndoj&dl=0). 
   Please, download data files into ./ICARO/data/ folder.

B) Script files:
After performing the changes previously indicated and properly installing and setting up the environment, these scripts should simply run without requiring any changes.

 0) **```icaro_resources.py```** - Includes several variables and functions that are called throughout the pipeline.
 1) **```1_dataset_split.py```** - Dataset split script, creating training, test, and validation sets. Each subset was written as a txt file with unique proteins, ligands, and interactions included in each set. It creates a new ./results/ and ./feature/ folder.

```bash
python3 1_dataset_split.py
```

 2) **```2_FeatureExtraction_ifeature.py```** - Retrieves iFeature protein features. A FASTA file with protein sequences from a dataset exports signature protein features. This requires the installation of an iFeature within . /feature/ folder.

```bash
cd ../features/
git clone https://github.com/Superzchen/iFeature
cd ../scripts/
python3 2_FeatureExtraction_ifeature.py
```
 3) **```2_Normalization_ifeature.py```** - - iFeature feature normalization script that receives a txt file of training, test, and validation proteins and writes an h5 file with iFeature descriptors normalized by the training features.

```bash
python3 2_Normalization_ifeature.py
```

 4) **```3_FeatureExtraction_mordred.py```** - Retrieves MORDRED ligand features. From ligand's smiles available at the "Canonical Smiles" column from the main dataset, it exports MORDRED features in a h5 file.

```bash
python3 3_FeatureExtraction_mordred.py
```

 5) **```4_fp_admet.sh```** - Retrieves fp-ADMET ligand descriptors. You will need to copy this file to the folder ./features/fpadmet/ after instalation and run with
```bash
cd ../features/
git clone https://github.com/jcheminform/fpadmet.git
cd ../scripts/
cp 4_fp_admet.sh ../features/fpadmet/
cd ../features/
bash 4_fp_admet.sh
cd ../scripts/
```

 7) **```4_Normalization_fp_admet.py```** - fp-ADMET descriptors normalization. It receives a txt file of ligands in the training, testing, and validation sets and writes an h5 file with ADMET descriptors normalized by training features.

```bash
python3 4_Normalization_fp_admet.py
```

 8) **```5_FeatureExtraction_prt_family.py```** - Constructs a one-hot encoding feature regarding protein family based on ChEMBL's classification. The script accepts two csv files: one with the protein ID ("single_proteins_id_noduplicate.csv") and other with the ChEMBL's family classification ("targets_family.csv"). It writes a .txt for each family name, and an h5 file to save the corresponding values.

```bash
python3 5_FeatureExtraction_prt_family.py
```

 9)  **```6_model_xtrees.py```** - Trains an Ensemble of Extreme Randomized Trees from h5 feature files from training data and tests them on the test and validation sets.
```bash
python3 6_model_xtrees.py
```

### If you use our predictor, please cite the following.

[Ana T. Gaspar, Catarina Marques-Pereira, António J. Preto and Irina S. Moreira - ICARO: IC<sub>50</sub> binding Affinity Regression Opti-mized] PENDING CITATION
