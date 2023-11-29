#!/bin/bash
# DESCRIPTION: Runs fp-admet descriptors for ligand smiles
# NOTE: after fp-admet installation in "./features/" folder with:
# > git clone https://github.com/jcheminform/fpadmet.git
# Please copy this file to the directory "./features/fpadmet/"
# and run: > bash fp_admet.sh

for param in {1..58}
do
  bash runadmet.sh -f ../../ligand_smiles_predict.smi -p $param -a
  mv RESULTS/fps.txt "RESULTS/fps_${param}.txt"
  mv RESULTS/predicted.txt "RESULTS/predicted_${param}.txt"
  echo "Parameter $param processed"
done
