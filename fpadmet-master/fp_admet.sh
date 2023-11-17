#!/bin/bash
features=(27 31 37 39 40 43 44 55)
for num in {1..58} #${features[@]}
do
  bash runadmet.sh -f smiles.smi -p $num -a #
done
