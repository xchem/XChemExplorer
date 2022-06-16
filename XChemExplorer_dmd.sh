#!/bin/bash

if [ -d "/dls/labxchem" ]
  then
    module load buster/20211020
    module load phenix/1.20
    module load ccp4/7.1.018
    export XChemExplorer_DIR="/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer"
    export BDG_TOOL_MOGUL=/dls_sw/apps/ccdc/CSD_2020/bin/mogul
fi

ccp4-python $XChemExplorer_DIR/XChemExplorer.py
