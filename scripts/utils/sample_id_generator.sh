#!/bin/bash


#==================================================================================================
# Created on: 2017-06-28
# Usage: ./sample_id_generator.sh
# Author: Javier Quilez (GitHub: jaquol)
# Goal: generates unique sample identifier (SAMPLE_ID) based on the samples metadata
#==================================================================================================



#==================================================================================================
# CONFIGURATION VARIABLES AND PATHS
#==================================================================================================

PSL=/users/GR/mb/jquilez/projects/parallel_sequencing_lives

# variables
mode=all_samples

# paths 
url=https://zenodo.org/record/817549/files/metadata.tsv
sample_id_generator=$PSL/scripts/utils/sample_id_generator.py
python=`which python`



#==================================================================================================
# COMMANDS
#==================================================================================================

# *** BEFORE PUBLICATION ***:
# uncomment `wget --no-check-certificate -q -O - $url > $spreadsheet`
# remove `cp data/metadata.tsv $spreadsheet`

# download input metadata from the above specified URL
spreadsheet=$PSL/metadata/metadata.tsv
wget --no-check-certificate -q -O - $url > $spreadsheet

# generate sample IDs
echo
python $sample_id_generator $spreadsheet $mode
echo