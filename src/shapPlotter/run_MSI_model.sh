#!/bin/bash

#$ -pe smp 8
## The following lines load my own conda environment.
## Fingers crossed, they may still work when another user runs them.
#unset PYTHONPATH
#. /Volumes/igm/home/jbg001/miniconda3/etc/profile.d/conda.sh
#conda activate rnaenv

## Now we put in the parameters. I would set "moddir" t
## The model assumes "infile" is a matrix formatted in Andrei's style.
## The outfile is a tiny csv with one row and two columns ("samp" and "yprob").
## The plotfile gives the shap feature importance plot.
## You can set moddir to wherever you store the directory "mods_071019_frozen"

infile=$1
sample_name=$2
output_dir=$3
script_dir=$4

python "$script_dir/shapPlotter/Main_apply_MSI_model_to_formatted_matrix.py" \
--infile  "$infile" \
--outfile  "$output_dir/$sample_name.csv" \
--plotfile  "$output_dir/$sample_name.shap_plot.pdf" \
--moddir "$script_dir/shapPlotter/mod_071019_frozen"
