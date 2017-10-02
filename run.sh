#!/usr/bin/env bash
CONFIG_FILE=$1
OUTPUT_DIR=$2
PATH="$HOME/miniconda3/bin:/APPS/bin/:$PATH"
source $HOME/miniconda3/bin/activate exome
snakemake -pF --cores 10 --snakefile $HOME/Devel/exome/Snakefile --configfile $CONFIG_FILE --unlock
snakemake -pF --cores 10 --snakefile $HOME/Devel/exome/Snakefile --configfile $CONFIG_FILE -d $OUTPUT_DIR
