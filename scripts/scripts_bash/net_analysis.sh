#!/bin/bash
# Bash launcher for the snakemake pipeline that generates and analyses the networks for the biowide project.
# This script was made using mamba v1.1.0 and conda v22.11.1.

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
REPO_DIR=$SCRIPT_DIR/../..

MAMBA_INSTALLATION_DIR=`mamba info | awk '/base environment/{print $4}'`

ENV_DIR=$REPO_DIR/envs
INPUT_DIR=$REPO_DIR/data
CONFIG_DIR=$REPO_DIR/config
PYTHON_DIR=$REPO_DIR/scripts/scripts_python

source $MAMBA_INSTALLATION_DIR/etc/profile.d/mamba.sh
source $MAMBA_INSTALLATION_DIR/etc/profile.d/conda.sh

if [[ ! -e $ENV_DIR/snakemake ]]; then

    	mamba env create -f $ENV_DIR/snakemake.yml -p $ENV_DIR/snakemake

fi

mamba activate $ENV_DIR/snakemake

cd $REPO_DIR

snakemake -c 15 \
	--use-conda --conda-prefix $ENV_DIR \
	--configfile $CONFIG_DIR/main.yaml \
	-s $PYTHON_DIR/nets.snakefile





