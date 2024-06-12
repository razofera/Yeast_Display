#!/bin/bash
# Setup

# create new conda enviroment
conda create --name yeast_screen
conda activate yeast_screen
conda install -c conda-forge biopython
conda deactivate

