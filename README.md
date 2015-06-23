# Enzyme Analysis Pipeline
This repository contains all scripts and data required to reproduce the analysis in the following paper:


I have tried to automate as much of the pipeline as possible. However, since many portions of the pipeline would need to be run on a cluster to be practical, the pipeline is broken into several steps that need to be manually completed. 

## Part I: Site-wise Evolutionary Rates

## Part II: Structural Metrics
1. Clean the data
    - *scripts/clean_pdb.py <input file> <output directory>*
        Begin with biological assemblies for all proteins in the directory ~/data/structures/raw.

