# Evaluation of MERFISHtools

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.6.1-brightgreen.svg?style=flat-square)](http://snakemake.bitbucket.org)

This Snakemake workflow generates the entire analysis of the forthcoming manuscript
"A Bayesian model for single cell transcript expression analysis on MERFISH data".

## Requirements

Any 64-bit Linux installation with GLIBC 2.5 or newer.

## Setup

### Step 1: Install Miniconda

To run this workflow, you need to install the Conda package manager which is most
easily obtained via the [Miniconda Python distribution](http://conda.pydata.org/miniconda.html).
Miniconda can be installed to your home directory without admin priviledges.

### Step 2: Download the workflow

To download this workflow, clone the git repository into a suitable working directory

    git clone https://github.com/merfishtools/merfishtools-evaluation.git

Then, change to the created directory

    cd merfishtools-evaluation

### Step 3: Install software dependencies

Using the conda package manager, we install the required software into an isolated
environment:

    conda env create --file environment.yaml --name merfishtools-evaluation

The file `environment.yaml` that is shipped with this repository ensures that all
software packages are installed in exactly the same versions as used for the paper,
thereby ensuring reproducibility.

### Step 4: Run the workflow

To run the workflow, you first activate the conda environment created above via

    source activate merfishtools-evaluation

Then you can perfom a dry-run of the workflow with

    snakemake -n

and execute the workflow using, e.g., 8 cores with

    snakemake --cores 8

For further possibilities, see

    snakemake --help

and http://snakemake.bitbucket.org.
