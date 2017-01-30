# Evaluation of MERFISHtools

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.10.1-brightgreen.svg)](http://snakemake.bitbucket.org)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.264016.svg)](https://doi.org/10.5281/zenodo.264016)

This Snakemake workflow generates the entire analysis of the forthcoming manuscript
"A Bayesian model for single cell transcript expression analysis on MERFISH data".

## Requirements

* Any 64-bit Linux installation with GLIBC 2.5 or newer (i.e. any Linux distribution that is newer than CentOS 6).

## Setup

### Step 1: Setup Snakemake

To run this workflow, you need to
[setup Snakemake via the Conda package manager](http://snakemake.readthedocs.io/en/latest/getting_started/installation.html).
This does not require admin priviledges.

### Step 2: Download the workflow

To download this workflow, clone the git repository into a suitable working directory

    git clone https://github.com/merfishtools/merfishtools-evaluation.git

Then, change to the created directory

    cd merfishtools-evaluation

### Step 4: Run the workflow

Then you can perfom a dry-run of the workflow with

    snakemake -n

and execute the workflow using, e.g., 24 cores with

    snakemake --cores 24 --use-conda

Note that the argument `--use-conda` is mandatory in such that Snakemake
can deploy the software packages delivered together with the workflow.
For further possibilities, see

    snakemake --help

and http://snakemake.bitbucket.org.

# Author

[Johannes Köster](http://johanneskoester.bitbucket.org)

# License

[MIT](LICENSE.md)
