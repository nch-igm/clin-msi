# Clin-MSI

Source code for MSI Detection in tumor samples.

## Installation

Install as python package

~~~~bash
python -m venv .venv
source .venv/bin/activate
pip install .
~~~~

Install using Bioconda

    conda install -c bioconda clin-msi

## How to use clin-msi?
### Usage:
    msisensor-pro <command> [options]

### Key Commands:
* **train**


    Train your model by supplying BAM files with their associated MSI status
    
* **predict**


    Predict the MSI status of a BAM by comparing it to your trained model
