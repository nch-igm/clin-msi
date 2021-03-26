# Clin-MSI

Source code for MSI Detection in tumor samples.

[![PyPI version](https://badge.fury.io/py/clin-msi.svg)](https://badge.fury.io/py/clin-msi)
![Build Status](https://fljy3qhm3i.execute-api.us-east-2.amazonaws.com/prod/?projectName=clin-msi-unittest)

## Installation

### pip

```bash
pip install clin-msi
```

### conda from Bioconda channel

```bash
conda install -c bioconda clin-msi
```

## How to use clin-msi?
### Usage:
    clin-msi <command> [options]

### Key Commands:
* **train**


    Train your model by supplying BAM files with their associated MSI status
    
* **predict**


    Predict the MSI status of a BAM by comparing it to your trained model


# Development

The package is managed using 

## Install Poetry

```bash
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -
source $HOME/.poetry/env
```

## Install projects dependencies

```bash
poetry install -v
```
