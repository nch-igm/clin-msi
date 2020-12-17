# IGM Churchill Variant Calling Tools

Source code for variant calling related processes in Churchill.

- Mutect2
- Haplotype Caller
- Genotype process for generating multi-sample VCFs

## Installation

Install as python package

~~~~bash
python -m venv .venv
source .venv/bin/activate
pip install .
~~~~

Build docker image

~~~~bash
docker-compose build
~~~~