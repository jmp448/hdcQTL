#!/bin/bash

module load python

chain_file="$1"
bed_input="$2"
bed_output="$3"

python $CONDA_PREFIX/bin/CrossMap.py bed $chain_file $bed_input $bed_output