#!/bin/bash

module load plink

genotypes="$1"
snps="$2"

plink --bfile ${genotypes} \
      --show-tags ${snps} \
      --tag-r2 0.8 \
      --tag-kb 1000
