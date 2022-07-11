import numpy as np
import pandas as pd

rule plink_rewrite_keepers:
    input:
	      "data/static/{annotation}/{aggregation}/{type}/individuals.tsv"
    output:
	      "data/static/{annotation}/{aggregation}/{type}/individuals_plink.tsv"
    shell:
	      "awk -v OFS='\t' '{{ $2=$1; print}}' {input} > {output}"

rule plink_expression_reformat:
    input:
        exp="data/static/{annotation}/{aggregation}/{type}/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv"
    output:
        exp="data/static/{annotation}/{aggregation}/{type}/expression.bed.gz"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/mashr/reformat_expression.R"

rule plink_genotype_reformat:
    resources:
        mem_mb=100000
    input:
	      genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
	      inds="data/static/{annotation}/{aggregation}/{type}/individuals_plink.tsv"
    output:
    	  expand("data/static/{{annotation}}/{{aggregation}}/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/static/{annotation}/{aggregation}/{type}/genotypes_filtered_plink"
    shell:
	      "code/mashr/plink_genotype_reformat.sh {input.genotypes} {input.inds} {params.prefix}"