import numpy as np
import pandas as pd

rule plink_rewrite_keepers:
    input:
        "data/static/{annotation}/{aggregation}/{type}/individuals.tsv"
    output:
        "data/static/{annotation}/{aggregation}/{type}/individuals_plink.tsv"
    shell:
        "awk -v OFS='\t' '{{ $2=$1; print}}' {input} > {output}"

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

rule tensorqtl_nominal:
    resources:
        mem_mb=100000,
        partition="gpu2",
        gres="gpu:1",
        nodes=1
    input:
        genotypes=expand("data/static/{{annotation}}/{{aggregation}}/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/static/{annotation}/{aggregation}/{type}/expression.bed.gz",
        cov="data/static/{annotation}/{aggregation}/{type}/covariates.tsv"
    output:
        expand("data/static/{{annotation}}/{{aggregation}}/{{type}}/tensorqtl.cis_qtl_pairs.chr{i}.parquet", i=range(1, 23))
    params:
        plink_prefix="data/static/{annotation}/{aggregation}/{type}/genotypes_filtered_plink",
        output_prefix="data/static/{annotation}/{aggregation}/{type}/tensorqtl"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/mashr/tensorqtl_nominal.py"
        
rule tensorqtl_mtc:
    resources:
        mem_mb=100000,
        partition="gpu2",
        gres="gpu:1",
        nodes=1
    input:
        genotypes=expand("data/static/{{annotation}}/{{aggregation}}/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/static/{annotation}/{aggregation}/{type}/expression.bed.gz",
        cov="data/static/{annotation}/{aggregation}/{type}/covariates.tsv"
    output:
        "data/static/{annotation}/{aggregation}/{type}/tensorqtl_permutations.tsv"
    params:
        plink_prefix="data/static/{annotation}/{aggregation}/{type}/genotypes_filtered_plink"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/mashr/tensorqtl_mtc.py"
