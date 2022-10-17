import numpy as np
import pandas as pd

def list_tensorqtl_outputs(wildcards):
    pb_file = f'data/single_cell_objects/{wildcards.annotation}.{wildcards.aggregation}.tsv'
    samples = list(pd.read_csv(pb_file, sep='\t', nrows=0).columns)[1:]
    pb_clusters = list(np.unique([s.split("_")[1] for s in samples]))
    if wildcards.decomp=="fastgxc":
        pb_clusters.append("Shared")
    return [f"results/static/{wildcards.annotation}/{wildcards.aggregation}/{wildcards.decomp}/{c}/tensorqtl_permutations.tsv"
            for c in pb_clusters]

rule plink_rewrite_keepers:
    input:
        "data/static/{annotation}/{aggregation}/{decomp}/{type}/individuals.tsv"
    output:
        "data/static/{annotation}/{aggregation}/{decomp}/{type}/individuals_plink.tsv"
    shell:
        "awk -v OFS='\t' '{{ $2=$1; print}}' {input} > {output}"

rule plink_genotype_reformat:
    resources:
        mem_mb=100000
    input:
        genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
        inds="data/static/{annotation}/{aggregation}/{decomp}/{type}/individuals_plink.tsv"
    output:
        expand("data/static/{{annotation}}/{{aggregation}}/{{decomp}}/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes_filtered_plink"
    shell:
        "code/mashr/plink_genotype_reformat.sh {input.genotypes} {input.inds} {params.prefix}"
          
rule plink_expression_reformat:
    input:
        exp="data/static/{annotation}/{aggregation}/{decomp}/{type}/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv"
    output:
        exp="data/static/{annotation}/{aggregation}/{decomp}/{type}/expression.bed.gz"
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
        genotypes=expand("data/static/{{annotation}}/{{aggregation}}/{{decomp}}/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/static/{annotation}/{aggregation}/{decomp}/{type}/expression.bed.gz",
        cov="data/static/{annotation}/{aggregation}/{decomp}/{type}/covariates.tsv"
    output:
        expand("results/static/{{annotation}}/{{aggregation}}/{{decomp}}/{{type}}/tensorqtl.cis_qtl_pairs.chr{i}.parquet", i=range(1, 23))
    params:
        plink_prefix="data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes_filtered_plink",
        output_prefix="results/static/{annotation}/{aggregation}/{decomp}/{type}/tensorqtl"
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
        genotypes=expand("data/static/{{annotation}}/{{aggregation}}/{{decomp}}/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/static/{annotation}/{aggregation}/{decomp}/{type}/expression.bed.gz",
        cov="data/static/{annotation}/{aggregation}/{decomp}/{type}/covariates.tsv"
    output:
        "results/static/{annotation}/{aggregation}/{decomp}/{type}/tensorqtl_permutations.tsv"
    params:
        plink_prefix="data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes_filtered_plink"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/mashr/tensorqtl_mtc.py"

rule run_tensorqtl:
    input:
        unpack(list_tensorqtl_outputs)
    output:
        "temp/{annotation}.{aggregation}.{decomp}.tensorqtl.done"
    shell:
        "echo booyah > {output}"
