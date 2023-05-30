#TODO loop in the analysis files where pseudobulking is actually applied
#TODO move tensorqtl list all cell types
#TODO

import numpy as np
import pandas as pd
from itertools import product

### HELPERS
def get_prefix(path, c=1):
    return '/'.join(path.split('/')[:-c])
    
def get_normalization_type(dc):
    return dc.split('_')[0]

def list_matrixeqtl_outputs(wildcards, npc="8pcs"):
    datasets=["vqtl_ipsc", "ebqtl_ipsc", "ebqtl_ipscnofilt"]
    aggregations=["pseudobulk_tmm", "pseudobulk_scran"]
    normalizations=["basic", "quantile"]
    errormodels=["homoscedastic", "heteroscedastic"]
    
    combos = list(product(datasets, aggregations, normalizations, errormodels))
    
    file_list = [f"results/benchmark_static_qtl_calling/{x[0]}/{x[1]}/{x[2]}/IPSC/{npc}/matrixeqtl_{x[3]}.cis_qtl_pairs.all.mtc.tsv" for x in combos]
    
    return file_list
    
def list_subsampling_outputs(wildcards, npc="5pcs"):
    datasets=["ebqtl_40sub", "ebqtl_30sub", "ebqtl_20sub"]
    file_list = [f"results/benchmark_static_qtl_calling/{s}/pseudobulk_tmm/basic/IPSC/{npc}/matrixeqtl_homoscedastic.cis_qtl_pairs.all.mtc.tsv" for s in datasets]
    return file_list

def list_tensorqtl_outputs(wildcards):
    pb_file = f'data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/eb_cellid.pseudobulk_tmm.tsv'
    samples = list(pd.read_csv(pb_file, sep='\t', nrows=0).columns)[1:]
    pb_clusters = list(np.unique([s.split("_")[1] for s in samples]))
    return [f"results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/basic/{c}/8pcs/tensorqtl_permutations.tsv"
            for c in pb_clusters]

### PSEUDOBULK PREPROCESSING
rule pseudobulk_qc_benchmark:
    resources:
        mem_mb=250000,
        time="1:30:00"
    input:
        pseudobulk="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{dataset}.{aggregation}.tsv",
        sample_summary="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/sample_summary.tsv"
    output:
        sample_summary_manual="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/sample_summary_manual.tsv"
    params:
        table_prefix = "data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}",
        fig_prefix = "figs/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}",
        normalization_type= get_normalization_type("{normalization}")
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/benchmarking_static_qtl_calling/{wildcards.aggregation}-{wildcards.normalization}-qc.R"
        
rule pseudobulk_agg_benchmark:
    resources:
        mem_mb=250000,
        time="1:30:00"
    input:
        pseudobulk="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{dataset}.{aggregation}.tsv",
        sample_summary_manual="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/sample_summary_manual.tsv",
        celltypes="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{dataset}.{aggregation}.tsv"
    output:
        all_expression="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/pseudobulk_all.tsv"
    params:
        table_prefix = "data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}",
        fig_prefix = "figs/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/benchmarking_static_qtl_calling/{wildcards.aggregation}-{wildcards.normalization}-agg.R"

### PREPROCESSING GENOTYPES
rule genotype_filter_benchmark:
    input:
	      genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
	      inds="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/individuals.tsv"
    output:
    	  "data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes_filtered.recode.vcf"
    params:
        prefix="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes_filtered"
    shell:
	      "code/benchmarking_static_qtl_calling/genotype_filter.sh {input.genotypes} {input.inds} {params.prefix}"

rule genotype_012_benchmark:
    input:
	      "data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes_filtered.recode.vcf"
    output:
    	  expand("data/benchmark_static_qtl_calling/{{dataset}}/{{aggregation}}/{{normalization}}/{{type}}/genotypes_filtered.{out}", out=['012', '012.indv', '012.pos'])
    params:
        prefix="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes_filtered"
    shell:
	      "code/benchmarking_static_qtl_calling/genotype_012.sh {input} {params.prefix}"

rule genotype_transpose_benchmark:
    resources:
        mem_mb=50000
    input:
	      "data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes_filtered.012"
    output:
	      "data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes_filtered.012.transpose"
    shell:
	      "code/benchmarking_static_qtl_calling/genotype_transpose.sh {input} {output}"

rule genotype_reformat_benchmark:
    resources:
        mem_mb=50000
    input:
        genotypes="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes_filtered.012.transpose",
        individuals="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes_filtered.012.indv",
        snp_locs="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes_filtered.012.pos"
    output:
        snp_locs="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/snp_locs.tsv",
        genotypes="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes.tsv"
    params:
        temp_loc="temp/genotype_reformat.{dataset}.{aggregation}.{normalization}.{type}"
    shell:
        "code/benchmarking_static_qtl_calling/genotype_reformat.sh {input.genotypes} {input.individuals} {input.snp_locs} {params.temp_loc} {output.snp_locs} {output.genotypes}" 

### PREPROCESSING JUST FOR TENSORQTL
rule plink_rewrite_keepers_benchmark:
    input:
        "data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/individuals.tsv"
    output:
        "data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/individuals_plink.tsv"
    shell:
        "awk -v OFS='\t' '{{ $2=$1; print}}' {input} > {output}"

rule plink_genotype_reformat_benchmark:
    resources:
        mem_mb=100000
    input:
        genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
        inds="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/individuals_plink.tsv"
    output:
        expand("data/benchmark_static_qtl_calling/{{dataset}}/{{aggregation}}/{{normalization}}/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes_filtered_plink"
    shell:
        "code/benchmarking_static_qtl_calling/plink_genotype_reformat.sh {input.genotypes} {input.inds} {params.prefix}"
          
rule plink_expression_reformat_benchmark:
    input:
        exp="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv"
    output:
        exp="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/expression.bed.gz"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/benchmarking_static_qtl_calling/reformat_expression.R"
        
### MATRIX EQTL RULES
rule matrixeqtl_nominal_benchmark:
    resources:
        mem_mb=75000,
        time="00:30:00"
    input:
        genotypes="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes.tsv",
        snp_locs="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/snp_locs.tsv",
        expression="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv",
        covariates="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/covariates.tsv",
        sample_summary="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/sample_summary_manual.tsv"
    output:
        eqtls="results/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/{npcs}/matrixeqtl_{errormodel}.cis_qtl_pairs.all.tsv",
        df="results/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/{npcs}/matrixeqtl_{errormodel}.df.tsv"
    conda: "../slurmy/r-matrixEQTL.yml"
    script:
        "../code/benchmarking_static_qtl_calling/matrixeqtl_nominal_{wildcards.errormodel}.R"

rule matrixeqtl_mtc_benchmark:
    resources:
        mem_mb=75000,
        time="00:15:00"
    input:
        eqtls="results/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/{npcs}/matrixeqtl_{errormodel}.cis_qtl_pairs.all.tsv",
        df="results/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/{npcs}/matrixeqtl_{errormodel}.df.tsv"
    output:
        all_tests="results/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/{npcs}/matrixeqtl_{errormodel}.cis_qtl_pairs.all.mtc.tsv",
        top_tests="results/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/{npcs}/matrixeqtl_{errormodel}.cis_qtl_pairs.tophits.tsv",
        n_hits="results/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/{npcs}/matrixeqtl_{errormodel}.cis_qtl_pairs.nhits.tsv"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/benchmarking_static_qtl_calling/matrixeqtl_mtc.R"

### TENSORQTL RULES
rule tensorqtl_nominal_benchmark:
    resources:
        mem_mb=100000,
        partition="gpu2",
        gres="gpu:1",
        nodes=1
    input:
        genotypes=expand("data/benchmark_static_qtl_calling/{{dataset}}/{{aggregation}}/{{normalization}}/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/expression.bed.gz",
        cov="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/covariates.tsv"
    output:
        expand("results/benchmark_static_qtl_calling/{{dataset}}/{{aggregation}}/{{normalization}}/{{type}}/tensorqtl.cis_qtl_pairs.chr{i}.parquet", i=range(1, 23))
    params:
        plink_prefix="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes_filtered_plink",
        output_prefix="results/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/tensorqtl"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/benchmarking_static_qtl_calling/tensorqtl_nominal.py"
        
rule tensorqtl_mtc_benchmark:
    resources:
        mem_mb=100000,
        partition="gpu2",
        gres="gpu:1",
        nodes=1,
        time="02:00:00"
    input:
        genotypes=expand("data/benchmark_static_qtl_calling/{{dataset}}/{{aggregation}}/{{normalization}}/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/expression.bed.gz",
        cov="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/covariates.tsv"
    output:
        "results/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/{npcs}/tensorqtl_permutations.tsv"
    params:
        plink_prefix="data/benchmark_static_qtl_calling/{dataset}/{aggregation}/{normalization}/{type}/genotypes_filtered_plink"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/benchmarking_static_qtl_calling/tensorqtl_mtc.py"

rule run_tensorqtl_benchmark:
    input:
        unpack(list_tensorqtl_outputs)
    output:
        "results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/basic/8pcs/tensorqtl.cis_qtl_pairs.tophits.tsv"
    conda:
        "../slurmy/mashr.yml"
    script:
        "../code/benchmarking_static_qtl_calling/tensorqtl_merge.R"

### SNAKEMAKE HELPERS 
ruleorder: run_matrixeqtl_subsampling_benchmark > run_matrixeqtl_npcs_benchmark

rule run_matrixeqtl_npcs_benchmark:
    input:
        unpack(list_matrixeqtl_outputs)
    output:
        "temp/benchmark_static_qtl_calling.{npcs}.done"
    shell:
        "echo booyah > {output}"
        
rule run_matrixeqtl_subsampling_benchmark:
    input:
        unpack(list_subsampling_outputs)
    output:
        "temp/benchmark_static_qtl_calling.{npc}.subsampling.done"
    shell:
        "echo booyah > {output}"
