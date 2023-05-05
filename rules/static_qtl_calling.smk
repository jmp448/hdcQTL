#TODO loop in the analysis files where pseudobulking is actually applied (right now analysis/annotation/assign_cellid, analysis/trajectory_inference/assign_cmstages)
#TODO move all finalized code to the static_qtl_calling folder
#TODO address discrepancy with elorbany data getting SCT counts for dynamic, RNA counts for static
#TODO look into whether identity matrix should actually be included with mash 

import numpy as np
import pandas as pd
from itertools import product

def list_tensorqtl_nominal_outputs(wildcards):
    pb_file = f'data/static_qtl_calling/{wildcards.annotation}/pseudobulk_tmm/{wildcards.annotation}.pseudobulk_tmm.tsv'
    samples = list(pd.read_csv(pb_file, sep='\t', nrows=0).columns)[1:]
    celltypes = list(np.unique([s.split("_")[1] for s in samples]))
    chromosomes = list(range(1, 23))
    return [f"results/static_qtl_calling/{wildcards.annotation}/pseudobulk_tmm/basic/{x[0]}/tensorqtl.cis_qtl_pairs.chr{x[1]}.parquet" for x in product(celltypes, chromosomes)]

def list_tensorqtl_outputs(wildcards):
    pb_file = f'data/static_qtl_calling/{wildcards.annotation}/pseudobulk_tmm/{wildcards.annotation}.pseudobulk_tmm.tsv'
    samples = list(pd.read_csv(pb_file, sep='\t', nrows=0).columns)[1:]
    pb_clusters = list(np.unique([s.split("_")[1] for s in samples]))
    return [f"results/static_qtl_calling/{wildcards.annotation}/pseudobulk_tmm/basic/{c}/8pcs/tensorqtl_permutations.tsv" for c in pb_clusters]

def list_tensorqtl_independent_outputs(wildcards):
    pb_file = f'data/static_qtl_calling/{wildcards.annotation}/pseudobulk_tmm/{wildcards.annotation}.pseudobulk_tmm.tsv'
    samples = list(pd.read_csv(pb_file, sep='\t', nrows=0).columns)[1:]
    pb_clusters = list(np.unique([s.split("_")[1] for s in samples]))
    return [f"results/static_qtl_calling/{wildcards.annotation}/pseudobulk_tmm/basic/{c}/8pcs/tensorqtl_independent.tsv" for c in pb_clusters]

### PSEUDOBULK PREPROCESSING
rule pseudobulk_qc:
    resources:
        mem_mb=250000,
        time="1:30:00"
    input:
        pseudobulk="data/static_qtl_calling/{annotation}/pseudobulk_tmm/{annotation}.pseudobulk_tmm.tsv",
        sample_summary="data/static_qtl_calling/{annotation}/pseudobulk_tmm/sample_summary.tsv"
    output:
        sample_summary_manual="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/sample_summary_manual.tsv"
    params:
        table_prefix = "data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic",
        fig_prefix = "figs/static_qtl_calling/{annotation}/pseudobulk_tmm/basic"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/static_qtl_calling/pseudobulk_tmm-basic-qc.R"
        
rule pseudobulk_agg:
    resources:
        mem_mb=250000,
        time="1:30:00"
    input:
        pseudobulk="data/static_qtl_calling/{annotation}/pseudobulk_tmm/{annotation}.pseudobulk_tmm.tsv",
        sample_summary_manual="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/sample_summary_manual.tsv",
        celltypes="data/static_qtl_calling/{annotation}/pseudobulk_tmm/{annotation}.pseudobulk_tmm.tsv"
    output:
        all_expression="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/pseudobulk_all.tsv"
    params:
        table_prefix = "data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic",
        fig_prefix = "figs/static_qtl_calling/{annotation}/pseudobulk_tmm/basic"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/static_qtl_calling/pseudobulk_tmm-basic-agg.R"

### PREPROCESSING GENOTYPES
rule list_all_individuals:
    input:
	      "data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz"
    output:
    	  "data/genotypes/all_individuals.tsv"
    shell:
	      "code/static_qtl_calling/list_all_individuals.sh {input} {output}"

rule genotype_filter:
    input:
	      genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
	      inds="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/individuals.tsv"
    output:
    	  "data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered.recode.vcf"
    params:
        prefix="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered"
    shell:
	      "code/static_qtl_calling/genotype_filter.sh {input.genotypes} {input.inds} {params.prefix}"

rule compute_af_all:
    resources:
        mem_mb=10000,
        time="30:00"
    input:
        genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
        inds="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/individuals.tsv"
	  output:
	      "data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/af_all.frq"
	  params:
	      prefix="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/af_all"
	  shell:
	      "code/static_qtl_calling/compute_af_all.sh {input.genotypes} {input.inds} {params.prefix}"

rule plink_rewrite_keepers:
    input:
        "data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/individuals.tsv"
    output:
        "data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/individuals_plink.tsv"
    shell:
        "awk -v OFS='\t' '{{ $2=$1; print}}' {input} > {output}"

rule plink_genotype_reformat:
    resources:
        mem_mb=100000
    input:
        genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
        inds="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/individuals_plink.tsv"
    output:
        expand("data/static_qtl_calling/{{annotation}}/pseudobulk_tmm/basic/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered_plink"
    shell:
        "code/static_qtl_calling/plink_genotype_reformat.sh {input.genotypes} {input.inds} {params.prefix}"
          
rule plink_expression_reformat:
    input:
        exp="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv"
    output:
        exp="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/expression.bed.gz"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/static_qtl_calling/reformat_expression.R"
        
### TENSORQTL RULES
rule tensorqtl_nominal:
    resources:
        mem_mb=100000,
        partition="gpu2",
        gres="gpu:1",
        nodes=1
    input:
        genotypes=expand("data/static_qtl_calling/{{annotation}}/pseudobulk_tmm/basic/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/expression.bed.gz",
        cov="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/covariates.tsv"
    output:
        expand("results/static_qtl_calling/{{annotation}}/pseudobulk_tmm/basic/{{type}}/tensorqtl.cis_qtl_pairs.chr{i}.parquet", i=range(1, 23))
    params:
        plink_prefix="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered_plink",
        output_prefix="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/tensorqtl"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/static_qtl_calling/tensorqtl_nominal.py"

rule tensorqtl_merge_nominal:
    resources:
        mem_mb=200000,
        time="30:00"
    input:
        unpack(list_tensorqtl_nominal_outputs)
    output:
        merged_df="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/tensorqtl_nominal.all.tsv",
        beta_df="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/tensorqtl_nominal.betas.tsv",
        se_df="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/tensorqtl_nominal.standard_errors.tsv"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/static_qtl_calling/tensorqtl_merge_nominal.py"
        
rule tensorqtl_permutations:
    resources:
        mem_mb=100000,
        partition="gpu2",
        gres="gpu:1",
        nodes=1,
        time="02:00:00"
    input:
        genotypes=expand("data/static_qtl_calling/{{annotation}}/pseudobulk_tmm/basic/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/expression.bed.gz",
        cov="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/covariates.tsv"
    output:
        cis_df="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}pcs/tensorqtl_permutations.tsv"
    params:
        plink_prefix="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered_plink"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/static_qtl_calling/tensorqtl_permutations.py"

rule tensorqtl_merge:
    resources:
        mem_mb=100000,
        time="30:00"
    input:
        unpack(list_tensorqtl_outputs)
    output:
        all_qtls="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/tensorqtl_permutations.all.tsv",
        top_qtls="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/tensorqtl_permutations.sighits.tsv"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/static_qtl_calling/tensorqtl_merge.R"

rule tensorqtl_fdr:
    resources:
        mem_mb=50000,
        time="30:00"
    input:
        "results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}pcs/tensorqtl_permutations.tsv"
    output:
        "results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}pcs/tensorqtl_permutations.fdr.tsv"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/static_qtl_calling/tensorqtl_fdr.R"

rule tensorqtl_independent:
    resources:
        mem_mb=100000,
        partition="gpu2",
        gres="gpu:1",
        nodes=1,
        time="02:00:00"
    input:
        genotypes=expand("data/static_qtl_calling/{{annotation}}/pseudobulk_tmm/basic/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/expression.bed.gz",
        cov="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/covariates.tsv",
        cis_df="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}pcs/tensorqtl_permutations.fdr.tsv"
    output:
        indep_df="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}pcs/tensorqtl_independent.tsv"
    params:
        plink_prefix="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered_plink"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/static_qtl_calling/tensorqtl_independent.py"
        
rule tensorqtl_merge_independent:
    resources:
        mem_mb=100000,
        time="30:00"
    input:
        unpack(list_tensorqtl_independent_outputs)
    output:
        indep_df="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/tensorqtl_independent.all.tsv"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/static_qtl_calling/tensorqtl_merge_independent.R"
 
rule list_significant_variant_gene_pairs:
    resources:
        mem_mb=100000,
        time="30:00"
    input:
        permutations="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/tensorqtl_permutations.all.tsv",
        nominal="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/tensorqtl_nominal.all.tsv"
    output:
        hit_list="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/signif_variant_gene_pairs.tsv"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/static_qtl_calling/list_significant_tests.R"

### MASH
rule mash_prep:
    resources:
        mem_mb=250000,
        time="06:00:00"
    input:
        beta_df="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/tensorqtl_nominal.betas.tsv",
        se_df="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/tensorqtl_nominal.standard_errors.tsv",
        sample_summary="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/sample_summary_manual.tsv"
    output:
        mash_inputs="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/mash_inputs.Rdata"
    params:
        min_contexts=3
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/mash/mash_prep.R"
        
rule mash:
    resources:
        mem_mb=250000,
        time="05:00:00"
    input:
        mash_inputs="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/mash_inputs.Rdata"
    output:
        trained_model="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/mash_trained_model.rds",
        tophits_fitted_model="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/mash_fitted_model.tophits.rds"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/mash/mash.R"