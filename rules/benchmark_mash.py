#TODO make assignment line up here, especially since the data inputs are currently just manually duplicated from benchmarking on the static part
#TODO that assignment file currently saves stuff to an unused 'static' subdirectory in 'data'
#TODO handle the fact that TreeBH has to be manually added to R environment (I added it to the mash package)
#TODO consider making qtl master loc be created in a separate file, not with the treebh 
#TODO I had to manually install flashier into my conda environment, fix that https://rdrr.io/github/willwerscheid/flashier/

import numpy as np
import pandas as pd
from itertools import product

### HELPERS
def get_decomp_type(dc):
    return dc.split('_')[0]

def list_tensorqtl_nominal_outputs(wildcards):
    pb_file = 'data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/eb_cellid.pseudobulk_tmm.tsv'
    samples = list(pd.read_csv(pb_file, sep='\t', nrows=0).columns)[1:]
    celltypes = list(np.unique([s.split("_")[1] for s in samples]))
    chromosomes = list(range(1, 23))
    return [f"results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{wildcards.decomp}/{x[0]}/{wildcards.tensorqtl}.cis_qtl_pairs.chr{x[1]}.parquet" for x in product(celltypes, chromosomes)]

### PSEUDOBULK PREPROCESSING
rule pseudobulk_qc_benchmark_specificity:
    resources:
        mem_mb=250000,
        time="1:30:00"
    input:
        pseudobulk="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/eb_cellid.pseudobulk_tmm.tsv",
        sample_summary="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/sample_summary.tsv"
    output:
        sample_summary_manual="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/sample_summary_manual.tsv"
    params:
        table_prefix = "data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}",
        fig_prefix = "figs/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}",
        decomp_type= get_decomp_type("{decomp}")
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/benchmark_specificity_methods/pseudobulk_tmm-{wildcards.decomp}-qc.R"
        
rule pseudobulk_agg_benchmark_specificity:
    resources:
        mem_mb=250000,
        time="1:30:00"
    input:
        pseudobulk="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/eb_cellid.pseudobulk_tmm.tsv",
        sample_summary_manual="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/sample_summary_manual.tsv",
        celltypes="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/eb_cellid.pseudobulk_tmm.tsv"
    output:
        all_expression="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/pseudobulk_all.tsv"
    params:
        table_prefix = "data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}",
        fig_prefix = "figs/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/benchmark_specificity_methods/pseudobulk_tmm-{wildcards.decomp}-agg.R"

rule plink_expression_reformat_benchmark_specificity:
    input:
        exp="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv"
    output:
        exp="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/expression.bed.gz"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods/reformat_expression.R"
        
### PREPROCESSING GENOTYPES
rule plink_rewrite_keepers_benchmark_specificity:
    input:
        "data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/individuals.tsv"
    output:
        "data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/individuals_plink.tsv"
    shell:
        "awk -v OFS='\t' '{{ $2=$1; print}}' {input} > {output}"

rule plink_genotype_reformat_benchmark_specificity:
    resources:
        mem_mb=100000
    input:
        genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
        inds="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/individuals_plink.tsv"
    output:
        expand("data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{{decomp}}/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes_filtered_plink"
    shell:
        "code/benchmark_specificity_methods/plink_genotype_reformat.sh {input.genotypes} {input.inds} {params.prefix}"

### TENSORQTL RULES
rule tensorqtl_nominal_benchmark_specificity:
    resources:
        mem_mb=100000,
        partition="gpu2",
        gres="gpu:1",
        nodes=1
    input:
        genotypes=expand("data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{{decomp}}/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/expression.bed.gz",
        cov="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/covariates.tsv"
    output:
        expand("results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{{decomp}}/{{type}}/{{tensorqtl}}.cis_qtl_pairs.chr{i}.parquet", i=range(1, 23))
    params:
        plink_prefix="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes_filtered_plink",
        output_prefix="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/{tensorqtl}"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/benchmark_specificity_methods/{wildcards.tensorqtl}_nominal.py"

rule tensorqtl_merge_nominal_benchmark_specificity:
    resources:
        mem_mb=200000,
        time="03:00:00"
    input:
        unpack(list_tensorqtl_nominal_outputs)
    output:
        merged_df="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_nominal.all.tsv",
        beta_df="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_nominal.betas.tsv",
        se_df="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_nominal.standard_errors.tsv"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/benchmark_specificity_methods/tensorqtl_merge_nominal.py"

"""       
rule tensorqtl_permutations:
    resources:
        mem_mb=100000,
        partition="gpu2",
        gres="gpu:1",
        nodes=1,
        time="02:00:00"
    input:
        genotypes=expand("data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/expression.bed.gz",
        cov="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/covariates.tsv"
    output:
        cis_df="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/{npcs}pcs/tensorqtl_permutations.tsv"
    params:
        plink_prefix="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes_filtered_plink"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/benchmark_specificity_methods/tensorqtl_permutations.py"

rule tensorqtl_merge:
    resources:
        mem_mb=100000,
        time="30:00"
    input:
        unpack(list_tensorqtl_outputs)
    output:
        all_qtls="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/tensorqtl_permutations.all.tsv",
        top_qtls="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/tensorqtl_permutations.sighits.tsv"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods/tensorqtl_merge.R"

rule tensorqtl_fdr:
    resources:
        mem_mb=50000,
        time="30:00"
    input:
        "results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/{npcs}pcs/tensorqtl_permutations.tsv"
    output:
        "results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/{npcs}pcs/tensorqtl_permutations.fdr.tsv"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods/tensorqtl_fdr.R"
"""

### MASH
rule mash_train_test_split:
    resources:
        partition="bigmem2",
        mem_mb=250000,
        time="06:00:00"
    input:
        beta_df="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_nominal.betas.tsv",
        se_df="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_nominal.standard_errors.tsv",
        gene_locs="/project2/gilad/jpopp/ebQTL/data/gencode/gencode.hg38.filtered.tss.tsv"
    output:
        mash_inputs="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_mash_inputs.train_test_split.Rdata"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods/mashr_train_test_split.R"
        
rule mash_train_test_split_datadriven:
    resources:
        partition="gilad",
        mem_mb=250000,
        time="06:00:00"
    input:
        beta_df="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_nominal.betas.tsv",
        se_df="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_nominal.standard_errors.tsv",
        gene_locs="/project2/gilad/jpopp/ebQTL/data/gencode/gencode.hg38.filtered.tss.tsv"
    output:
        mash_inputs="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_mash_inputs.train_test_split_datadriven.Rdata"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods/mashr_train_test_split_datadriven.R"

        
rule mashr_corr_eval:
    resources:
        partition="bigmem2",
        mem_mb=250000,
        time="05:00:00"
    input:
        data="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_mash_inputs.train_test_split.Rdata"
    output:
        "results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_mash_corr_eval.Rdata"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods/mashr_correlation_evaluation.R"

rule mashr_alpha_eval:
    resources:
        partition="bigmem2",
        mem_mb=250000,
        time="05:00:00"
    input:
        data="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_mash_inputs.train_test_split.Rdata"
    output:
        "results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_mash_alpha_eval.Rdata"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods/mashr_alpha_evaluation.R"

rule mashr_U_eval:
    resources:
        mem_mb=250000,
        time="05:00:00"
    input:
        data="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_mash_inputs.train_test_split_datadriven.Rdata"
    output:
        "results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/{tensorqtl}_mash_U_eval.Rdata"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods/mashr_U_evaluation.R"

rule mash_benchmark_specificity:
    resources:
        mem_mb=250000,
        time="05:00:00"
    input:
        mash_inputs="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/mash_inputs.Rdata"
    output:
        trained_model="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/mash_trained_model.rds",
        tophits_fitted_model="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{npcs}pcs/mash_fitted_model.tophits.rds"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/mash/mash.R"
