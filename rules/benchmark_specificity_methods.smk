#TODO make assignment line up here, especially since the data inputs are currently just manually duplicated from benchmarking on the static part
#TODO that assignment file currently saves stuff to an unused 'static' subdirectory in 'data'
#TODO handle the fact that TreeBH has to be manually added to R environment (I added it to the mash package)
#TODO consider making qtl master loc be created in a separate file, not with the treebh 

import numpy as np
import pandas as pd
from itertools import product

### HELPERS
def get_decomp_type(dc):
    return dc.split('_')[0]
    
def list_matrixeqtl_outputs(wildcards, npc=8):
    pseudobulk_file = f'data/benchmark_specificity_methods/eb_cellid/eb_cellid.pseudobulk_tmm.tsv'
    samples = list(pd.read_csv(pseudobulk_file, sep='\t', nrows=0).columns)[1:]
    celltypes = list(np.unique([s.split("_")[1] for s in samples]))
    file_list = [f"results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/basic/{ct}/{npc}pcs/matrixeqtl.cis_qtl_pairs.all.tsv" for ct in celltypes]
    return file_list
    
### PSEUDOBULK
    
rule pseudobulk_qc_benchmarking:
    resources:
        mem_mb=250000,
        time="1:30:00"
    input:
        pseudobulk="data/benchmark_specificity_methods/eb_cellid/eb_cellid.pseudobulk_tmm.tsv",
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
        
rule pseudobulk_agg_benchmarking:
    resources:
        mem_mb=250000,
        time="1:30:00"
    input:
        pseudobulk="data/benchmark_specificity_methods/eb_cellid/eb_cellid.pseudobulk_tmm.tsv",
        sample_summary_manual="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/sample_summary_manual.tsv",
        celltypes="data/benchmark_specificity_methods/eb_cellid/eb_cellid.pseudobulk_tmm.tsv"
    output:
        all_expression="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/pseudobulk_all.tsv"
    params:
        table_prefix = "data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}",
        fig_prefix = "figs/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/benchmark_specificity_methods/pseudobulk_tmm-{wildcards.decomp}-agg.R"

rule genotype_filter_benchmarking:
    input:
	      genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
	      inds="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/individuals.tsv"
    output:
    	  "data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes_filtered.recode.vcf"
    params:
        prefix="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes_filtered"
    shell:
	      "code/benchmark_specificity_methods/genotype_filter.sh {input.genotypes} {input.inds} {params.prefix}"

rule genotype_012_benchmarking:
    input:
	      "data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes_filtered.recode.vcf"
    output:
    	  expand("data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{{decomp}}/{{type}}/genotypes_filtered.{out}", out=['012', '012.indv', '012.pos'])
    params:
        prefix="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes_filtered"
    shell:
	      "code/benchmark_specificity_methods/genotype_012.sh {input} {params.prefix}"

rule genotype_transpose_benchmarking:
    resources:
        mem_mb=50000
    input:
	      "data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes_filtered.012"
    output:
	      "data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes_filtered.012.transpose"
    shell:
	      "code/benchmark_specificity_methods/genotype_transpose.sh {input} {output}"

rule genotype_reformat_benchmarking:
    resources:
        mem_mb=50000
    input:
        genotypes="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes_filtered.012.transpose",
        individuals="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes_filtered.012.indv",
        snp_locs="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes_filtered.012.pos"
    output:
        snp_locs="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/snp_locs.tsv",
        genotypes="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes.tsv"
    params:
        temp_loc="temp/genotype_reformat.eb_cellid.pseudobulk_tmm.{decomp}.{type}"
    shell:
        "code/benchmark_specificity_methods/genotype_reformat.sh {input.genotypes} {input.individuals} {input.snp_locs} {params.temp_loc} {output.snp_locs} {output.genotypes}" 

rule matrixeqtl_nominal_benchmarking:
    resources:
        mem_mb=75000,
        time="00:30:00"
    input:
        genotypes="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/genotypes.tsv",
        snp_locs="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/snp_locs.tsv",
        expression="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv",
        covariates="data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/covariates.tsv"
    output:
        eqtls="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.all.tsv",
        df="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/{npcs}/matrixeqtl.df.tsv"
    conda: "../slurmy/r-matrixEQTL.yml"
    script:
        "../code/benchmark_specificity_methods/matrixeqtl_nominal.R"

rule run_matrixeqtl_all:
    input:
        unpack(list_matrixeqtl_outputs)
    output:
        "temp/eb_cellid.pseudobulk_tmm.basic.matrixeqtl.{npcs}pcs.done"
    shell:
        "echo booyah > {output}"

rule tree_bh_benchmark:
    resources:
        mem_mb=150000,
        time="60:00:00"
    input:
        unpack(list_matrixeqtl_outputs)
    output:
        all_qtls="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/all/{npcs}/matrixeqtl.cis_qtl_pairs.{hierarchy}.tsv",
        qtl_significance="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/all/{npcs}/matrixeqtl.cis_qtl_pairs.{hierarchy}.treebh_significance.tsv"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods/treebh.R"
        
"""
rule mashr:
    resources:
        mem_mb=60000,
        time="05:00:00",
        partition="broadwl"
    input:
        unpack(list_matrixeqtl_outputs_npcs)
    output:
        trained="results/benchmark_specificity_methods/{annotation}/pseudobulk_tmm/{decomp}/{npcs}/mashr.trained.rds",
        tophits="results/benchmark_specificity_methods/{annotation}/pseudobulk_tmm/{decomp}/{npcs}/mashr.tophits.rds",
        random="results/benchmark_specificity_methods/{annotation}/pseudobulk_tmm/{decomp}/{npcs}/mashr.random_tests.tsv",
        top="results/benchmark_specificity_methods/{annotation}/pseudobulk_tmm/{decomp}/{npcs}/mashr.top_hits.tsv",
        datasets="results/benchmark_specificity_methods/{annotation}/pseudobulk_tmm/{decomp}/{npcs}/mashr.training_data.RData",
        all_qtls="temp/{annotation}.pseudobulk_tmm.{decomp}.{npcs}.qtls_combined.rds"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods/mashr.R"
"""

rule mashr_corr_eval:
    resources:
        mem_mb=75000,
        time="05:00:00"
    input:
        data="results/benchmark_specificity_methods/highpass_cellid_all/pseudobulk-scran/basic/mashr.training_data.RData",
        gene_locs="/project2/gilad/jpopp/ebQTL/data/gencode/gencode.hg38.filtered.tss.tsv"
    output:
        "results/benchmark_specificity_methods/highpass_cellid_all/pseudobulk-scran/basic/mashr_corr_eval.RData"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods/mashr_correlation_evaluation.R"
        
rule ashr:
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
        eqtls="results/benchmark_specificity_methods/{annotation}/pseudobulk_tmm/{decomp}/Shared/matrixeqtl.cis_qtl_pairs.all.mtc.tsv",
        random="results/benchmark_specificity_methods/{annotation}/pseudobulk_tmm/{decomp}/mashr.random_tests.tsv",
        top="results/benchmark_specificity_methods/{annotation}/pseudobulk_tmm/{decomp}/mashr.top_hits.tsv",
        df="results/benchmark_specificity_methods/{annotation}/pseudobulk_tmm/{decomp}/Shared/matrixeqtl.df.tsv"
    output:
        trained="results/benchmark_specificity_methods/{annotation}/pseudobulk_tmm/{decomp}/ashr.trained.rds",
        tophits="results/benchmark_specificity_methods/{annotation}/pseudobulk_tmm/{decomp}/ashr.tophits.rds"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods/ashr.R"

"""
## Maybe graveyard
rule matrixeqtl_mtc:
    resources:
        mem_mb=75000,
        time="00:15:00"
    input:
        eqtls="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.all.tsv",
        df="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/{npcs}/matrixeqtl.df.tsv"
    output:
        all_tests="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.all.mtc.tsv",
        top_tests="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.tophits.tsv",
        n_hits="results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.nhits.tsv"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods/matrixeqtl_mtc.R"

rule run_matrixeqtl_allpcs:
    input:
        unpack(list_matrixeqtl_outputs)
    output:
        "temp/{annotation}.pseudobulk_tmm.{decomp}.matrixeqtl.done"
    shell:
        "echo booyah > {output}"
        

def list_matrixeqtl_outputs_npcs(wildcards):
    pseudobulk_file = f'data/single_cell_objects/{wildcards.annotation}.{wildcards.aggregation}.tsv'
    samples = list(pd.read_csv(pseudobulk_file, sep='\t', nrows=0).columns)[1:]
    celltypes = list(np.unique([s.split("_")[1] for s in samples]))
    file_list = [f"results/benchmark_specificity_methods/{wildcards.annotation}/{wildcards.aggregation}/{wildcards.decomp}/{c}/{wildcards.npcs}/matrixeqtl.cis_qtl_pairs.all.mtc.tsv"
          for c in celltypes]
    return file_list

def get_prefix(path, c=1):
    return '/'.join(path.split('/')[:-c])
"""

