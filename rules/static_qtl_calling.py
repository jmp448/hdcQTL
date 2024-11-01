#TODO loop in the analysis files where pseudobulking is actually applied (right now analysis/annotation/assign_cellid, analysis/trajectory_inference/assign_cmstages)
#TODO move all finalized code to the static_qtl_calling folder
#TODO address discrepancy with elorbany data getting SCT counts for dynamic, RNA counts for static
#TODO find more stable software solution than installing flashier from github for mash https://github.com/willwerscheid/flashier

import numpy as np
import pandas as pd
from itertools import product

def list_tensorqtl_nominal_outputs(wildcards):
    pb_file = f'data/static_qtl_calling/{wildcards.annotation}/pseudobulk_tmm/{wildcards.annotation}.pseudobulk_tmm.tsv'
    samples = list(pd.read_csv(pb_file, sep='\t', nrows=0).columns)[1:]
    celltypes = list(np.unique([s.split("_")[1] for s in samples]))
    chromosomes = list(range(1, 23))
    return [f"results/static_qtl_calling/{wildcards.annotation}/pseudobulk_tmm/basic/{x[0]}/{wildcards.npcs}pcs/tensorqtl.cis_qtl_pairs.chr{x[1]}.parquet" for x in product(celltypes, chromosomes)]

def list_tensorqtl_outputs(wildcards):
    pb_file = f'data/static_qtl_calling/{wildcards.annotation}/pseudobulk_tmm/{wildcards.annotation}.pseudobulk_tmm.tsv'
    samples = list(pd.read_csv(pb_file, sep='\t', nrows=0).columns)[1:]
    pb_clusters = list(np.unique([s.split("_")[1] for s in samples]))
    return [f"results/static_qtl_calling/{wildcards.annotation}/pseudobulk_tmm/basic/{c}/{wildcards.npcs}pcs/tensorqtl_permutations.tsv" for c in pb_clusters]

### PSEUDOBULK PREPROCESSING
## Initial pseudobulk aggregation is implemented in analysis/annotation/assign_cellid.ipynb
rule pseudobulk_qc:
    resources:
        mem_mb=10000,
        time="30:00"
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
        
# Pseudobulk QC is implemented in `analysis/static_qtl_calling/pseudobulk_qc.Rmd`

rule pseudobulk_agg:
    resources:
        mem_mb=10000,
        time="30:00"
    input:
        pseudobulk="data/static_qtl_calling/{annotation}/pseudobulk_tmm/{annotation}.pseudobulk_tmm.tsv",
        sample_summary_manual="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/sample_summary_manual.tsv",
        metadata="/project2/gilad/katie/ebQTL/CombinedFormationAndCollectionMetadata_102andPilot_SWAPSANDCONTAMINATIONADDED_012522.csv",
        celltypes="data/static_qtl_calling/{annotation}/pseudobulk_tmm/{annotation}.pseudobulk_tmm.tsv"
    output:
        all_expression="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/pseudobulk_all.tsv"
    params:
        table_prefix = "data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic",
        fig_prefix = "figs/static_qtl_calling/{annotation}/pseudobulk_tmm/basic"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/static_qtl_calling/pseudobulk_tmm-basic-agg.R"

"""
This just gives TMM-normalized gene expression across all cell types, whereas the previous don't actually save
this intermediate (they additionally inverse normal transform the expression)

It performs TMM normalization separately within each cell type, then it aggregates so all
samples are in one dataframe, and it also takes the median across individuals from same cell type
"""
rule pseudobulk_tmm_only:
    resources:
        mem_mb=250000,
        time="30:00"
    input:
        pseudobulk="data/static_qtl_calling/{annotation}/pseudobulk_tmm/{annotation}.pseudobulk_tmm.tsv",
        sample_summary_manual="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/sample_summary_manual.tsv",
        celltypes="data/static_qtl_calling/{annotation}/pseudobulk_tmm/{annotation}.pseudobulk_tmm.tsv"
    output:
        all="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/pseudobulk_tmm_all.tsv",
        median="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/pseudobulk_tmm_median.tsv"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/static_qtl_calling/tmm_normalization_only.R"

rule outlier_detection:
    resources:
        mem_mb=250000,
        time="30:00"
    input:
        pseudobulk="data/static_qtl_calling/{annotation}/pseudobulk_tmm/{annotation}.pseudobulk_tmm.tsv",
        sample_summary="data/static_qtl_calling/{annotation}/pseudobulk_tmm/sample_summary.tsv"
    output:
        dstats="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/d_statistics.tsv"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/static_qtl_calling/outlier_stats.R"

### DATA WRANGLING
#### Genotypes
rule list_study_individuals:
    input: 
        "/project2/gilad/jpopp/ebQTL/data/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/sample_summary.tsv"
    output: 
        "data/genotypes/all_individuals_53.tsv"
    shell: 
        "cat {input} | cut -f 3 | tail -n +2 | sort -u | sed 's/^/NA/' > {output}"
    
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
        inds="data/genotypes/all_individuals_53.tsv"
	  output:
	      "data/genotypes/af_all.frq"
	  params:
	      prefix="data/genotypes/af_all"
	  shell:
	      "code/static_qtl_calling/compute_af_all.sh {input.genotypes} {input.inds} {params.prefix}"

rule compute_af_all_celltypespecific:
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
        
### TENSORQTL 
rule tensorqtl_permutations:
    resources:
        mem_mb=10000,
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
        mem_mb=25000,
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

        
## Nominal analysis (to get all significant eQTLs, not just top per gene)
rule tensorqtl_nominal:
    resources:
        mem_mb=10000,
        partition="gpu2",
        gres="gpu:1",
        nodes=1
    input:
        genotypes=expand("data/static_qtl_calling/{{annotation}}/pseudobulk_tmm/basic/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/expression.bed.gz",
        cov="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/covariates.tsv"
    output:
        expand("results/static_qtl_calling/{{annotation}}/pseudobulk_tmm/basic/{{type}}/{{npcs}}pcs/tensorqtl.cis_qtl_pairs.chr{i}.parquet", i=range(1, 23))
    params:
        plink_prefix="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered_plink",
        output_prefix="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}pcs/tensorqtl"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/static_qtl_calling/tensorqtl_nominal.py"

rule tensorqtl_merge_nominal:
    resources:
        mem_mb=25000,
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
        
# List all significant variant-gene pairs
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
        
rule tensorqtl_summary_to_bed_allsigtests:
    resources:
        mem_mb=30000,
        time="15:00"
    input:
        qtl_summary="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/signif_variant_gene_pairs.tsv",
        bim_file="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/all_celltypes_combined/genotypes_filtered_plink.bim",
        gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/tensorqtl/original/signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/tensorqtl_summary_to_bed_allsigtests.R"

## GSEA
rule gsea:
    resources:
        mem_mb=30000,
        time="15:00"
    input:
        eb_bed="results/static_eqtl_followup/qtl_sets/tensorqtl/original/signif_variant_gene_pairs.bed",
        gtex_bed="results/static_eqtl_followup/qtl_sets/tensorqtl/original/signif_variant_gene_pairs.all_tissue_overlap.bed",
        gmt="data/gene_sets/c5.go.bp.v2022.1.Hs.symbols.gmt"
    output:
        gsea_results="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/gsea_results.tsv"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/gsea.R"


# rule tensorqtl_summary_to_bed_alltests:
#     resources:
#         mem_mb=30000,
#         time="15:00"
#     input:
#         qtl_summary="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/tensorqtl_nominal.all.tsv",
#         bim_file="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/PNS-glia/genotypes_filtered_plink.bim",
#         gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
#     output:
#         bedfile="results/static_eqtl_followup/qtl_sets/tensorqtl/original/tensorqtl-all_variant_gene_pairs.bed"
#     conda: "../slurmy/r-mashr.yml"
#     script:
#         "../code/static_eqtl_followup/tensorqtl_summary_to_bed_alltests.R"
        
# rule tensorqtl_fdr:
#     resources:
#         mem_mb=50000,
#         time="30:00"
#     input:
#         "results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}pcs/tensorqtl_permutations.tsv"
#     output:
#         "results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}pcs/tensorqtl_permutations.fdr.tsv"
#     conda:
#         "../slurmy/r-mashr.yml"
#     script:
#         "../code/static_qtl_calling/tensorqtl_fdr.R"
        
# rule list_all_individuals:
#     input:
# 	      "data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz"
#     output:
#     	  "data/genotypes/all_individuals_120.tsv"
#     shell:
# 	      "code/static_qtl_calling/list_all_individuals.sh {input} {output}"
