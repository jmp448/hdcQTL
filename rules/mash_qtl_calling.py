#TODO loop in the analysis files where pseudobulking is actually applied (right now analysis/annotation/assign_cellid, analysis/trajectory_inference/assign_cmstages)
#TODO move all finalized code to the static_qtl_calling folder
#TODO address discrepancy with elorbany data getting SCT counts for dynamic, RNA counts for static
#TODO find more stable software solution than installing flashier from github for mash https://github.com/willwerscheid/flashier

import numpy as np
import pandas as pd
from itertools import product

def list_bims(wildcards):
    pb_file = f'data/static_qtl_calling/{wildcards.annotation}/pseudobulk_tmm/{wildcards.annotation}.pseudobulk_tmm.tsv'
    samples = list(pd.read_csv(pb_file, sep='\t', nrows=0).columns)[1:]
    celltypes = list(np.unique([s.split("_")[1] for s in samples]))
    return [f"data/static_qtl_calling/{{annotation}}/pseudobulk_tmm/basic/{c}/genotypes_filtered_plink.bim" for c in celltypes]

### MASH
rule mash_prep:
    resources:
        mem_mb=20000,
        time="01:00:00"
    input:
        beta_df="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/tensorqtl_nominal.betas.tsv",
        se_df="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/tensorqtl_nominal.standard_errors.tsv",
        sample_summary="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/sample_summary_manual.tsv"
    output:
        mash_inputs="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/mash_inputs.Rdata"
    params:
        min_contexts=10
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/mash_qtl_calling/mash_prep.R"
        
rule mash:
    resources:
        mem_mb=10000,
        time="01:00:00"
    input:
        mash_inputs="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/mash_inputs.Rdata"
    output:
        trained_model="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/mash_trained_model.rds",
        tophits_fitted_model="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/mash_fitted_model.tophits.rds"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/mash_qtl_calling/mash.R"

rule mash_full:
    resources:
        mem_mb=100000,
        time="06:00:00"
    input:
        mash_inputs="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/mash_inputs.Rdata",
        trained_model="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/mash_trained_model.rds"
    output:
        full_fitted_model="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/mash_fitted_model.full.rds"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/mash_qtl_calling/mash_full.R"

rule merge_bims_mash:
    resources:
        mem_mb=50000,
        disk_mb=50000
    input:
        unpack(list_bims)
    output:
        "data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/all_celltypes_combined/genotypes_filtered_plink.bim"
    shell:
        """
        cat {input} | sort -u -k 1,1 -k 4,4n > {output}
        """

rule mash_to_bed_allsigtests:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        mash_model="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/mash_fitted_model.full.rds",
        bim_file="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/all_celltypes_combined/genotypes_filtered_plink.bim",
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/mash/original/mash-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/mash_qtl_calling/mash_to_bed_allsigtests.R"


