#TODO analyze different weighting options
#TODO think about whether the inverse normal transformation's tie-breaking approach is a problem
#TODO do something better about handling duplicated SNP ID's in tensorqtl interaction test
#TODO submit bug report to tensorqtl to fix the chromosome indexing
#TODO find more stable solution than to manually overwrite cis.py (/project2/gilad/jpopp/ebQTL/.snakemake/conda/e138ad8ec0b845ac547107b3c3d4cf50/lib/python3.10/site-packages/tensorqtl/cis.py) with updates from https://github.com/broadinstitute/tensorqtl/commit/6879d887db13d880fc5fc0619906c4fbdbbb2843
#TODO in pseudobulk_assign_dynamic, get rid of the 25-sample cutoff for pseudotime bin inclusion, which is a holdout from cellid analysis

rule pseudobulk_assign_dynamic:
    resources:
        mem_mb=250000,
        time="01:30:00"
    input:
        raw="data/single_cell_objects/highpass/eb_raw.qc.h5ad",
        pseudotimed="data/trajectory_inference/{trajectory}_lineage/{trajectory}_lineage.{nbins}_pseudotime.adata"
    output:
        celltype_summary="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/samples_per_celltype.tsv",
        sample_summary="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/sample_summary.tsv",
        pseudobulk="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{trajectory}_{nbins}.pseudobulk_tmm.tsv"
    params:
        drops_allowed=2,
        min_cells_per_bin=5
    conda: "../slurmy/scvi.yml"
    script:
        "../code/dynamic_qtl_calling/pseudobulk_tmm-assign-dynamic.py"

rule pseudobulk_qc_dynamic:
    resources:
        mem_mb=250000,
        time="30:00"
    input:
        pseudobulk="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{trajectory}_{nbins}.pseudobulk_tmm.tsv",
        sample_summary="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/sample_summary.tsv"
    output:
        sample_summary_manual="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/sample_summary_manual.tsv"
    params:
        table_prefix = "data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/",
        fig_prefix = "figs/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/dynamic_qtl_calling/pseudobulk_tmm-qc-dynamic.R"

rule pseudobulk_agg_dynamic:
    resources:
        mem_mb=250000,
        time="30:00"
    input:
        pseudobulk="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{trajectory}_{nbins}.pseudobulk_tmm.tsv",
        sample_summary_manual="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/sample_summary_manual.tsv"
    output:
        all_expression="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/pseudobulk_all.tsv"
    params:
        table_prefix = "data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/",
        fig_prefix = "figs/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/dynamic_qtl_calling/pseudobulk_tmm-agg-dynamic.R"
        
rule plink_rewrite_keepers_dynamic:
    input:
        "data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/individuals.tsv"
    output:
        "data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/individuals_plink.tsv"
    shell:
        "awk -v OFS='\t' '{{ $2=$1; print}}' {input} > {output}"

rule plink_genotype_reformat_dynamic:
    resources:
        mem_mb=100000
    input:
        genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
        inds="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/individuals_plink.tsv"
    output:
        expand("data/dynamic_qtl_calling/{{trajectory}}_{{nbins}}/pseudobulk_tmm/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/genotypes_filtered_plink"
    shell:
        "code/dynamic_qtl_calling/plink_genotype_reformat.sh {input.genotypes} {input.inds} {params.prefix}"

rule plink_expression_reformat_dynamic:
    input:
        exp="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv"
    output:
        exp="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/expression.bed.gz"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/dynamic_qtl_calling/reformat_expression.R"

rule tensorqtl_interaction_prep:
    resources:
        mem_mb=250000,
        time="01:30:00"
    input:
        genotypes=expand("data/dynamic_qtl_calling/{{trajectory}}_{{nbins}}/pseudobulk_tmm/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/expression.bed.gz",
        cov="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/cell_line_pcs.tsv"
    output:
        covariate_df="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{n_cl_pcs}clpcs/covariate_df.tsv",
        genotype_df="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{n_cl_pcs}clpcs/genotype_df.tsv",
        variant_df="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{n_cl_pcs}clpcs/variant_df.tsv"
    params:
        plink_prefix="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/genotypes_filtered_plink",
        output_prefix="results/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{n_cl_pcs}clpcs/tensorqtl_interactions"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/dynamic_qtl_calling/tensorqtl_interaction_prep.py"

rule tensorqtl_nominal_interaction:
    resources:
        mem_mb=120000,
        partition="gpu2",
        gres="gpu:1",
        nodes=1
    input:
        exp="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/expression.bed.gz",
        covariate_df="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{n_cl_pcs}clpcs/covariate_df.tsv",
        genotype_df="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{n_cl_pcs}clpcs/genotype_df.tsv",
        variant_df="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{n_cl_pcs}clpcs/variant_df.tsv",
        pseudotime="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/pseudotime.tsv"
    output:
        expand("results/dynamic_qtl_calling/{{trajectory}}_{{nbins}}/pseudobulk_tmm/{{n_cl_pcs}}clpcs/tensorqtl_interactions.cis_qtl_pairs.chr{i}.parquet", i=range(1, 23))
    params:
        plink_prefix="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/genotypes_filtered_plink",
        output_prefix="results/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{n_cl_pcs}clpcs/tensorqtl_interactions"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/dynamic_qtl_calling/tensorqtl_interaction_nominal.py"

rule tensorqtl_merge_interaction:
    resources:
        mem_mb=200000,
        time="30:00"
    input:
        expand("results/dynamic_qtl_calling/{{trajectory}}_{{nbins}}/pseudobulk_tmm/{{n_cl_pcs}}clpcs/tensorqtl_interactions.cis_qtl_pairs.chr{i}.parquet", i=range(1, 23))
    output:
        merged_df="results/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{n_cl_pcs}clpcs/tensorqtl_interactions.all.tsv"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/dynamic_qtl_calling/tensorqtl_merge_interaction.py"
