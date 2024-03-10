# manually overwrite cis.py (/project2/gilad/jpopp/ebQTL/.snakemake/conda/e138ad8ec0b845ac547107b3c3d4cf50/lib/python3.10/site-packages/tensorqtl/cis.py) with updates from https://github.com/broadinstitute/tensorqtl/commit/6879d887db13d880fc5fc0619906c4fbdbbb2843
# similar issue, /project2/gilad/jpopp/ebQTL/.snakemake/conda/633df945cb6942e683de1035f6270d68_/lib/python3.11/site-packages/tensorqtl/eigenmt.py had to change np.float (deprecated) to np.float64

# Generate pseudobulk expression data and latent covariates for QTL calling
rule pseudobulk_assign_dynamic:
    resources:
        mem_mb=250000,
        time="01:30:00"
    input:
        raw="data/single_cell_objects/highpass/eb_raw.h5ad",
        pseudotimed="data/trajectory_inference/{trajectory}_lineage/{trajectory}_lineage.{nbins}_pseudotime.adata",
        metadata="/project2/gilad/katie/ebQTL/CombinedFormationAndCollectionMetadata_102andPilot_SWAPSANDCONTAMINATIONADDED_012522.csv"
    output:
        donor_summary="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/samples_per_donor.tsv",
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
        mem_mb=500000,
        time="05:00:00"
    input:
        pseudobulk="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{trajectory}_{nbins}.pseudobulk_tmm.tsv",
        sample_summary="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/sample_summary.tsv"
    output:
        sample_summary_manual="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/sample_summary_manual.tsv"
    params:
        table_prefix = "data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/",
        fig_prefix = "figs/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/dynamic_qtl_calling/pseudobulk_tmm-qc-dynamic.R"

rule pseudobulk_agg_dynamic:
    resources:
        mem_mb=500000,
        time="05:00:00"
    input:
        pseudobulk="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{trajectory}_{nbins}.pseudobulk_tmm.tsv",
        metadata="/project2/gilad/katie/ebQTL/CombinedFormationAndCollectionMetadata_102andPilot_SWAPSANDCONTAMINATIONADDED_012522.csv",
        sample_summary_manual="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/sample_summary_manual.tsv"
    output:
        raw_expression="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/pseudobulk_raw.tsv",
        norm_expression="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/expression.tsv",
        covariates="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/covariates.tsv",
        individuals="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/individuals.tsv",
        pseudotime="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/pseudotime.tsv"
    params:
        table_prefix = "data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/",
        fig_prefix = "figs/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/dynamic_qtl_calling/pseudobulk_tmm-agg-dynamic.R"
        
# Data wrangling to prep for tensorqtl
rule plink_rewrite_keepers_dynamic:
    input:
        "data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/individuals.tsv"
    output:
        "data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/individuals_plink.tsv"
    shell:
        "awk -v OFS='\t' '{{ $2=$1; print}}' {input} > {output}"

rule plink_genotype_reformat_dynamic:
    resources:
        mem_mb=100000
    input:
        genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
        inds="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/individuals_plink.tsv"
    output:
        expand("data/dynamic_qtl_calling/{{trajectory}}_{{nbins}}/pseudobulk_tmm/{{pca}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/genotypes_filtered_plink"
    shell:
        "code/dynamic_qtl_calling/plink_genotype_reformat.sh {input.genotypes} {input.inds} {params.prefix}"

rule plink_expression_reformat_dynamic:
    input:
        exp="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv"
    output:
        exp="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/expression.bed.gz"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/dynamic_qtl_calling/reformat_expression.R"

rule tensorqtl_interaction_prep:
    resources:
        mem_mb=250000,
        time="02:30:00"
    input:
        genotypes=expand("data/dynamic_qtl_calling/{{trajectory}}_{{nbins}}/pseudobulk_tmm/{{pca}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/expression.bed.gz",
        cov="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/covariates.tsv",
        pseudotime="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/pseudotime.tsv"
    output:
        covariate_df="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/{n_cl_pcs}clpcs/covariate_df.tsv",
        genotype_df="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/{n_cl_pcs}clpcs/genotype_df.tsv",
        variant_df="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/{n_cl_pcs}clpcs/variant_df.tsv"
    params:
        plink_prefix="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/genotypes_filtered_plink",
        output_prefix="results/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/{n_cl_pcs}clpcs/tensorqtl_interactions"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/dynamic_qtl_calling/tensorqtl_interaction_prep.py"

# Run tensorqtl
rule tensorqtl_interaction:
    resources:
        mem_mb=80000,
        partition="gpu2",
        gres="gpu:1",
        nodes=1
    input:
        exp="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/expression.bed.gz",
        covariate_df="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/{n_cl_pcs}clpcs/covariate_df.tsv",
        genotype_df="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/{n_cl_pcs}clpcs/genotype_df.tsv",
        variant_df="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/{n_cl_pcs}clpcs/variant_df.tsv",
        pseudotime="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/pseudotime.tsv"
    output:
        expand("results/dynamic_qtl_calling/{{trajectory}}_{{nbins}}/pseudobulk_tmm/{{pca}}/{{n_cl_pcs}}clpcs/tensorqtl_interactions.cis_qtl_pairs.chr{i}.parquet", i=range(1, 23)),
        "results/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/{n_cl_pcs}clpcs/tensorqtl_interactions.cis_qtl_top_assoc.txt.gz"
    params:
        plink_prefix="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/genotypes_filtered_plink",
        output_prefix="results/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/{n_cl_pcs}clpcs/tensorqtl_interactions"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/dynamic_qtl_calling/tensorqtl_interaction.py"

rule tensorqtl_merge_interaction:
    resources:
        mem_mb=200000,
        time="30:00"
    input:
        expand("results/dynamic_qtl_calling/{{trajectory}}_{{nbins}}/pseudobulk_tmm/{{pca}}/{{n_cl_pcs}}clpcs/tensorqtl_interactions.cis_qtl_pairs.chr{i}.parquet", i=range(1, 23))
    output:
        merged_df="results/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/{n_cl_pcs}clpcs/tensorqtl_interactions.all.tsv"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/dynamic_qtl_calling/tensorqtl_merge_interaction.py"

rule list_all_dynamic_eqtls:
    resources:
        mem_mb=200000,
        time="30:00"
    input:
        egenes="results/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/{n_cl_pcs}clpcs/tensorqtl_interactions.cis_qtl_top_assoc.txt.gz",
        alltests="results/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/{n_cl_pcs}clpcs/tensorqtl_interactions.all.tsv"
    output:
        sighits="results/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/{n_cl_pcs}clpcs/tensorqtl_interactions.cis_qtl_all_signif.fdr{fdr}.tsv"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/dynamic_qtl_calling/list_all_dynamic_eqtls.R"
        
# Classify effects into early/late/switch and write to BED
rule classify_dynamic_trajhits:
    resources:
        mem_mb=65000,
        time="1:00:00"
    input:
        genotypes="data/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/genotype_df.tsv",
        qtls="results/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.cis_qtl_all_signif.fdr0.1.tsv",
        pseudotime="data/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/pseudotime.tsv",
        bim_file="data/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim",
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
    output:
        early="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/earlydynamic-{trajectory}-signif_variant_gene_pairs.bed",
        late="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/latedynamic-{trajectory}-signif_variant_gene_pairs.bed",
        switch="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/switchdynamic-{trajectory}-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/dynamic_qtl_calling/classify_dynamic_eqtls.R"

# Write to BED file
rule merge_bims:
    input:
        cm_bim="data/dynamic_qtl_calling/eb-cm_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim",
        hep_bim="data/dynamic_qtl_calling/eb-hep_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim",
        neur_bim="data/dynamic_qtl_calling/eb-neur_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim"
    output:
        "data/dynamic_qtl_calling/all_trajectories_combined/genotypes_filtered_plink.bim"
    conda: "../slurmy/r-mashr.yml"
    shell:
        """
        cat {input} | sort -u -k 1,1 -k 4,4n > {output}
        """
        
# rule dynamic_to_bed_alltests:
#     resources:
#         mem_mb=50000,
#         time="15:00"
#     input:
#         neur_tests="results/dynamic_qtl_calling/eb-neur_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.all.tsv",
#         cm_tests="results/dynamic_qtl_calling/eb-cm_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.all.tsv",
#         hep_tests="results/dynamic_qtl_calling/eb-hep_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.all.tsv",
#         bim_file="data/dynamic_qtl_calling/all_trajectories_combined/genotypes_filtered_plink.bim",
#         gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
#     output:
#         bedfile="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/dynamic-all_variant_gene_pairs.bed"
#     conda: "../slurmy/r-mashr.yml"
#     script:
#         "../code/static_eqtl_followup/dynamic_to_bed_alltests.R"
        
rule dynamic_to_bed_allhits:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        cm_eqtls="results/dynamic_qtl_calling/eb-cm_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.cis_qtl_all_signif.fdr0.1.tsv",
        neur_eqtls="results/dynamic_qtl_calling/eb-neur_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.cis_qtl_all_signif.fdr0.1.tsv",
        hep_eqtls="results/dynamic_qtl_calling/eb-hep_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.cis_qtl_all_signif.fdr0.1.tsv",
        bim_file="data/dynamic_qtl_calling/all_trajectories_combined/genotypes_filtered_plink.bim",
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/original/dynamic-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/dynamic_to_bed_allhits.R"

rule dynamic_to_bed_trajhits:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        eqtls="results/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.cis_qtl_all_signif.fdr0.1.tsv",
        bim_file="data/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim",
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/original/dynamic-{trajectory}-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/dynamic_to_bed_hits.R"
        
