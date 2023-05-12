#TODO download_fca
#TODO incorporate pull_differentiation_genesets.R, maybe wrangle_geneset.R

rule list_pcgs:
    resources:
        partition="broadwl",
        mem_mb=5000
    input:
        "/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:
        "data/fca/protein_coding_genes.tsv"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/list_pcgs.R"

rule subset_fca:
    resources:
        mem_mb=50000,
        partition="broadwl"
    input:
        counts="data/fca/counts.sampled.rds",
        cell_metadata="data/fca/cell_metadata.rds",
        gene_metadata="data/fca/gene_metadata.rds",
        pc_genes="data/fca/protein_coding_genes.tsv"
    output:
        fca_subsampled="data/fca/counts.subsampled.sce"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/subset_fca.R"

rule preprocess_fca:
    resources:
        mem_mb=50000,
        partition="broadwl"
    input:
        fca_subsampled="data/fca/counts.subsampled.sce"
    output:
        hvgs="data/fca/fca_subsampled_hvg.tsv",
        fca_subsampled_lognorm="data/fca/counts.subsampled.lognorm.sce"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/preprocess_fca.R"

rule subset_ebs:
    resources:
        mem_mb=450000,
        partition="bigmem2",
        time="02:00:00"
    input:
        "data/single_cell_objects/eb_pflog1ppfnorm.h5ad",
        "data/fca/fca_subsampled_hvg.tsv"
    output:
        "data/single_cell_objects/eb_pflog1ppfnorm.fca_hvg.h5ad"
    conda:
        "../slurmy/scvi.yml"
    script:
        "../code/annotation/eb_subset_to_fca_hvgs.py"

rule h5ad_to_sce:
    resources:
        mem_mb=250000,
        partition="bigmem2"
    input:
        h5ad="data/single_cell_objects/eb_pflog1ppfnorm.fca_hvg.h5ad"
    output:
        sce="data/single_cell_objects/eb_pflog1ppfnorm.fca_hvg.sce"
    params:
        X_name="logcounts"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/convert_h5ad_to_sce.R"
        
rule learn_fca_signatures:
    resources:
        mem_mb=200000,
        partition="bigmem2"
    input:
        fca="data/fca/counts.subsampled.lognorm.sce",
        hvg="data/fca/fca_subsampled_hvg.tsv"
    output:
        fca_embedded="data/fca/counts.subsampled.lognorm.mca.sce",
        signatures="data/fca/fca_signatures.rds"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/learn_fca_signatures.R"
    
rule classify_ebs:
    resources:
        mem_mb=250000,
        partition="bigmem2"
    input:
        eb_hvgs="data/single_cell_objects/eb_pflog1ppfnorm.fca_hvg.sce",
        signatures="data/fca/fca_signatures.rds"
    output:
        mca_embedding="data/fca/eb_mca_embedding.tsv",
        cellid_labels="data/fca/eb_cellid_labels.tsv",
        cellid_all_enrichments="data/fca/eb_cellid_enrichments.tsv"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/classify_ebs.R"

rule sce_to_h5ad:
    resources:
        mem_mb=250000,
        partition="bigmem2"
    input:
        sce="data/fca/counts.subsampled.lognorm.mca.sce"
    output:
        h5ad="data/fca/counts.subsampled.lognorm.mca.h5ad"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/convert_sce_to_h5ad.R"
        
rule scdrs_preprocessing_fca:
    resources:
        mem_mb=250000,
        partition="bigmem2"
    input:
        fca="data/fca/counts.subsampled.lognorm.mca.h5ad"
    output:
        preprocessed="data/fca/counts.subsampled.lognorm.mca.scdrs_preprocessed.h5ad"
    conda:
        "../slurmy/scvi.yml"
    script:
        "../code/annotation/scdrs_preprocessing_fca.py"
        
rule scdrs_preprocessing_eb:
    resources:
        mem_mb=250000
    input:
        eb="data/single_cell_objects/eb_pflog1ppfnorm.fca_hvg.h5ad",
        mca="data/fca/eb_mca_embedding.tsv",
        celltypes="data/fca/eb_cellid_labels.tsv"
    output:
        eb_preprocessed="data/single_cell_objects/eb_pflog1ppfnorm.fca_hvg.mca.subsampled.scdrs_preprocessed.h5ad"
    conda:
        "../slurmy/scvi.yml"
    script:
        "../code/annotation/scdrs_preprocessing_eb.py"
        
rule scdrs_diff_scores:
    resources:
        mem_mb=250000,
        time="05:00:00"
    input:
        diff_gs="data/gene_sets/c5.go.bp.Hs.symbols.differentiation.gs",
        fca="data/fca/counts.subsampled.lognorm.mca.scdrs_preprocessed.h5ad",
        eb="data/single_cell_objects/eb_pflog1ppfnorm.fca_hvg.mca.subsampled.scdrs_preprocessed.h5ad"
    params:
        table_prefix = "data/module_scoring/",
        fig_prefix = "figs/module_scoring/fca_eb_comp/"
    output:
        "figs/module_scoring/fca_eb_comp/GOBP_NEUROEPITHELIAL_CELL_DIFFERENTIATION.png"
    conda:
        "../slurmy/scvi.yml"
    script:
        "../code/annotation/score_differentiation_modules.py"
