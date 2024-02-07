
rule normalize_counts_wenhe:
    resources:
        mem_mb=450000,
        time="03:00:00"
    input:
        "data/single_cell_objects/highpass/eb_raw.h5ad"
    output:
        "data/fca_wenhe/eb_pflog1ppfnorm.h5ad",
        "data/fca_wenhe/eb_pflog1ppfnorm.hvg.h5ad"
    conda: "../slurmy/scvi.yml"
    script:
        "../code/annotation_wenhe/normalization.py"


rule preprocess_fca_wenhe:
    resources:
        mem_mb=50000,
        partition="broadwl"
    input:
        fca_subsampled="data/fca/counts.subsampled.sce"
    output:
        hvgs="data/fca_wenhe/fca_subsampled_hvg.tsv",
        fca_subsampled_lognorm="data/fca_wenhe/counts.subsampled.lognorm.sce"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation_wenhe/preprocess_fca.R"

rule subset_ebs_wenhe:
    resources:
        mem_mb=500000,
        partition="bigmem2",
        time="02:00:00"
    input:
        "data/fca_wenhe/eb_pflog1ppfnorm.h5ad",
        "data/fca_wenhe/fca_subsampled_hvg.tsv"
    output:
        "data/fca_wenhe/eb_pflog1ppfnorm.fca_hvg.h5ad"
    conda:
        "../slurmy/scvi.yml"
    script:
        "../code/annotation_wenhe/eb_subset_to_fca_hvgs.py"

rule h5ad_to_sce_wenhe:
    resources:
        mem_mb=250000,
        partition="bigmem2"
    input:
        h5ad="data/fca_wenhe/eb_pflog1ppfnorm.fca_hvg.h5ad"
    output:
        sce="data/fca_wenhe/eb_pflog1ppfnorm.fca_hvg.sce"
    params:
        X_name="logcounts"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation_wenhe/convert_h5ad_to_sce.R"

rule classify_ebs_wenhe:
    resources:
        mem_mb=250000,
        partition="bigmem2"
    input:
        eb_hvgs="data/fca_wenhe/eb_pflog1ppfnorm.fca_hvg.sce",
        signatures="/project2/gilad/wenhe/EB/batch3/data/fca/fca_signatures_8KHVG_100plus50DE_merged.rds"
    output:
        mca_embedding="data/fca_wenhe/eb_mca_embedding.tsv",
        cellid_labels="data/fca_wenhe/eb_cellid_labels.tsv",
        cellid_all_enrichments="data/fca_wenhe/eb_cellid_enrichments.tsv"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation_wenhe/classify_ebs.R"

