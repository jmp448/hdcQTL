rule scdrs_preprocessing:
    resources:
        mem_mb=500000,
        time="01:30:00"
    input:
        adata="data/single_cell_objects/eb_pflog1ppfnorm.h5ad"
    output:
        adata_preprocessed="data/trajectory_inference/eb_scdrs_preprocessed.adata"
    conda: "../slurmy/scvi.yml"
    script:
        "../code/trajectory_inference/scdrs_preprocessing.py"

rule trajectory_isolation:
    resources:
        mem_mb=500000,
        time="01:30:00"
    input:
        adata="data/trajectory_inference/eb_scdrs_preprocessed.adata",
        gs="data/gene_sets/{lineage}_lineage.tsv"
    output:
        scdrs="data/trajectory_inference/eb-{lineage}_lineage/scdrs_scores.tsv",
        adata_filtered="data/trajectory_inference/eb-{lineage}_lineage/eb-{lineage}_lineage.adata"
    conda: "../slurmy/scvi.yml"
    script:
        "../code/trajectory_inference/trajectory_isolation.py"

rule infer_pseudotime_trim_ends:
    resources:
        mem_mb=250000,
        time="30:00"
    input:
        adata="data/trajectory_inference/{trajectory}_lineage/{trajectory}_lineage.adata"
    output:
        adata_pseudotime="data/trajectory_inference/{trajectory}_lineage/{trajectory}_lineage.{nbins}trimmed_pseudotime.adata"
    conda: "../slurmy/scvi.yml"
    script:
        "../code/trajectory_inference/infer_pseudotime_trim_ends.py"


