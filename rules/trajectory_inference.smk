#TODO analyze different weighting options

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
        gs="data/gene_sets/{trajectory}_lineage.tsv"
    output:
        scdrs="data/trajectory_inference/{trajectory}_lineage/scdrs_scores.tsv",
        adata_filtered="data/trajectory_inference/{trajectory}_lineage/eb_{trajectory}_lineage.adata"
    conda: "../slurmy/scvi.yml"
    script:
        "../code/trajectory_inference/trajectory_isolation.py"

rule infer_pseudotime:
    resources:
        mem_mb=250000,
        time="30:00"
    input:
        adata="data/trajectory_inference/{trajectory}_lineage/eb_{trajectory}_lineage.adata",
        stage_labels="data/trajectory_inference/{trajectory}_lineage/{trajectory}_lineage_stages.tsv"
    output:
        adata_pseudotime="data/trajectory_inference/{trajectory}_lineage/eb_{trajectory}_lineage.pseudotime.adata"
    conda: "../slurmy/scvi.yml"
    script:
        "../code/trajectory_inference/infer_pseudotime.py"


