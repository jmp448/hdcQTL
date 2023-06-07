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
    # NOTE - calling it lineage here just to exclude 'eb-' from the wildcard
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

rule infer_pseudotime:
    # NOTE - {nbins} wildcard here EXCLUDES 'trimmed', which is different from the dynamic_qtl_calling rules file
    # this allows this rule to be differentiated from the next one
    resources:
        mem_mb=250000,
        time="30:00"
    input:
        adata="data/trajectory_inference/{trajectory}_lineage/{trajectory}_lineage.adata",
        stage_labels="data/trajectory_inference/{trajectory}_lineage/{trajectory}_lineage_stages.tsv"
    output:
        adata_pseudotime="data/trajectory_inference/{trajectory}_lineage/{trajectory}_lineage.pseudotime.adata"
    conda: "../slurmy/scvi.yml"
    script:
        "../code/trajectory_inference/infer_pseudotime.py"
        
rule infer_pseudotime_trim_ends:
    # NOTE - {nbins} wildcard here EXCLUDES 'trimmed', which is different from the dynamic_qtl_calling rules file
    # this allows this rule to be differentiated from the previous one
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



