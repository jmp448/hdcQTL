#TODO tidy up QC
#TODO look at batch effects
#TODO update the text file outputs to include cell names
#TODO go back to kenneth and katie's pipelines

rule create_anndata:
    resources:
        mem_mb=450000,
        time="03:00:00"
    input:
        "data/single_cell_objects/highpass/eb_raw.h5ad"
    output:
        "data/single_cell_objects/highpass/eb_raw.qc.h5ad"
    conda: "../slurmy/scvi.yml"
    script:
        "../code/sc_preprocessing/quality_control.py"

# rule quality_control:
#     resources:
#         mem_mb=450000,
#         time="03:00:00"
#     input:
#         "data/single_cell_objects/highpass/eb_raw.h5ad"
#     output:
#         "data/single_cell_objects/highpass/eb_raw.qc.h5ad"
#     conda: "../slurmy/scvi.yml"
#     script:
#         "../code/sc_preprocessing/quality_control.py"

rule normalize_counts:
    resources:
        mem_mb=450000,
        time="03:00:00"
    input:
        "data/single_cell_objects/highpass/eb_raw.h5ad"
    output:
        "data/single_cell_objects/eb_pflog1ppfnorm.h5ad",
        "data/single_cell_objects/eb_pflog1ppfnorm.hvg.h5ad"
    conda: "../slurmy/scvi.yml"
    script:
        "../code/sc_preprocessing/normalization.py"
        
rule scvi_embedding:
    resources:
        partition="gpu2",
        mem_mb=100000,
        time="03:00:00",
        gres="gpu:2",
        nodes=2
    input:
        "data/single_cell_objects/eb_pflog1ppfnorm.hvg.h5ad"
    output:
        "data/single_cell_objects/eb_pflog1ppfnorm.hvg.scvi_embedding.txt",
        "data/single_cell_objects/eb_pflog1ppfnorm.hvg.scvi_embedding.h5ad"
    params:
        "models/scvi_ldvae"
    conda: "../slurmy/scvi.yml"
    script:
        "../code/sc_preprocessing/scvi_embedding.py"
        
rule umap_embedding:
    resources:
        mem_mb=450000,
        time="04:00:00"
    input:
        "data/single_cell_objects/eb_pflog1ppfnorm.hvg.scvi_embedding.h5ad"
    output:
        "data/single_cell_objects/eb_pflog1ppfnorm.hvg.umap_embedding.txt",
        "data/single_cell_objects/eb_pflog1ppfnorm.hvg.umap_embedding.h5ad"
    conda: "../slurmy/scvi.yml"
    script:
        "../code/sc_preprocessing/umap_embedding.py"

