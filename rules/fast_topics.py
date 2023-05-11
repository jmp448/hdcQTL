rule subsample_for_pseudocells:
    resources:
        mem_mb=200000,
        time="30:00"
    input:
        adata="data/single_cell_objects/eb_pflog1ppfnorm.hvg.scvi_embedding.h5ad"
    output:
        adata="data/single_cell_objects/eb_pflog1ppfnorm.hvg.scvi_embedding.subset_for_pseudocells.h5ad"
    conda: 
      "../slurmy/scvi.yml"
    script:
        "../code/fast_topics/subsample_for_pseudocells.py"

rule pseudocell_mapping:
    resources:
        mem_mb=400000,
        time="01:00:00"
    input:
        adata_subsampled="data/single_cell_objects/eb_pflog1ppfnorm.hvg.scvi_embedding.subset_for_pseudocells.h5ad"
    output:
        pseudocell_mapping="data/fast_topics/cell_pseudocell_mapping.tsv"
    conda: 
      "../slurmy/scvi.yml"
    script:
        "../code/fast_topics/pseudocell_mapping.py"

rule pseudocell_aggregation:
    resources:
        mem_mb=450000,
        time="05:00:00"
    input:
        adata_raw="data/single_cell_objects/highpass/eb_raw.h5ad",
        pseudocell_mapping="data/fast_topics/cell_pseudocell_mapping.tsv"
    output:
        adata_pseudocell="data/single_cell_objects/eb_pseudocells_raw.h5ad"
    conda: 
      "../slurmy/scvi.yml"
    script:
        "../code/fast_topics/pseudocell_aggregation.py"

rule h5ad2seurat:
    resources:
        mem_mb = 50000,
        time = "00:30:00"
    input:
        anndata="data/single_cell_objects/eb_pseudocells_raw.h5ad"
    output:
        seurat="data/single_cell_objects/eb_pseudocells_raw.seurat.rds"
    conda:
        "../slurmy/r-fasttopics.yml"
    envmodules:
        ["hdf5_hl"]
    script:
        "../code/fast_topics/h5ad2seurat.R"
        
rule prepare_data_for_fasttopics:
    resources:
        partition = "gilad",
        mem_mb = 30000,
        time = "00:30:00"
    input:
        input_file="/project2/gilad/mli/eb_project/snakemake/0.test_dataset/test_data.rds"
    output:
        output_file="/project2/gilad/mli/eb_project/snakemake/0.test_dataset/test_data.RData"
    conda:
        "/project2/gilad/mli/eb_project/snakemake/slurmy/r-package.yml"
    script:
        "/project2/gilad/mli/eb_project/snakemake/3.prepare_data_for_fasttopics/prepare_data_for_fasttopics.R"
  
