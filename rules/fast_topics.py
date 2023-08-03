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
        adata_pseudocell="data/single_cell_objects/eb_pseudocells_raw.h5ad",
        adata_pseudocell_names = "data/fast_topics/eb_pseudocells_raw.pseudocell_names.tsv",
        adata_gene_names = "data/fast_topics/eb_pseudocells_raw.gene_names.tsv"
    conda: 
      "../slurmy/scvi.yml"
    script:
        "../code/fast_topics/pseudocell_aggregation.py"

rule fasttopics_h5ad_to_sce:
    resources:
        mem_mb=250000,
        partition="bigmem2"
    input:
        h5ad="data/single_cell_objects/eb_pseudocells_raw.h5ad"
    output:
        sce="data/single_cell_objects/eb_pseudocells_raw.sce"
    params:
        X_name="counts"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/convert_h5ad_to_sce.R"

rule fast_topics:
    resources:
        mem_mb=100000,
        partition="gilad",
        time="48:00:00",
        ntasks_per_node=28
    input:
        pseudocells_sce="data/single_cell_objects/eb_pseudocells_raw.sce",
        adata_pseudocell_names = "data/fast_topics/eb_pseudocells_raw.pseudocell_names.tsv",
        adata_gene_names = "data/fast_topics/eb_pseudocells_raw.gene_names.tsv"
    output:
        fasttopics_fit="results/fast_topics/fasttopics_{k}topics_fit.rds"
    conda:
        "../slurmy/r-fasttopics.yml"
    script:
        "../code/fast_topics/fast_topics.R"

rule topic_de_analysis:
    resources:
        mem_mb=100000,
        partition="gilad",
        time="48:00:00",
        ntasks_per_node=28
    input:
        pseudocells_sce="data/single_cell_objects/eb_pseudocells_raw.sce",
        fasttopics_fit="results/fast_topics/fasttopics_{k}topics_fit.rds"
    output:
        de_analysis="results/fast_topics/fasttopics_{k}topics_de_analysis.Rdata"
    conda:
        "../slurmy/r-fasttopics.yml"
    script:
        "../code/fast_topics/topic_de_analysis.R"

rule get_context_loading_file:
    resources:
        mem_mb=10000,
        partition="gilad"
    input:
        de_result="results/fast_topics/fasttopics_{k}topics_de_analysis.Rdata"
    output:
        context_loading="results/fast_topics/fasttopics_{k}topics_loadings.tsv"
    conda:
        "../slurmy/r-fasttopics.yml"
    script:
        "../code/fast_topics/get_loadings.R"
