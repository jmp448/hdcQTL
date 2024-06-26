# PREPROCESSING & HELPER FUNCTIONS
rule process_gtf:
    # gtf from https://www.10xgenomics.com/support/software/cell-ranger/downloads
    resources:
        mem_mb=10000,
        time="30:00"
    input:
        gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf",
        tss_loc="data/gencode/gencode.hg38.filtered.tss.tsv",
        bed_loc="data/gencode/gencode.hg38.filtered.tss.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/benchmarking_static_qtl_calling/gene_locs.R"
  
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

rule h5ad_to_sce:
    resources:
        mem_mb=250000,
        partition="bigmem2"
    input:
        h5ad=ancient("data/single_cell_objects/eb_pflog1ppfnorm.fca_hvg.h5ad")
    output:
        sce="data/single_cell_objects/eb_pflog1ppfnorm.fca_hvg.sce"
    params:
        X_name="logcounts"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/convert_h5ad_to_sce.R"

# 77 CELL TYPE CLASSIFIER
## Classification based on full data
rule learn_signatures_77:
    resources:
        mem_mb=200000,
        time="06:00:00"
    input:
        counts="data/fca/counts.sampled.rds",
        cell_metadata="data/fca/cell_metadata.rds",
        gene_metadata="data/fca/gene_metadata.rds",
        pc_genes="data/fca/protein_coding_genes.tsv",
        hv_genes="data/fca/fca_subsampled_hvg.tsv"
    output:
        fca_77="data/fca/counts.subsampled.lognorm.mca.77celltypes.sce",
        signatures_77="data/fca/fca_signatures.77celltypes.rds"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/learn_signatures_77.R"

rule classify_ebs_77:
    resources:
        mem_mb=250000,
        time="06:00:00"
    input:
        eb_hvgs="data/single_cell_objects/eb_pflog1ppfnorm.fca_hvg.sce",
        signatures="data/fca/fca_signatures.77celltypes.rds"
    output:
        mca_embedding="data/fca/eb_mca_embedding.77celltypes.tsv",
        cellid_labels="data/fca/eb_cellid_labels.77celltypes.tsv",
        cellid_all_enrichments="data/fca/eb_cellid_enrichments.77celltypes.tsv"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/classify_ebs.R"

## 77 cell type classifier assessment based on 50/50 train/test split
rule assess_signatures_77:
    resources:
        mem_mb=200000,
        time="06:00:00"
    input:
        counts="data/fca/counts.sampled.rds",
        cell_metadata="data/fca/cell_metadata.rds",
        gene_metadata="data/fca/gene_metadata.rds",
        pc_genes="data/fca/protein_coding_genes.tsv"
    output:
        train_77="data/fca/fca_train.lognorm.mca.77celltypes.sce",
        test_77="data/fca/fca_test.lognorm.mca.77celltypes.sce",
        signatures_77="data/fca/assessment_signatures.77celltypes.rds"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/assess_signatures_77.R"
        
rule compare_signatures_77:
    resources:
        mem_mb=200000,
        time="06:00:00"
    input:
        test_77="data/fca/fca_test.lognorm.mca.77celltypes.sce",
        signatures_77="data/fca/assessment_signatures.77celltypes.rds"
    output:
        labeled_test_77="data/fca/fca_test.lognorm.labeled.77celltypes.sce"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/compare_signatures_77.R"

# 33 CELL TYPE CLASSIFIER [** THIS IS THE ONE USED FOR DOWNSTREAM ANALYSES **]
# Classifier development based on full dataset
rule subset_fca:
    resources:
        mem_mb=50000,
        partition="broadwl"
    input:
        counts="data/fca/counts.sampled.rds",
        cell_metadata="data/fca/cell_metadata.rds",
        gene_metadata="data/fca/gene_metadata.rds",
        pc_genes="data/gencode/gencode.hg38.filtered.gtf"
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

## Classifier assessment based on 50/50 train/test split
rule assess_signatures_subset:
    resources:
        mem_mb=200000,
        time="06:00:00"
    input:
        sce="data/fca/counts.subsampled.sce",
        pc_genes="data/fca/protein_coding_genes.tsv"
    output:
        train_subset="data/fca/fca_train.lognorm.mca.subset_celltypes.sce",
        test_subset="data/fca/fca_test.lognorm.mca.subset_celltypes.sce",
        signatures_subset="data/fca/assessment_signatures.subset_celltypes.rds"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/assess_signatures_subset.R"
        
rule compare_signatures_subset:
    resources:
        mem_mb=200000,
        time="06:00:00"
    input:
        test_subset="data/fca/fca_test.lognorm.mca.subset_celltypes.sce",
        signatures_subset="data/fca/assessment_signatures.subset_celltypes.rds"
    output:
        labeled_test_subset="data/fca/fca_test.lognorm.labeled.subset_celltypes.sce"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/compare_signatures_subset.R"

rule save_signatures_to_tsv:
    resources:
        mem_mb=10000
    input:
        fca_signatures="data/fca/fca_signatures.rds"
    output:
        fca_signatures_tsv="data/fca/fca_signatures.tsv"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/wrangle_gene_sets.R"
        
# 34 CELL TYPE CLASSIFIER (33 FCA CELL TYPES + IPSC SIGNATURE FROM PLURINET)
rule merge_fca_ipsc_signatures:
    resources:
        mem_mb=10000
    input:
        plurinet_signatures="data/gene_sets/MUELLER_PLURINET.v2023.2.Hs.gmt",
        fca_signatures="data/fca/fca_signatures.rds"
    output:
        all_signatures="data/fca/fca_signatures.with_ipsc.rds"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/learn_signatures_with_ipsc.R"

rule classify_ebs_with_ipsc:
    resources:
        mem_mb=250000,
        time="06:00:00"
    input:
        eb_hvgs="data/single_cell_objects/eb_pflog1ppfnorm.fca_hvg.sce",
        signatures="data/fca/fca_signatures.with_ipsc.rds" 
    output:
        mca_embedding="data/fca/eb_mca_embedding.with_ipsc.tsv",
        cellid_labels="data/fca/eb_cellid_labels.with_ipsc.tsv",
        cellid_all_enrichments="data/fca/eb_cellid_enrichments.with_ipsc.tsv"
    conda:
        "../slurmy/r-fca.yml"
    script:
        "../code/annotation/classify_ebs.R"
