import pandas as pd

def list_cellregmap_output_files(wildcards):
    test_genes = set(pd.read_csv(f"results/static_qtl_calling/{wildcards.annotation}/pseudobulk_tmm/basic/8pcs/{wildcards.variant_group}_variant_gene_pairs.bed", sep="\t")['EB_HGNC'])
    return [f"results/cellregmap_eqtl_calling/{wildcards.annotation}/pseudobulk_tmm/basic/{wildcards.variant_group}.fasttopics_{wildcards.k}topics.{g}.cellregmap.tsv" for g in test_genes]

rule make_genotype_bed:
    resources:
        mem_mb=100000
    input:
        genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
        inds="data/genotypes/all_individuals_53.tsv"
    output:
        expand("data/genotypes/yri_maf0.1_all.hg38.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/genotypes/yri_maf0.1_all.hg38"
    shell:
        "code/cellregmap_eqtl_calling/make_genotype_bed.sh {input.genotypes} {input.inds} {params.prefix}"

rule make_kinship:
    resources:
        mem_mb=30000,
        time="1:00:00"
    input:
        genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
        inds="data/genotypes/all_individuals_53.tsv"
    output:
        expand("data/genotypes/yri_kinship.{ext}", ext=['rel', 'rel.id', 'log', 'nosex'])#,
        #pruned_bed=expand("data/genotypes/pruned_bed.{ext}", ext=['bed', 'log'])
    params:
        kinship_prefix="data/genotypes/yri_kinship",
        bed_prefix="data/genotypes/pruned_bed"
    shell:
        "code/cellregmap_eqtl_calling/make_kinship.sh {input.genotypes} {input.inds} {params.kinship_prefix} {params.bed_prefix}"

rule wrangle_kinship:
    resources:
        mem_mb=10000,
        time="1:00:00"
    input:
        kinship="data/genotypes/yri_kinship.rel",
        kinship_id="data/genotypes/yri_kinship.rel.id"
    output:
        kinship_tsv="data/genotypes/yri_kinship.tsv"
    conda: 
        "../slurmy/cellregmap.yml"
    script:
        "../code/cellregmap_eqtl_calling/wrangle_kinship.py"

rule compress_expression:
    resources:
        mem_mb=30000,
        time="1:00:00"
    input:
        pseudocell_adata="data/single_cell_objects/eb_pseudocells_raw.h5ad"
    output:
        exp="data/single_cell_objects/eb_pseudocells_normalized.nc"
    conda: 
        "../slurmy/cellregmap.yml"
    script:
        "../code/cellregmap_eqtl_calling/compress_expression.py"

rule make_pseudocell_metadata:
    resources:
        mem_mb=30000,
        time="1:00:00"
    input:
        pseudocell_map="data/fast_topics/cell_pseudocell_mapping.tsv",
        metadata="/project2/gilad/katie/ebQTL/CombinedFormationAndCollectionMetadata_102andPilot_SWAPSANDCONTAMINATIONADDED_012522.csv"
    output:
        pseudocell_metadata="data/cellregmap/pseudocell_metadata.tsv"
    conda: 
        "../slurmy/cellregmap.yml"
    script:
        "../code/cellregmap_eqtl_calling/make_pseudocell_metadata.py"

rule filter_tests:
    resources:
        mem_mb=50000,
        time="3:00:00"
    input:
        test_eqtl_file="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/8pcs/{variant_group}_variant_gene_pairs.bed",
        genotype_file="data/genotypes/yri_maf0.1_all.hg38.bed"
    output:
        filtered_eqtl_file="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/8pcs/{variant_group}_variant_gene_pairs.maf0.1.bed"
    conda: 
        "../slurmy/cellregmap.yml"
    script:
        "../code/cellregmap_eqtl_calling/filter_tests_fullpanel_maf_0.1.py"

rule run_interaction_test_fasttopics:
    resources:
        partition = "gilad",
        mem_mb = 20000,
        time = "1:00:00"
    input:
        test_eqtl_file="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/8pcs/{variant_group}_variant_gene_pairs.maf0.1.bed",
        sample_mapping_file = "data/cellregmap/pseudocell_metadata.tsv",
        genotype_file="data/genotypes/yri_maf0.1_all.hg38.bed" ,
        kinship_file = "data/genotypes/yri_kinship.tsv",
        exp = "data/single_cell_objects/eb_pseudocells_normalized.nc", 
        cell_contexts = "results/fast_topics/fasttopics_{k}topics_loadings.tsv"
    output:
        out="results/cellregmap_eqtl_calling/{annotation}/pseudobulk_tmm/basic/{variant_group}.fasttopics_{k}topics.{g}.cellregmap.tsv"
    conda:
        "../slurmy/cellregmap.yml"
    script:
        "../code/cellregmap_eqtl_calling/single_gene_interaction_test.py"
        
rule merge_interaction_test_fasttopics:
    resources:
        partition = "gilad",
        mem_mb = 10000,
        time = "01:00:00"
    input:
        unpack(list_cellregmap_output_files)
    output:
        "results/cellregmap_eqtl_calling/{annotation}/pseudobulk_tmm/basic/all_genes_merged.{variant_group}.fasttopics_{k}topics.cellregmap.tsv"
    params:
        tempfile="temp/cellregmap_tempfile.txt"
    shell:
        """
        # Combine input files into a single file
        cat {input} > {params.tempfile}
        
        # Filter lines to unique rows and save to output
        echo -e "EB_HGNC\tEB_VARIANT_ID\tP_CELLREGMAP" > {output}
        sort -u {params.tempfile} >> {output}
        
        # Remove temp file
        rm {params.tempfile}
        """
  
