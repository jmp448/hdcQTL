import pandas as pd

def list_cellregmap_output_files(wildcards):
    test_genes = set(pd.read_csv(f"results/static_qtl_calling/{wildcards.annotation}/pseudobulk_tmm/basic/8pcs/{wildcards.variant_group}_variant_gene_pairs.bed", sep="\t")['EB_HGNC'])
    return [f"results/cellregmap_eqtl_calling/{wildcards.annotation}/pseudobulk_tmm/basic/{wildcards.variant_group}.fasttopics_{wildcards.k}topics.{g}.cellregmap.tsv" for g in test_genes]

def list_cellregmap_beta_files(wildcards):
    iqtl_genes = set(pd.read_csv(f"results/cellregmap_eqtl_calling/{wildcards.annotation}/pseudobulk_tmm/basic/all_genes_merged.{wildcards.variant_group}.fasttopics_{wildcards.k}topics.cellregmap.sighits.tsv", sep="\t")['EB_HGNC'])
    return [f"results/cellregmap_eqtl_calling/{wildcards.annotation}/pseudobulk_tmm/basic/{wildcards.variant_group}.fasttopics_{wildcards.k}topics.{g}.cellregmap.betas.tsv" for g in iqtl_genes]

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
        exp="data/single_cell_objects/eb_pseudocells_normalized.nc",
        tabular_exp="data/single_cell_objects/eb_pseudocells_normalized.tsv"
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
        test_eqtl_file="results/static_eqtl_followup/qtl_sets/mash/original/mash-signif_variant_gene_pairs.bed",
        genotype_file="data/genotypes/yri_maf0.1_all.hg38.bed"
    output:
        filtered_eqtl_file="data/cellregmap/mash-signif_variant_gene_pairs.maf0.1.bed"
    conda: 
        "../slurmy/cellregmap.yml"
    script:
        "../code/cellregmap_eqtl_calling/filter_tests_fullpanel_maf_0.1.py"

rule run_interaction_test_fasttopics:
    resources:
        partition = "gilad",
        mem_mb = 20000,
        time = "24:00:00"
    input:
        test_eqtl_file="data/cellregmap/{variant_group}_variant_gene_pairs.maf0.1.bed",
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
        time = "10:00:00"
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

rule cellregmap_mtc:
    resources:
        mem_mb = 10000,
        time = "10:00"
    input:
        eqtls="results/cellregmap_eqtl_calling/{annotation}/pseudobulk_tmm/basic/all_genes_merged.{variant_group}.fasttopics_{k}topics.cellregmap.tsv"
    output:
        all_tests="results/cellregmap_eqtl_calling/{annotation}/pseudobulk_tmm/basic/all_genes_merged.{variant_group}.fasttopics_{k}topics.cellregmap.mtc.tsv",
        top_tests="results/cellregmap_eqtl_calling/{annotation}/pseudobulk_tmm/basic/all_genes_merged.{variant_group}.fasttopics_{k}topics.cellregmap.tophits.tsv",
        sig_hits="results/cellregmap_eqtl_calling/{annotation}/pseudobulk_tmm/basic/all_genes_merged.{variant_group}.fasttopics_{k}topics.cellregmap.sighits.tsv"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/cellregmap_eqtl_calling/cellregmap_mtc.R"

rule estimate_effect_sizes:
    resources:
        partition = "gilad",
        mem_mb = 20000,
        time = "24:00:00"
    input:
        interaction_eqtl_file="results/cellregmap_eqtl_calling/{annotation}/pseudobulk_tmm/basic/all_genes_merged.{variant_group}.fasttopics_{k}topics.cellregmap.sighits.tsv",
        sample_mapping_file = "data/cellregmap/pseudocell_metadata.tsv",
        genotype_file="data/genotypes/yri_maf0.1_all.hg38.bed" ,
        kinship_file = "data/genotypes/yri_kinship.tsv",
        exp = "data/single_cell_objects/eb_pseudocells_normalized.nc", 
        cell_contexts = "results/fast_topics/fasttopics_{k}topics_loadings.tsv"
    output:
        out="results/cellregmap_eqtl_calling/{annotation}/pseudobulk_tmm/basic/{variant_group}.fasttopics_{k}topics.{g}.cellregmap.betas.tsv"
    conda:
        "../slurmy/cellregmap.yml"
    script:
        "../code/cellregmap_eqtl_calling/single_gene_effect_estimates.py"

rule merge_effect_size_estimates:
    resources:
        partition = "gilad",
        mem_mb = 10000,
        time = "10:00:00"
    input:
        unpack(list_cellregmap_beta_files)
    output:
        "temp/cellregmap_eqtl_calling.{annotation}.pseudobulk_tmm.basic.all_genes_merged.{variant_group}.fasttopics_{k}topics.cellregmap.DUMMYDONEFILE.tsv"
    shell:
        """
        echo done > {output}
        """

rule crm_to_bed:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        crm_hits="results/cellregmap_eqtl_calling/eb_cellid/pseudobulk_tmm/basic/all_genes_merged.mash-signif.fasttopics_10topics.cellregmap.sighits.tsv",
        bim_file="data/genotypes/yri_maf0.1_all.hg38.bim",
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/crm-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/crm_to_bed.R"
        

rule cellregmap_mtc_all_signif:
    resources:
        mem_mb = 10000,
        time = "10:00"
    input:
        eqtls="results/cellregmap_eqtl_calling/{annotation}/pseudobulk_tmm/basic/all_genes_merged.{variant_group}.fasttopics_{k}topics.cellregmap.tsv"
    output:
        all_qtls="results/cellregmap_eqtl_calling/{annotation}/pseudobulk_tmm/basic/all_genes_merged.{variant_group}.fasttopics_{k}topics.cellregmap.all_signif.fdr0.1.tsv"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/cellregmap_eqtl_calling/cellregmap_mtc_all.R"

rule crm_full_to_bed:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        crm_hits="results/cellregmap_eqtl_calling/eb_cellid/pseudobulk_tmm/basic/all_genes_merged.mash-signif.fasttopics_10topics.cellregmap.all_signif.fdr0.1.tsv",
        bim_file="data/genotypes/yri_maf0.1_all.hg38.bim",
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/crm-all-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/cellregmap_eqtl_calling/crm_all_to_bed.R"


# ## TROAP
# rule filter_tests_troap:
#     resources:
#         mem_mb=50000,
#         time="3:00:00"
#     input:
#         test_eqtl_file="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/troap-region_variant_gene_pairs.bed",
#         genotype_file="data/genotypes/yri_maf0.1_all.hg38.bed"
#     output:
#         filtered_eqtl_file="data/cellregmap/troap-region_variant_gene_pairs.maf0.1.bed"
#     conda: 
#         "../slurmy/cellregmap.yml"
#     script:
#         "../code/cellregmap_eqtl_calling/filter_tests_fullpanel_maf_0.1.py"

