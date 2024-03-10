import glob

## DYNAMIC EQTL ANALYSIS
rule remove_gtex_overlap_dynamic:
    resources:
      mem_mb=100000
    input:
      dynamics_bed="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/original/dynamic-signif_variant_gene_pairs.bed",
      gtex_bed="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/original/dynamic-signif_variant_gene_pairs.all_tissue_overlap.bed"
    output:
      non_gtex_dynamics="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.snplist.txt",
      non_gtex_dynamics_bed="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.bed"
    conda:
        "../slurmy/r-sva.yml"
    script:
        "../code/complex_trait_analysis/remove_gtex_overlap_dynamic.R"
        
# to get the opentargets overlap, run opentargets_analysis/opentargets_query_dynamic.R
# #inputs
# dynamic_qtl_bed_loc <- "results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.bed"
# bim_loc <- "data/dynamic_qtl_calling/all_trajectories_combined/genotypes_filtered_plink.bim"
# 
# #outputs
# snplist_loc <- "results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.snplist.txt"
# dynamic_qtl_gwas_bed_loc <- "results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.bed"

rule list_yri_tagged_opentargets_eqtls:
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
	      genotypes=expand("data/genotypes/yri_all.hg38.deduplicated.{out}", out=['bed', 'bim', 'fam']),
	      snps="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.snplist.txt"
    output:
    	  expand("results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.snplist.{ext}", ext=['tags.list', 'tags', 'nosex', 'log'])
    params:
        genotypes="data/genotypes/yri_all.hg38.deduplicated",
        outfiles="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.snplist"
    shell:
	      """
	      module load plink
	      plink --bfile {params.genotypes} --show-tags {input.snps} --tag-r2 0.5 --tag-kb 1000 --list-all --out {params.outfiles}
	      """

rule expand_bed_to_tag_snps_opentargets:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        eqtls="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/original/dynamic-signif_variant_gene_pairs.bed",
        tags="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.snplist.tags.list",
        bim_file="data/dynamic_qtl_calling/all_trajectories_combined/genotypes_filtered_plink.bim"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/tagging-dynamic-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/complex_trait_analysis/expand_bed_to_tag_snps.R"
        
rule refine_opentargets_dynamic_qtls:
    input:
        non_gtex_gwas_eqtls="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.bed",
        gtex_overlap_tag_variants="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/tagging-dynamic-signif_variant_gene_pairs.all_tissue_overlap.bed"
    output:
        filtered_bed="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.tagged_gtex_removed.bed"
    conda:
        "../slurmy/r-sva.yml"
    script:
        "../code/complex_trait_analysis/refine_gwas_qtl_snplist.R"

### CRM OVERLAP
rule remove_gtex_overlap_crm:
    resources:
      mem_mb=100000
    input:
      crm_bed="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/crm-all-signif_variant_gene_pairs.bed",
      gtex_bed="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/crm-all-signif_variant_gene_pairs.all_tissue_overlap.bed"
    output:
      non_gtex_crm="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/crm-all-signif_variant_gene_pairs.gtex_removed.snplist.txt",
      non_gtex_crm_bed="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/crm-all-signif_variant_gene_pairs.gtex_removed.bed"
    conda:
        "../slurmy/r-sva.yml"
    script:
        "../code/complex_trait_analysis/remove_gtex_overlap_crm.R"
        
# # to get the opentargets overlap, run opentargets_analysis/opentargets_query_crm.R
# #inputs
# crm_qtl_bed_loc <- "results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/crm-all-signif_variant_gene_pairs.gtex_removed.bed"
# bim_loc <- "data/genotypes/yri_maf0.1_all.hg38.bim"
# 
# #outputs
# snplist_loc <- "results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/crm-all-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.snplist.txt"
# crm_qtl_gwas_bed_loc <- "results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/crm-all-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.bed"

rule list_yri_tagged_opentargets_eqtls_crm:
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
	      genotypes=expand("data/genotypes/yri_all.hg38.deduplicated.{out}", out=['bed', 'bim', 'fam']),
	      snps="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/crm-all-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.snplist.txt"
    output:
    	  expand("results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/crm-all-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.snplist.{ext}", ext=['tags.list', 'tags', 'nosex', 'log'])
    params:
        genotypes="data/genotypes/yri_all.hg38.deduplicated",
        outfiles="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/crm-all-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.snplist"
    shell:
	      """
	      module load plink
	      plink --bfile {params.genotypes} --show-tags {input.snps} --tag-r2 0.5 --tag-kb 1000 --list-all --out {params.outfiles}
	      """

rule expand_bed_to_tag_snps_opentargets_crm:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        eqtls="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/crm-all-signif_variant_gene_pairs.bed",
        tags="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/crm-all-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.snplist.tags.list",
        bim_file="data/genotypes/yri_maf0.1_all.hg38.bim"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/tagging-crm-all-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/complex_trait_analysis/expand_bed_to_tag_snps.R"

rule refine_opentargets_dynamic_qtls_crm:
    input:
        non_gtex_gwas_eqtls="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/crm-all-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.bed",
        gtex_overlap_tag_variants="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/tagging-crm-all-signif_variant_gene_pairs.all_tissue_overlap.bed"
    output:
        filtered_bed="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/crm-all-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.tagged_gtex_removed.bed"
    conda:
        "../slurmy/r-sva.yml"
    script:
        "../code/complex_trait_analysis/refine_gwas_qtl_snplist.R"


