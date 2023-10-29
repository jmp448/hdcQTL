import glob

gtex_tissues = set(glob_wildcards("/project2/gilad/jpopp/GTEx_Analysis_v8_eQTL/{tissues}.v8.signif_variant_gene_pairs.txt").tissues)

rule reformat_GTEx_eqtls:
# convert GTEx significant eGene-eVariant pairs to BED format
    resources:
        mem_mb=150000,
        disk_mb=20000
    input:
        "/project2/gilad/jpopp/GTEx_Analysis_v8_eQTL/{tissue}.v8.signif_variant_gene_pairs.txt"
    output:
        "results/static_eqtl_followup/gtex/{tissue}.signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/reformat_gtex_snps.R"
        
rule sort_bed:
    input:
        "{fname}.bed"
    output:
        "{fname}.sorted.bed"
    shell:
        """
        module load bedtools
        sort -k1,1 -k2,2n {input} > {output}
        """

# note that I'm just copying some snpsets like the dynamic ones into this folder        
rule list_overlap_snps_sigtests_per_tissue:
    """
    This command will:
    1. find the SNP intersection of two bedfiles,
    2. filter to SNPs which were tested against the same gene
    """
    resources:
        mem_mb=10000
    input:
        #ebqtl="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/{snpset}_variant_gene_pairs.sorted.bed",
        ebqtl="results/static_eqtl_followup/qtl_sets/{analysis}/{snpset}_variant_gene_pairs.sorted.bed",
        gtex="results/static_eqtl_followup/gtex/{tissue}.signif_variant_gene_pairs.sorted.bed"
    output:
        "results/static_eqtl_followup/qtl_sets/{analysis}/{snpset}_variant_gene_pairs.{tissue}.overlap.bed"
    shell:
        """
        module load bedtools
        bedtools intersect -a {input.ebqtl} -b {input.gtex} -wa -wb | awk '$4 == $11' > {output}
        """

rule merge_overlap_sigtests:
    resources:
        mem_mb=200000,
        time="30:00"
    input:
        expand("results/static_eqtl_followup/qtl_sets/{{analysis}}/{{snpset}}_variant_gene_pairs.{tissue}.overlap.bed", tissue=gtex_tissues)
    output:
        "results/static_eqtl_followup/qtl_sets/{analysis}/{snpset}_variant_gene_pairs.all_tissue_overlap.bed"
    params:
        tempfile="temp/gtex_overlap_tempfile.txt"
    shell:
        """
        # Combine input files into a single file
        cat {input} > {params.tempfile}
        
        # Filter lines to unique rows
        sort -u {params.tempfile} > {output}
        
        # Remove temporary file
        rm {params.tempfile}
        """

# rule list_window_snps_sigtests_per_tissue:
#     """
#     This command will:
#     1. find the SNP intersection of two bedfiles,
#     2. filter to SNPs which were tested against the same gene
#     """
#     resources:
#         mem_mb=10000
#     input:
#         ebqtl="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/{snpset}_variant_gene_pairs.sorted.bed",
#         gtex="results/static_eqtl_followup/gtex/{tissue}.signif_variant_gene_pairs.sorted.bed"
#     output:
#         "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/{snpset}_variant_gene_pairs.{tissue}.window_overlap.bed"
#     shell:
#         """
#         module load bedtools
#         bedtools window -a {input.ebqtl} -b {input.gtex} -w 5000 | awk '$4 == $11' > {output}
#         """
# 
# rule merge_window_sigtests:
#     resources:
#         mem_mb=200000,
#         time="30:00"
#     input:
#         unpack(list_sigtest_window_files)
#     output:
#         "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/{snpset}_variant_gene_pairs.full_gtex_window_overlap.bed"
#     params:
#         tempfile="temp/gtex_overlab_tempfile.txt"
#     shell:
#         """
#         # Combine input files into a single file
#         cat {input} > {params.tempfile}
#         
#         # Filter lines to unique rows
#         sort -u {params.tempfile} > {output}
        # 
        # # Remove temporary file
        # rm {params.tempfile}
        # """

# ### REGULATORY OVERLAP WITH GTEX
# rule filter_and_tidy_gtex_allpairs:
#     """This command filters to entries that have <50kb dist to TSS, and fixes the ENSG ID to drop the version number (after the .)"""
#     input:
#         "../GTEx_Analysis_v8_eQTL_all_associations/{tissue}.allpairs.txt"
#     output:
#         "results/static_eqtl_followup/gtex/allpairs_filtered/{tissue}.allpairs.filtered_50kb_dist2tss.txt"
#     shell:
#         """
#         awk -F"\t" 'BEGIN {{OFS="\t"}} NR == 1 {{print; next}} {{gsub(/\\..*/,"",$1); if (sqrt($3^2) <= 50000) print}}' {input} > {output}
#         """
#         

# rule filter_and_tidy_gtex_sigtests:
#     """This command filters to entries that have <50kb dist to TSS, and fixes the ENSG ID to drop the version number (after the .)"""
#     """Note that gene names were in column 1 for the allpairs files but are in column 2 here""" 
#     input:
#         "../GTEx_Analysis_v8_eQTL/{tissue}.v8.signif_variant_gene_pairs.txt"
#     output:
#         "results/static_eqtl_followup/gtex/sigtests_filtered/{tissue}.signif_variant_gene_pairs.filtered_50kb_dist2tss.txt"
#     shell:
#         """
#         awk -F"\t" 'BEGIN {{OFS="\t"}} NR == 1 {{print; next}} {{gsub(/\\..*/,"",$2); if (sqrt($3^2) <= 50000) print}}' {input} > {output}
#         """
# 
# rule list_all_gtex_tests:
#     resources:
#         mem_mb=100000
#     input:
#         expand("results/static_eqtl_followup/gtex/allpairs_filtered/{tissue}.allpairs.filtered_50kb_dist2tss.txt", 
#         tissue=['Liver', 'Heart_Left_Ventricle', 'Brain_Cortex', 'Cells_Cultured_fibroblasts'])
#     output:
#         "results/static_eqtl_followup/gtex/allpairs_filtered/gtex_all_tests.txt"
#     shell:
#         """
#         tail -q -n +2 {input} | cut -f1,2 | sort -u > {output}
#         """
# 
# rule list_variant_ids_for_gtex:
#     input:
#         "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/eb_gtex_harmonized_tests.txt"
#     output:
#         "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/tested_variants.txt"
#     shell:
#         """
#         tail -n +2 {input} | cut -f 3 > {output}
#         """
  
# rule list_joint_hits:
#   resources:
#       mem_mb=100000
#   input:
#       all_gtex="results/static_eqtl_followup/gtex/allpairs_filtered/{tissue}.allpairs.filtered_50kb_dist2tss.txt",
#       hits_gtex="results/static_eqtl_followup/gtex/sigtests_filtered/{tissue}.signif_variant_gene_pairs.filtered_50kb_dist2tss.txt",
#       all_eb="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/tensorqtl_nominal.all.tsv",
#       hits_eb="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif_variant_gene_pairs.tsv",
#       harmonized_tests="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/eb_gtex_harmonized_tests.txt"
#   output:
#       joint_hit_positions="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif_positions_joint.{tissue}.{celltype}.pos",
#       joint_effects="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/regulatory_overlap.{tissue}.{celltype}.tsv"
#   conda: 
#       "../slurmy/r-mashr.yml"
#   script:
#       "../code/static_eqtl_followup/regulatory_overlap_gtex.R"
