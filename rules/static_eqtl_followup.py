#TODO update list_yri_tagged_snps so input file isn't from benchmarking, and isn't IPSC
#TODO update tensorqtl_summary_to_bed_allsigtests to remove the manual removal of NAs
# grab chain file from https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
# write a bash script for this `ls ../GTEx_Analysis_v8_eQTL/*.txt | cut -d '/' -f 3- | cut -d '.' -f 1 | sort | uniq > temp/all_gtex_tissues.txt`
#TODO formalize the following command: `wget -O GTEx_gene_median_tpm.gct.gz https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`
#TODO update the scripts where we list significant hits to NOT subset to GTEx tested variants
#TODO write new scripts for subsetting bed files to gtex variants, and then conducting gtex overlap analysis
#TODO add a clumping step for each time we list the top hits - they all have huge LD blocks
#TODO tidy up so the background for static and dynamic are the same rule
#TODO don't forget the GTEx overlap and everything done on BED files needs to be updated to reflect file reorganization (qtl_sets folder)
#TODO note that ldetect data is from https://bitbucket.org/nygcresearch/ldetect-data/src/master/AFR/
#TODO to run a different trajectory need to update the bim file in the make 1k rule, and need to manually set the trajectory at the bottom of the assoc script
#TODO address use of ancient in crm to bed 

import pandas as pd
import glob

def list_sigtest_overlap_files(wildcards):
    tissues = list(pd.read_csv("temp/all_gtex_tissues.txt", names=['tissue'])['tissue'])
    return [f"results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{wildcards.npcs}/{wildcards.snpset}_variant_gene_pairs.{t}.overlap.bed" for t in tissues]
 
def list_sigtest_window_files(wildcards):
    tissues = list(pd.read_csv("temp/all_gtex_tissues.txt", names=['tissue'])['tissue'])
    return [f"results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{wildcards.npcs}/{wildcards.snpset}_variant_gene_pairs.{t}.window_overlap.bed" for t in tissues]

rule map_rsid_to_snpinfo:
    input:
        vcf_all="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz"
    output:
        rsid_snpinfo_map="data/genotypes/human.YRI.hg38.all.AF.gencode.rsid_snpinfo_map.tsv"
    shell:
        """
        module load bcftools
        bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' {input.vcf_all} | awk -F'\t' '{{OFS="\t"; print $1, $2"_"$3"_"$4"_"$5"_b38"}}' > {output.rsid_snpinfo_map}
        """

rule harmonize_eb_to_gtex_alltests:
    resources:
        mem_mb=100000
    input:
        all_gtex="results/static_eqtl_followup/gtex/allpairs_filtered/gtex_all_tests.txt",
        all_eb="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/tensorqtl_nominal.all.tsv",
        gtf="data/gencode/gencode.hg38.filtered.gtf",
        rsid_map="data/genotypes/human.YRI.hg38.all.AF.gencode.rsid_snpinfo_map.tsv"
    output:
        harmonized_tests="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/eb_gtex_harmonized_tests.txt"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/harmonize_ebs_to_gtex_all.R"

rule tensorqtl_summary_to_bed_alltests:
    resources:
        mem_mb=30000,
        time="15:00"
    input:
        qtl_summary="results/static_qtl_calling/eb-cellid/pseudobulk_tmm/basic/8pcs/tensorqtl_nominal.all.tsv",
        bim_file="data/static_qtl_calling/eb-cellid/pseudobulk_tmm/basic/PNS-glia/genotypes_filtered_plink.bim",
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/tensorqtl-all_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/tensorqtl_summary_to_bed_alltests.R"
        
rule tensorqtl_summary_to_bed_allsigtests:
    resources:
        mem_mb=30000,
        time="15:00"
    input:
        qtl_summary="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/signif_variant_gene_pairs.tsv",
        bim_file="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/Retinal-cells/genotypes_filtered_plink.bim",
        gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/tensorqtl_summary_to_bed_allsigtests.R"

rule mash_to_bed_alltests:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        mash_model="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/mash_fitted_model.full.rds",
        bim_file="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/PNS-glia/genotypes_filtered_plink.bim",
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/mash-all_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/mash_to_bed_alltests.R"

# rule mash_to_bed_tophits:
# This is for if you fit the mash model to just the top hits per gene rather than all
#     resources:
#         mem_mb=10000,
#         time="05:00"
#     input:
#         mash_model="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/mash_fitted_model.tophits.rds",
#         harmonized_tests="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/eb_gtex_harmonized_tests.txt"
#     output:
#         bedfile="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/mash-tophits_variant_gene_pairs.bed"
#     conda: "../slurmy/r-mashr.yml"
#     script:
#         "../code/static_eqtl_followup/mash_to_bed_allsigtests.R"
    
rule mash_to_bed_allsigtests:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        mash_model="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/mash_fitted_model.full.rds",
        harmonized_tests="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/eb_gtex_harmonized_tests.txt"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/mash-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/mash_to_bed_allsigtests.R"

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

rule dynamic_to_bed_alltests:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        neur_tests="results/dynamic_qtl_calling/eb-neur_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.all.tsv",
        cm_tests="results/dynamic_qtl_calling/eb-cm_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.all.tsv",
        hep_tests="results/dynamic_qtl_calling/eb-hep_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.all.tsv",
        bim_file="data/dynamic_qtl_calling/eb-neur_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim",
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/dynamic-all_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/dynamic_to_bed_alltests.R"
        
rule dynamic_to_bed_allhits:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        cm_eqtls="results/dynamic_qtl_calling/eb-cm_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.cis_qtl_all_signif.fdr0.1.tsv",
        neur_eqtls="results/dynamic_qtl_calling/eb-neur_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.cis_qtl_all_signif.fdr0.1.tsv",
        hep_eqtls="results/dynamic_qtl_calling/eb-hep_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.cis_qtl_all_signif.fdr0.1.tsv",
        bim_file="data/dynamic_qtl_calling/eb-neur_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim",
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/original/dynamic-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/dynamic_to_bed_allhits.R"

rule dynamic_to_bed_trajhits:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        eqtls="results/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.cis_qtl_all_signif.fdr0.1.tsv",
        bim_file="data/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim",
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/original/dynamic-{trajectory}-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/dynamic_to_bed_hits.R"
        
# This is currently incorrect, it classifies only the top hits
# rule classify_dynamic_trajhits:
#     resources:
#         mem_mb=65000,
#         time="30:00"
#     input:
#         genotypes="data/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/genotype_df.tsv",
#         qtls="results/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.cis_qtl_top_assoc.txt.gz",
#         pseudotime="data/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/pseudotime.tsv",
#         bim_file="data/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim",
#         gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
#     output:
#         early="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/earlydynamic-{trajectory}-signif_variant_gene_pairs.bed",
#         late="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/latedynamic-{trajectory}-signif_variant_gene_pairs.bed"
#     conda: "../slurmy/r-mashr.yml"
#     script:
#         "../code/static_eqtl_followup/classify_dynamic_eqtls.R"

rule classify_dynamic_trajhits_top:
    resources:
        mem_mb=65000,
        time="30:00"
    input:
        genotypes="data/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/genotype_df.tsv",
        qtls="results/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/nipals/10clpcs/tensorqtl_interactions.cis_qtl_top_assoc.txt.gz",
        pseudotime="data/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/pseudotime.tsv",
        bim_file="data/dynamic_qtl_calling/{trajectory}_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim",
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf"
    output:
        combined="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/original/dynamic-{trajectory}-top_variant_gene_pairs.bed",
        early="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/original/earlydynamic-{trajectory}-top_variant_gene_pairs.bed",
        late="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/original/latedynamic-{trajectory}-top_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/classify_dynamic_eqtls_top.R"

rule get_matched_background:
    resources:
        mem_mb=50000,
        time="01:00:00"
    input:
        tests_bed="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/{tested_snpset}_variant_gene_pairs.bed",
        hits_bed="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/{hit_snpset}_variant_gene_pairs.bed",
        afs="data/genotypes/af_all.frq",
        gtf="data/gencode/gencode.hg38.filtered.gtf"
    output:
        hits_bed="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/{hit_snpset}_variant_gene_pairs.matchable_hits.{tested_snpset}_pool.bed",
        background_bed="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/{hit_snpset}_variant_gene_pairs.matched_background.{tested_snpset}_pool.bed",
        match_details="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/matcher.{hit_snpset}_variant_gene_pairs.{tested_snpset}_pool.rds"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/get_matched_background.R"

rule merge_bims:
    input:
        cm_bim="data/dynamic_qtl_calling/eb-cm_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim",
        hep_bim="data/dynamic_qtl_calling/eb-hep_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim",
        neur_bim="data/dynamic_qtl_calling/eb-neur_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim"
    output:
        "data/dynamic_qtl_calling/all_trajectories_combined/genotypes_filtered_plink.bim"
    conda: "../slurmy/r-mashr.yml"
    shell:
        """
        cat {input} | sort -u -k 1,1 -k 4,4n > {output}
        """
        
rule get_matched_background_1K:
    resources:
        mem_mb=50000,
        time="30:00",
        partition="broadwl"
    input:
        tests_bed="results/static_eqtl_followup/qtl_sets/{analysis}/original/{full_snpset}_variant_gene_pairs.bed",
        hits_bed="results/static_eqtl_followup/qtl_sets/{analysis}/original/{snpset}_variant_gene_pairs.bed",
        afs="data/genotypes/af_all.frq",
        gtf="data/gencode/gencode.hg38.filtered.gtf",
        bim="data/dynamic_qtl_calling/all_trajectories_combined/genotypes_filtered_plink.bim"
    output:
        expand("results/static_eqtl_followup/qtl_sets/{{analysis}}/background_matched/{{snpset}}_{i}.{{full_snpset}}_background.bed", i=range(1000))
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/get_matched_background_1K.R"

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
        ebqtl="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/{snpset}_variant_gene_pairs.sorted.bed",
        gtex="results/static_eqtl_followup/gtex/{tissue}.signif_variant_gene_pairs.sorted.bed"
    output:
        "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/{snpset}_variant_gene_pairs.{tissue}.overlap.bed"
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
        unpack(list_sigtest_overlap_files)
    output:
        "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/{snpset}_variant_gene_pairs.full_gtex_overlap.bed"
    params:
        tempfile="temp/gtex_overlab_tempfile.txt"
    shell:
        """
        # Combine input files into a single file
        cat {input} > {params.tempfile}
        
        # Filter lines to unique rows
        sort -u {params.tempfile} > {output}
        
        # Remove temporary file
        rm {params.tempfile}
        """

rule list_window_snps_sigtests_per_tissue:
    """
    This command will:
    1. find the SNP intersection of two bedfiles,
    2. filter to SNPs which were tested against the same gene
    """
    resources:
        mem_mb=10000
    input:
        ebqtl="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/{snpset}_variant_gene_pairs.sorted.bed",
        gtex="results/static_eqtl_followup/gtex/{tissue}.signif_variant_gene_pairs.sorted.bed"
    output:
        "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/{snpset}_variant_gene_pairs.{tissue}.window_overlap.bed"
    shell:
        """
        module load bedtools
        bedtools window -a {input.ebqtl} -b {input.gtex} -w 5000 | awk '$4 == $11' > {output}
        """

rule merge_window_sigtests:
    resources:
        mem_mb=200000,
        time="30:00"
    input:
        unpack(list_sigtest_window_files)
    output:
        "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/{snpset}_variant_gene_pairs.full_gtex_window_overlap.bed"
    params:
        tempfile="temp/gtex_overlab_tempfile.txt"
    shell:
        """
        # Combine input files into a single file
        cat {input} > {params.tempfile}
        
        # Filter lines to unique rows
        sort -u {params.tempfile} > {output}
        
        # Remove temporary file
        rm {params.tempfile}
        """

### REGULATORY OVERLAP WITH GTEX
rule filter_and_tidy_gtex_allpairs:
    """This command filters to entries that have <50kb dist to TSS, and fixes the ENSG ID to drop the version number (after the .)"""
    input:
        "../GTEx_Analysis_v8_eQTL_all_associations/{tissue}.allpairs.txt"
    output:
        "results/static_eqtl_followup/gtex/allpairs_filtered/{tissue}.allpairs.filtered_50kb_dist2tss.txt"
    shell:
        """
        awk -F"\t" 'BEGIN {{OFS="\t"}} NR == 1 {{print; next}} {{gsub(/\\..*/,"",$1); if (sqrt($3^2) <= 50000) print}}' {input} > {output}
        """
        
rule list_all_gtex_tests:
    resources:
        mem_mb=100000
    input:
        expand("results/static_eqtl_followup/gtex/allpairs_filtered/{tissue}.allpairs.filtered_50kb_dist2tss.txt", 
        tissue=['Liver', 'Heart_Left_Ventricle', 'Brain_Cortex', 'Cells_Cultured_fibroblasts'])
    output:
        "results/static_eqtl_followup/gtex/allpairs_filtered/gtex_all_tests.txt"
    shell:
        """
        tail -q -n +2 {input} | cut -f1,2 | sort -u > {output}
        """

rule filter_and_tidy_gtex_sigtests:
    """This command filters to entries that have <50kb dist to TSS, and fixes the ENSG ID to drop the version number (after the .)"""
    """Note that gene names were in column 1 for the allpairs files but are in column 2 here""" 
    input:
        "../GTEx_Analysis_v8_eQTL/{tissue}.v8.signif_variant_gene_pairs.txt"
    output:
        "results/static_eqtl_followup/gtex/sigtests_filtered/{tissue}.signif_variant_gene_pairs.filtered_50kb_dist2tss.txt"
    shell:
        """
        awk -F"\t" 'BEGIN {{OFS="\t"}} NR == 1 {{print; next}} {{gsub(/\\..*/,"",$2); if (sqrt($3^2) <= 50000) print}}' {input} > {output}
        """

rule list_variant_ids_for_gtex:
    input:
        "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/eb_gtex_harmonized_tests.txt"
    output:
        "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/tested_variants.txt"
    shell:
        """
        tail -n +2 {input} | cut -f 3 > {output}
        """
  
rule list_joint_hits:
  resources:
      mem_mb=100000
  input:
      all_gtex="results/static_eqtl_followup/gtex/allpairs_filtered/{tissue}.allpairs.filtered_50kb_dist2tss.txt",
      hits_gtex="results/static_eqtl_followup/gtex/sigtests_filtered/{tissue}.signif_variant_gene_pairs.filtered_50kb_dist2tss.txt",
      all_eb="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/tensorqtl_nominal.all.tsv",
      hits_eb="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif_variant_gene_pairs.tsv",
      harmonized_tests="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/eb_gtex_harmonized_tests.txt"
  output:
      joint_hit_positions="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif_positions_joint.{tissue}.{celltype}.pos",
      joint_effects="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/regulatory_overlap.{tissue}.{celltype}.tsv"
  conda: 
      "../slurmy/r-mashr.yml"
  script:
      "../code/static_eqtl_followup/regulatory_overlap_gtex.R"

rule report_cor:
  resources:
      mem_mb=15000
  input:
      pos="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif_positions_joint.{tissue}.{celltype}.pos",
      effects="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/regulatory_overlap.{tissue}.{celltype}.tsv"
  output:
      cor="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/regulatory_cor.{tissue}.{celltype}.txt"
  conda: 
      "../slurmy/r-mashr.yml"
  script:
      "../code/static_eqtl_followup/report_regulatory_cor.R"

## GWAS OVERLAP
rule remove_duplicated_snps:
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
	      expand("data/dynamic_qtl_calling/{{trajectory}}_{{nbins}}/pseudobulk_tmm/{{pca}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam'])
    output:
    	  expand("data/dynamic_qtl_calling/{{trajectory}}_{{nbins}}/pseudobulk_tmm/{{pca}}/genotypes_filtered_plink.deduplicated.{out}", out=['bed', 'bim', 'fam'])
    params:
        genotypes="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/genotypes_filtered_plink",
        temp_allsnp="temp/all_snps",
        temp_dupsnp="temp/duplicated_snps"
    shell:
	      """
	      module load plink
	      plink --bfile {params.genotypes} --write-snplist --out {params.temp_allsnp}
	      cat {params.temp_allsnp}.snplist | sort | uniq -d > {params.temp_dupsnp}.snplist
	      plink --bfile {params.genotypes} --exclude {params.temp_dupsnp}.snplist --make-bed --out {params.genotypes}.deduplicated
	      """

rule list_yri_tagged_snps:
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
	      genotypes=expand("data/dynamic_qtl_calling/{{trajectory}}_{{nbins}}/pseudobulk_tmm/{{pca}}/genotypes_filtered_plink.deduplicated.{out}", out=['bed', 'bim', 'fam']),
	      snps="results/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/novel_dynamic_gwas_snplist.txt"
    output:
    	  expand("results/dynamic_qtl_calling/{{trajectory}}_{{nbins}}/pseudobulk_tmm/{{pca}}/novel_dynamic_gwas_snplist.{ext}", ext=['tags.list', 'tags', 'nosex', 'log']) 
    params:
        genotypes="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/genotypes_filtered_plink.deduplicated",
        outfiles="results/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/novel_dynamic_gwas_snplist"
    shell:
	      """
	      module load plink
	      plink --bfile {params.genotypes} --show-tags {input.snps} --tag-r2 0.5 --tag-kb 1000 --list-all --out {params.outfiles}
	      """

rule remove_duplicated_snps_crm:
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
	      expand("data/genotypes/yri_maf0.1_all.hg38.{out}", out=['bed', 'bim', 'fam'])
    output:
    	  expand("data/genotypes/yri_maf0.1_all.hg38.deduplicated.{out}", out=['bed', 'bim', 'fam'])
    params:
        genotypes="data/genotypes/yri_maf0.1_all.hg38",
        temp_allsnp="temp/all_snps",
        temp_dupsnp="temp/duplicated_snps"
    shell:
	      """
	      module load plink
	      plink --bfile {params.genotypes} --write-snplist --out {params.temp_allsnp}
	      cat {params.temp_allsnp}.snplist | sort | uniq -d > {params.temp_dupsnp}.snplist
	      plink --bfile {params.genotypes} --exclude {params.temp_dupsnp}.snplist --make-bed --out {params.genotypes}.deduplicated
	      """

rule list_yri_tagged_snps_crm:
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
	      genotypes=expand("data/genotypes/yri_maf0.1_all.hg38.deduplicated.{out}", out=['bed', 'bim', 'fam']),
	      snps="results/cellregmap_eqtl_calling/eb_cellid/pseudobulk_tmm/basic/novel_crm_gwas_snplist.txt"
    output:
    	  expand("results/cellregmap_eqtl_calling/eb_cellid/pseudobulk_tmm/basic/novel_crm_gwas_snplist.{ext}", ext=['tags.list', 'tags', 'nosex', 'log']) 
    params:
        genotypes="data/genotypes/yri_maf0.1_all.hg38.deduplicated",
        outfiles="results/cellregmap_eqtl_calling/eb_cellid/pseudobulk_tmm/basic/novel_crm_gwas_snplist"
    shell:
	      """
	      module load plink
	      plink --bfile {params.genotypes} --show-tags {input.snps} --tag-r2 0.5 --tag-kb 1000 --list-all --out {params.outfiles}
	      """

rule dynamic_to_bed_tags:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        eqtls="results/dynamic_qtl_calling/{trajectory}_{nbins}/dynamic-{trajectory}-signif_variant_gene_pairs.bed",
        tags="results/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/nipals/novel_dynamic_gwas_snplist.tags.list",
        bim_file="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim"
    output:
        bedfile="results/dynamic_qtl_calling/{trajectory}_{nbins}/dynamic-{trajectory}-signif-tags_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/dynamic_to_bed_tags.R"

rule crm_to_bed_tags:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        eqtls="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/crm-signif_variant_gene_pairs.bed",
        tags="results/cellregmap_eqtl_calling/eb_cellid/pseudobulk_tmm/basic/novel_crm_gwas_snplist.tags.list",
        bim_file="data/genotypes/yri_maf0.1_all.hg38.bim"
    output:
        bedfile="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/crm-signif-tags_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/dynamic_to_bed_tags.R"

rule silly_copy_hack:
  # make sure to fix this for both dynamic and crm eqtls
    input:
        #"results/dynamic_qtl_calling/{trajectory}_15binstrimmed/dynamic-{trajectory}-signif-tags_variant_gene_pairs.bed"
        "results/static_eqtl_followup/qtl_sets/dynamic-eqtls/crm-signif_variant_gene_pairs.bed"
    output:
        #"results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/dynamic-{trajectory}-signif-tags_variant_gene_pairs.bed"
        "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/crm-signif_variant_gene_pairs.bed"
    shell: "cp {input} {output}"

## BOOTSTRAPPING
rule liftover_ldetect:
    resources:
        mem_mb=50000,
        time="30:00"
    input:
        chain_file="data/genotypes/hg19ToHg38.over.chain.gz",
        bed_input="data/ldetect-data/AFR/fourier_ls-all.bed"
    output:
        bed_output="data/genotypes/ldetect_AFR.hg38.bed"
    conda: "../slurmy/genome-toolkit.yml"
    shell:
        """
        python $CONDA_PREFIX/bin/CrossMap.py region {input.chain_file} {input.bed_input} {output.bed_output}
        """

rule list_ld_blocks:
    input:
        eqtls="results/static_eqtl_followup/qtl_sets/{analysis}/original/{snpset}_variant_gene_pairs.sorted.bed",
        ld_blocks="data/genotypes/ldetect_AFR.hg38.sorted.bed"
    output:
        ld_blocks_implicated="results/static_eqtl_followup/qtl_sets/{analysis}/original/{snpset}_ldblocks.bed"
    shell:
        """
        module load bedtools
        bedtools intersect -sorted -a {input.ld_blocks} -b {input.eqtls} -wa -wb > {output}
        """

# rule resample_ld_blocks:
#     input:
#         ld_blocks="results/static_eqtl_followup/qtl_sets/ld_blocks/{snpset}_ldblocks.sorted.bed"
#     output:
#         ld_blocks_resampled=expand("results/static_eqtl_followup/qtl_sets/resampled_ld_blocks/{{snpset}}/{i}.bed", i=range(1000))
#     shell:
#         """
#         module load bedtools
#         bedtools intersect -sorted -a {input.ld_blocks} -b {input.eqtls} -wa | uniq > {output}
#         """

rule block_resample_snps:
    input:
        ld_block_map="results/static_eqtl_followup/qtl_sets/{analysis}/original/{snpset}_ldblocks.bed",
    output:
        expand("results/static_eqtl_followup/qtl_sets/{{analysis}}/resampled/{{snpset}}_{i}.bed", i=range(1000))
    conda: "../slurmy/scvi.yml"
    script:
        "../code/static_eqtl_followup/block_resample.py"

## ROADMAP Enrichments
rule count_tests:
    input:
        tests="results/static_eqtl_followup/qtl_sets/{analysis}/original/{full_snpset}_variant_gene_pairs.bed"
    output:
        "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{full_snpset}_n_tests.txt"
    shell: "awk 'END {{print NR}}' {input.tests} > {output}"

rule count_regulatory_element:
    input:
        annotation="data/roadmap_epigenomes/{annotation}/{annotation}.{tissue}_combined.hg38.bed",
        tests="results/static_eqtl_followup/qtl_sets/{analysis}/original/{full_snpset}_variant_gene_pairs.bed"
    output:
        "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/n_{annotation}_{full_snpset}_background.txt"
    shell: 
        """
        module load bedtools
        bedtools intersect -a {input.tests} -b {input.annotation} -wa | awk 'END {{print NR}}'> {output}
        """

rule count_original_qtls:
    input:
        qtls="results/static_eqtl_followup/qtl_sets/{analysis}/original/{snpset}_variant_gene_pairs.bed"
    output:
        "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_original_n_qtls.txt"
    shell: "awk 'END {{print NR}}' {input.qtls} > {output}"

rule count_resampled_qtls:
    input:
        qtls="results/static_eqtl_followup/qtl_sets/{analysis}/resampled/{snpset}_{i}.bed"
    output:
        "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_{i}_n_qtls.txt"
    shell: "awk 'END {{print NR}}' {input.qtls} > {output}"

# rule count_matched_qtls:
#     input:
#         full_snpset="results/static_eqtl_followup/qtl_sets/{analysis}/original/{full_snpset}_variant_gene_pairs.bed",
#         qtls="results/static_eqtl_followup/qtl_sets/{analysis}/background_matched/{snpset}_{i}.{full_snpset}_background.bed" 
#     output:
#         "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_{i}.{full_snpset}_background.n_qtls.txt"
#     shell: "awk 'END {{print NR}}' {input.qtls} > {output}"
    
rule count_original_qtl_enhancer_overlap:
    input:
        annotation="data/roadmap_epigenomes/{annotation}/{annotation}.{tissue}_combined.hg38.bed",
        qtls="results/static_eqtl_followup/qtl_sets/{analysis}/original/{snpset}_variant_gene_pairs.bed"
    output:
        "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_original_n_overlaps.txt"
    shell:
        """
        module load bedtools
        bedtools intersect -a {input.qtls} -b {input.annotation} -wa | awk 'END {{print NR}}'> {output}
        """

rule count_resampled_qtl_enhancer_overlap:
    input:
        annotation="data/roadmap_epigenomes/{annotation}/{annotation}.{tissue}_combined.hg38.bed",
        qtls="results/static_eqtl_followup/qtl_sets/{analysis}/resampled/{snpset}_{i}.bed"
    output:
        "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_{i}_n_overlaps.txt"
    shell:
        """
        module load bedtools
        bedtools intersect -a {input.qtls} -b {input.annotation} -wa | awk 'END {{print NR}}'> {output}
        """
        
rule count_matched_qtl_enhancer_overlap:
    input:
        annotation="data/roadmap_epigenomes/{annotation}/{annotation}.{tissue}_combined.hg38.bed",
        full_snpset="results/static_eqtl_followup/qtl_sets/{analysis}/original/{full_snpset}_variant_gene_pairs.bed",
        qtls="results/static_eqtl_followup/qtl_sets/{analysis}/background_matched/{snpset}_{i}.{full_snpset}_background.bed"
    output:
        "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_{i}.{full_snpset}_background.n_overlaps.txt"
    shell:
        """
        module load bedtools
        bedtools intersect -a {input.qtls} -b {input.annotation} -wa | awk 'END {{print NR}}'> {output}
        """

rule prepare_or_inputs:
    input:
        n_tests="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{full_snpset}_n_tests.txt",
        n_elements="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/n_{annotation}_{full_snpset}_background.txt",
        n_qtls="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_original_n_qtls.txt",
        n_overlap="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_original_n_overlaps.txt"
    output:
        "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_qtls_{full_snpset}_background.txt"
    shell:
        """
        echo ready to run > {output}
        """

rule combine_or_confidence_interval_inputs:
    input:
        n_qtls=expand("results/static_eqtl_followup/roadmap_enrichments/{{analysis}}/{{annotation}}_{{tissue}}/{{snpset}}_{i}_n_qtls.txt", i=range(1000)),
        n_overlap=expand("results/static_eqtl_followup/roadmap_enrichments/{{analysis}}/{{annotation}}_{{tissue}}/{{snpset}}_{i}_n_overlaps.txt", i=range(1000))
    output:
        qtl_counts_combined="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_all_qtl_counts.txt",
        overlap_counts_combined="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_all_overlap_counts.txt"
    shell:
        """
        cat {input.n_qtls} > {output.qtl_counts_combined}
        cat {input.n_overlap} > {output.overlap_counts_combined}
        """

rule compile_or_confidence_interval_inputs:
    input:
        n_tests="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{full_snpset}_n_tests.txt",
        n_elements="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/n_{annotation}_{full_snpset}_background.txt",
        n_qtls="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_all_qtl_counts.txt",
        n_overlap="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_all_overlap_counts.txt"
    output: 
      "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_qtls_{full_snpset}_background_bootstrapped.txt"
    shell: "echo ready for analysis > {output}"
    
rule combine_or_matched_ci_inputs:
    input:
        #n_qtls=expand("results/static_eqtl_followup/roadmap_enrichments/{{analysis}}/{{annotation}}_{{tissue}}/{{snpset}}_{i}.{{full_snpset}}_background.n_qtls.txt", i=range(1000)),
        n_overlap=expand("results/static_eqtl_followup/roadmap_enrichments/{{analysis}}/{{annotation}}_{{tissue}}/{{snpset}}_{i}.{{full_snpset}}_background.n_overlaps.txt", i=range(1000))
    output:
        #qtl_counts_combined="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_all_qtl_counts.{full_snpset}_background.txt",
        overlap_counts_combined="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_all_overlap_counts.{full_snpset}_background.txt"
    shell:
        """
        cat {input.n_overlap} > {output.overlap_counts_combined}
        """
        # removed line from above # cat {input.n_qtls} > {output.qtl_counts_combined}
        
rule compile_or_matched_ci_inputs:
    input:
        n_tests="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{full_snpset}_n_tests.txt",
        n_elements="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/n_{annotation}_{full_snpset}_background.txt",
        # n_qtls="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_all_qtl_counts.{full_snpset}_background.txt",
        n_overlap="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_all_overlap_counts.{full_snpset}_background.txt"
    output: 
      "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_qtls_{full_snpset}_background_matched.txt"
    shell: "echo ready for analysis > {output}"
       
# ## Covariate matching version of enrichment analysis
# rule count_bg:
#     input:
#         tests="results/static_eqtl_followup/qtl_sets/{analysis}/original/{full_snpset}_variant_gene_pairs.bed"
#     output:
#         "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{full_snpset}_n_tests.txt"
#     shell: "awk 'END {{print NR}}' {input.tests} > {output}"
# 
# rule count_regulatory_element:
#     input:
#         annotation="data/roadmap_epigenomes/{annotation}/{annotation}.{tissue}_combined.hg38.bed",
#         tests="results/static_eqtl_followup/qtl_sets/{analysis}/original/{full_snpset}_variant_gene_pairs.bed"
#     output:
#         "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/n_{annotation}_{full_snpset}_background.txt"
#     shell: 
#         """
#         module load bedtools
#         bedtools intersect -a {input.tests} -b {input.annotation} -wa | awk 'END {{print NR}}'> {output}
#         """
# 
# rule count_original_qtls:
#     input:
#         qtls="results/static_eqtl_followup/qtl_sets/{analysis}/original/{snpset}_variant_gene_pairs.bed"
#     output:
#         "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_original_n_qtls.txt"
#     shell: "awk 'END {{print NR}}' {input.qtls} > {output}"
# 
# rule count_resampled_qtls:
#     input:
#         qtls="results/static_eqtl_followup/qtl_sets/{analysis}/resampled/{snpset}_{i}.bed"
#     output:
#         "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_{i}_n_qtls.txt"
#     shell: "awk 'END {{print NR}}' {input.qtls} > {output}"
# 
# rule count_original_qtl_enhancer_overlap:
#     input:
#         annotation="data/roadmap_epigenomes/{annotation}/{annotation}.{tissue}_combined.hg38.bed",
#         qtls="results/static_eqtl_followup/qtl_sets/{analysis}/original/{snpset}_variant_gene_pairs.bed"
#     output:
#         "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_original_n_overlaps.txt"
#     shell:
#         """
#         module load bedtools
#         bedtools intersect -a {input.qtls} -b {input.annotation} -wa | awk 'END {{print NR}}'> {output}
#         """
# 
# rule count_resampled_qtl_enhancer_overlap:
#     input:
#         annotation="data/roadmap_epigenomes/{annotation}/{annotation}.{tissue}_combined.hg38.bed",
#         qtls="results/static_eqtl_followup/qtl_sets/{analysis}/resampled/{snpset}_{i}.bed"
#     output:
#         "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_{i}_n_overlaps.txt"
#     shell:
#         """
#         module load bedtools
#         bedtools intersect -a {input.qtls} -b {input.annotation} -wa | awk 'END {{print NR}}'> {output}
#         """
# 
# rule prepare_or_inputs:
#     input:
#         n_tests="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{full_snpset}_n_tests.txt",
#         n_elements="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/n_{annotation}_{full_snpset}_background.txt",
#         n_qtls="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_original_n_qtls.txt",
#         n_overlap="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_original_n_overlaps.txt"
#     output:
#         "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_qtls_{full_snpset}_background.txt"
#     shell:
#         """
#         echo ready to run > {output}
#         """
# 
# rule combine_or_confidence_interval_inputs:
#     input:
#         n_qtls=expand("results/static_eqtl_followup/roadmap_enrichments/{{analysis}}/{{annotation}}_{{tissue}}/{{snpset}}_{i}_n_qtls.txt", i=range(1000)),
#         n_overlap=expand("results/static_eqtl_followup/roadmap_enrichments/{{analysis}}/{{annotation}}_{{tissue}}/{{snpset}}_{i}_n_overlaps.txt", i=range(1000))
#     output:
#         qtl_counts_combined="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_all_qtl_counts.txt",
#         overlap_counts_combined="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_all_overlap_counts.txt"
#     shell:
#         """
#         cat {input.n_qtls} > {output.qtl_counts_combined}
#         cat {input.n_overlap} > {output.overlap_counts_combined}
#         """
# 
# rule compile_or_confidence_interval_inputs:
#     input:
#         n_tests="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{full_snpset}_n_tests.txt",
#         n_elements="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/n_{annotation}_{full_snpset}_background.txt",
#         n_qtls="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_all_qtl_counts.txt",
#         n_overlap="results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_all_overlap_counts.txt"
#     output: 
#       "results/static_eqtl_followup/roadmap_enrichments/{analysis}/{annotation}_{tissue}/{snpset}_qtls_{full_snpset}_background_bootstrapped.txt"
#     shell: "echo ready for analysis > {output}"

# Scripts used for the LD tagging analysis
# 
# rule matrixeqtl_summary_to_bed:
#     resources:
#         mem_mb=10000,
#         time="00:05:00"
#     input:
#         "results/benchmark_specificity_methods/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.tophits.tsv"
#     output:
#         "results/benchmark_specificity_methods/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.tophits.bed"
#     conda: "../slurmy/r-mashr.yml"
#     script:
#         "../code/benchmark_specificity_methods_eqtl_followup/matrixeqtl_summary_to_bed.R"
#         
# rule matrixeqtl_summary_to_positions:
#     resources:
#         mem_mb=10000,
#         time="00:05:00"
#     input:
#         "results/benchmark_specificity_methods/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.tophits.tsv"
#     output:
#         "results/benchmark_specificity_methods/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.tophits.pos"
#     conda: "../slurmy/r-mashr.yml"
#     script:
#         "../code/static_eqtl_followup/matrixeqtl_summary_to_positions.R"
# 
# rule tensorqtl_summary_to_bed:
#     resources:
#         mem_mb=10000,
#         time="00:05:00"
#     input:
#         qtl_summary="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}/tensorqtl_permutations.sighits.tsv",
#         bim_file="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered_plink.bim",
#         gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
#     output:
#         bedfile="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/tensorqtl_permutations.sighits.bed"
#     conda: "../slurmy/r-mashr.yml"
#     script:
#         "../code/static_eqtl_followup/tensorqtl_summary_to_bed.R"
#         
# rule tensorqtl_summary_to_positions:
#     resources:
#         mem_mb=10000,
#         time="00:05:00"
#     input:
#         qtl_summary="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}/tensorqtl_permutations.sighits.tsv",
#         bim_file="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered_plink.bim"
#     output:
#         posfile="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/tensorqtl_permutations.sighits.pos"
#     conda: "../slurmy/r-mashr.yml"
#     script:
#         "../code/static_eqtl_followup/tensorqtl_summary_to_positions.R"
# 
# rule list_yri_tagged_snps:
# # Expand the list of SNPs to all LD tagged SNPs   
#     resources:
#         mem_mb=50000,
#         time="03:00:00"
#     input:
# 	      genotypes="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered.recode.vcf",
# 	      snps="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/tensorqtl_permutations.sighits.pos"
#     output:
#     	  "results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/yri_tagging.list.geno.ld"
#     params:
#         "results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/yri_tagging"
#     shell:
# 	      "code/static_eqtl_followup/list_ld_tagged_snps.sh {input.genotypes} {input.snps} {params}"
# 
# rule reformat_yri_tagged_snps:
# # put the LD block in BED format and list the gene associated with each SNP by ENSG ID
#     input:
#         ld_block="results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/yri_tagging.list.geno.ld",
#         top_hits="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/tensorqtl_permutations.sighits.pos",
#         gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
#         bim_file="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered_plink.bim"
#     output:
#         "results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/yri_ldblock_tophits.bed"
#     conda: "../slurmy/r-mashr.yml"
#     script:
#         "../code/static_eqtl_followup/reformat_yri_ld_block.R"
#         
# rule reformat_sarkar_eqtls:
# # convert Sarkar's single-cell IPSC eGenes to BED format
#     input:
#         "/project2/gilad/jpopp/data_from_abhishek/mean.txt.gz"
#     output:
#         "results/static_eqtl_followup/gtex/Sarkar_IPSC.signif_variant_gene_pairs.bed"
#     conda: "../slurmy/r-mashr.yml"
#     script:
#         "../code/static_eqtl_followup/reformat_sarkar_snps.R"
#         
# rule reformat_i2qtl_eqtls:
# # convert I2QTL single-cell IPSC eGenes to BED format
#     input:
#         "data/i2qtl_paper/i2qtl_IPSC_summary_statistics_egene.tsv"
#     output:
#         "results/static_eqtl_followup/gtex/I2QTL_IPSC_hg19.signif_variant_gene_pairs.bed"
#     conda: "../slurmy/r-mashr.yml"
#     script:
#         "../code/static_eqtl_followup/reformat_i2qtl_snps.R"
# 
# rule crossmap_hg19_to_hg38:
#     input:
#         chain_file="data/genotypes/hg19ToHg38.over.chain.gz",
#         bed_input="results/static_eqtl_followup/gtex/I2QTL_IPSC_hg19.signif_variant_gene_pairs.sorted.bed"
#     output:
#         bed_output="results/static_eqtl_followup/gtex/I2QTL_IPSC_hg38.signif_variant_gene_pairs.sorted.bed"
#     conda: "../slurmy/genome-toolkit.yml"
#     shell:
#         """
#         python $CONDA_PREFIX/bin/CrossMap.py bed {input.chain_file} {input.bed_input} {output.bed_output}
#         """
# 
# rule estimate_tagged_GTEx_tests:
#     """
#     This command will:
#     1. find the GTEx QTLs within 100kb of an EB test,
#     2. filter to SNPs which were tested against the same gene,
#     3. filter to lines tagging a unique GTEx eQTL (so large LD blocks with multiple overlaps aren't multiply counted)
#     """
#     input:
#         ebqtl="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/tensorqtl_alltests.sorted.bed",
#         gtex="results/static_eqtl_followup/gtex/{tissue}.signif_variant_gene_pairs.sorted.bed"
#     output:
#         "results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/n_tests_tagging_gtex.{tissue}.bed"
#     shell:
#         """
#         module load bedtools
#         bedtools window -w 10000 -a {input.ebqtl} -b {input.gtex} | awk '$4 == $10' | awk '!seen[$4,$6]++' | awk 'END {{print NR}}' > {output}
#         """
# 
# Stuff from misunderstanding how harmonization would work
# rule list_rsids_alltests:
#     input:
#         eb_alltests="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/tensorqtl_nominal.all.tsv"
#     output:
#         rsid_list="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/rsids.alltests.tsv"
#     shell:
#         """
#         cut -f 2 {input.eb_alltests} | tail -n +2 | sort -u > {output.rsid_list}
#         """
# 
# rule list_snpinfo_alltests:
#     input:
#         rsid_list="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/rsids.alltests.tsv",
#         vcf_all="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz"
#     output:
#         snp_info="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}pcs/snp_info.alltests.tsv"
#     shell:
#         """
#         module load bcftools
#         bcftools view -R {input.rsid_list} -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' {input.vcf_all} > {output.snp_info}
#         """
# rule list_overlap_snps:
#     """
#     This command will:
#     1. find the SNP intersection of two bedfiles,
#     2. filter to SNPs which were tested against the same gene,
#     3. filter to lines tagging a unique EB eQTL (so large LD blocks with multiple overlaps aren't multiply counted)
#     """
#     input:
#         ebqtl="results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/yri_ldblock_tophits.sorted.bed",
#         gtex="results/static_eqtl_followup/gtex/{tissue}.signif_variant_gene_pairs.sorted.bed"
#     output:
#         "results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/eqtl_overlap.{tissue}.bed"
#     shell:
#         """
#         module load bedtools
#         bedtools intersect -a {input.ebqtl} -b {input.gtex} -wa -wb | awk '$4 == $10' | awk '!seen[$4,$6]++' > {output}
#         """
