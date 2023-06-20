#TODO update list_yri_tagged_snps so input file isn't from benchmarking, and isn't IPSC
#TODO update tensorqtl_summary_to_bed_allsigtests to remove the manual removal of NAs
# grab chain file from https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
# write a bash script for this `ls ../GTEx_Analysis_v8_eQTL/*.txt | cut -d '/' -f 3- | cut -d '.' -f 1 | sort | uniq > temp/all_gtex_tissues.txt`
#TODO formalize the following command: `wget -O GTEx_gene_median_tpm.gct.gz https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`

import pandas as pd

def list_sigtest_overlap_files(wildcards):
    tissues = list(pd.read_csv("temp/all_gtex_tissues.txt", names=['tissue'])['tissue'])
    return [f"results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{wildcards.npcs}/{wildcards.snpset}_variant_gene_pairs.{t}.overlap.bed" for t in tissues]
    
rule tensorqtl_summary_to_bed_alltests:
    resources:
        mem_mb=10000,
        time="00:05:00"
    input:
        qtl_summary="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}/tensorqtl_nominal.all.tsv",
        bim_file="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered_plink.bim",
        gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:
        bedfile="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/tensorqtl_alltests.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/tensorqtl_summary_to_bed_alltests.R"
        
rule tensorqtl_summary_to_bed_allsigtests:
    resources:
        mem_mb=10000,
        time="00:05:00"
    input:
        qtl_summary="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif_variant_gene_pairs.tsv",
        bim_file="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/Retinal-cells/genotypes_filtered_plink.bim",
        gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:
        bedfile="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/tensorqtl_summary_to_bed_allsigtests.R"
        
rule list_background_snps_eb_hits_all:
    resources:
        mem_mb=50000,
        time="01:00:00"
    input:
        test_list="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/eb_gtex_harmonized_tests.txt",
        afs="data/genotypes/af_all.frq",
        gtf="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
        all_hits="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif_variant_gene_pairs.tsv"
    output:
        hits_bed="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif-matched_variant_gene_pairs.bed",
        background_bed="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/background-matched_variant_gene_pairs.bed",
        match_details="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif-hits-matcher.rds"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/list_background_snps_eb_hits_all.R"

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
        
rule list_overlap_snps_sigtests_per_tissue:
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
        gtf="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
        rsid_map="data/genotypes/human.YRI.hg38.all.AF.gencode.rsid_snpinfo_map.tsv"
    output:
        harmonized_tests="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/eb_gtex_harmonized_tests.txt"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/harmonize_ebs_to_gtex_all.R"

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

### TRANS EQTL CALLING
rule list_trans_qtl_candidate_variants:
  resources:
      mem_mb=50000
  input:
      tests_list="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/eb_gtex_harmonized_tests.txt",
      afs="data/genotypes/gtex_maf_{tissue}.frq",
      eb_hits="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif_variant_gene_pairs.tsv",
      eb_gtex_overlap="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif_variant_gene_pairs.full_gtex_overlap.bed",
      gtf="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
      #gtf="data/gencode/gencode.hg38.filtered.gtf",
      gmt="data/gene_sets/c5.go.bp.v2022.1.Hs.symbols.gmt"
  output:
      candidate_info="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/trans_eqtl_variant_candidate_info.{gs}.{tissue}.tsv",
      candidates="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/trans_eqtl_variant_candidates.{gs}.{tissue}.txt"
      #match_details="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/trans_eqtl_variant_matchers.{gs}.{tissue}.Rdata"
  conda: 
      "../slurmy/r-mashr.yml"
  script:
      "../code/static_eqtl_followup/list_trans_qtl_candidate_variants.R"

rule list_trans_qtl_candidate_genes:
  resources:
      mem_mb=50000
  input:
      gtf="data/gencode/gencode.hg38.filtered.gtf",
      gmt="data/gene_sets/c5.go.bp.v2022.1.Hs.symbols.gmt"
  output:
      gene_set="data/gene_sets/{gs}.ensg.tsv"
  conda: 
      "../slurmy/r-mashr.yml"
  script:
      "../code/static_eqtl_followup/list_trans_qtl_candidate_genes.R"


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
