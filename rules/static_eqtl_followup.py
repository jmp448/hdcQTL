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
        bedfile="results/static_eqtl_followup/qtl_sets/dynamic-eqtls/crm-signif-tags_variant_gene_pairs.bed"
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

