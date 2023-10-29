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
