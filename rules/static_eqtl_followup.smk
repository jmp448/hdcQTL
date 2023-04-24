#TODO update list_yri_tagged_snps so input file isn't from benchmarking, and isn't IPSC
# grab chain file from https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

rule matrixeqtl_summary_to_bed:
    resources:
        mem_mb=10000,
        time="00:05:00"
    input:
        "results/benchmark_specificity_methods/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.tophits.tsv"
    output:
        "results/benchmark_specificity_methods/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.tophits.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/benchmark_specificity_methods_eqtl_followup/matrixeqtl_summary_to_bed.R"
        
rule matrixeqtl_summary_to_positions:
    resources:
        mem_mb=10000,
        time="00:05:00"
    input:
        "results/benchmark_specificity_methods/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.tophits.tsv"
    output:
        "results/benchmark_specificity_methods/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.tophits.pos"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/matrixeqtl_summary_to_positions.R"

rule tensorqtl_summary_to_bed:
    resources:
        mem_mb=10000,
        time="00:05:00"
    input:
        qtl_summary="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}/tensorqtl_permutations.sighits.tsv",
        bim_file="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered_plink.bim"
    output:
        bedfile="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/tensorqtl_permutations.sighits.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/tensorqtl_summary_to_bed.R"
        
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
        
rule tensorqtl_summary_to_positions:
    resources:
        mem_mb=10000,
        time="00:05:00"
    input:
        qtl_summary="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{npcs}/tensorqtl_permutations.sighits.tsv",
        bim_file="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered_plink.bim"
    output:
        posfile="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/tensorqtl_permutations.sighits.pos"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/tensorqtl_summary_to_positions.R"

rule list_yri_tagged_snps:
# Expand the list of SNPs to all LD tagged SNPs   
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
	      genotypes="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered.recode.vcf",
	      snps="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/tensorqtl_permutations.sighits.pos"
    output:
    	  "results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/yri_tagging.list.geno.ld"
    params:
        "results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/yri_tagging"
    shell:
	      "code/static_eqtl_followup/list_ld_tagged_snps.sh {input.genotypes} {input.snps} {params}"

rule reformat_yri_tagged_snps:
# put the LD block in BED format and list the gene associated with each SNP by ENSG ID
    input:
        ld_block="results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/yri_tagging.list.geno.ld",
        top_hits="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/tensorqtl_permutations.sighits.pos",
        gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
        bim_file="data/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered_plink.bim"
    output:
        "results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/yri_ldblock_tophits.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/reformat_yri_ld_block.R"
        
rule reformat_GTEx_eqtls:
# convert GTEx significant eGene-eVariant pairs to BED format
    input:
        "/project2/gilad/jpopp/GTEx_Analysis_v8_eQTL/{tissue}.v8.signif_variant_gene_pairs.txt"
    output:
        "results/static_eqtl_followup/gtex/{tissue}.signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/reformat_gtex_snps.R"
        
rule reformat_sarkar_eqtls:
# convert Sarkar's single-cell IPSC eGenes to BED format
    input:
        "/project2/gilad/jpopp/data_from_abhishek/mean.txt.gz"
    output:
        "results/static_eqtl_followup/gtex/Sarkar_IPSC.signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/reformat_sarkar_snps.R"
        
rule reformat_i2qtl_eqtls:
# convert I2QTL single-cell IPSC eGenes to BED format
    input:
        "data/i2qtl_paper/i2qtl_IPSC_summary_statistics_egene.tsv"
    output:
        "results/static_eqtl_followup/gtex/I2QTL_IPSC_hg19.signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_eqtl_followup/reformat_i2qtl_snps.R"
        
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

rule crossmap_hg19_to_hg38:
    input:
        chain_file="data/genotypes/hg19ToHg38.over.chain.gz",
        bed_input="results/static_eqtl_followup/gtex/I2QTL_IPSC_hg19.signif_variant_gene_pairs.sorted.bed"
    output:
        bed_output="results/static_eqtl_followup/gtex/I2QTL_IPSC_hg38.signif_variant_gene_pairs.sorted.bed"
    conda: "../slurmy/genome-toolkit.yml"
    shell:
        """
        python $CONDA_PREFIX/bin/CrossMap.py bed {input.chain_file} {input.bed_input} {output.bed_output}
        """

rule estimate_tagged_GTEx_tests:
    """
    This command will:
    1. find the GTEx QTLs within 100kb of an EB test,
    2. filter to SNPs which were tested against the same gene,
    3. filter to lines tagging a unique GTEx eQTL (so large LD blocks with multiple overlaps aren't multiply counted)
    """
    input:
        ebqtl="results/static_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/tensorqtl_alltests.sorted.bed",
        gtex="results/static_eqtl_followup/gtex/{tissue}.signif_variant_gene_pairs.sorted.bed"
    output:
        "results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/n_tests_tagging_gtex.{tissue}.bed"
    shell:
        """
        module load bedtools
        bedtools window -w 10000 -a {input.ebqtl} -b {input.gtex} | awk '$4 == $10' | awk '!seen[$4,$6]++' | awk 'END {{print NR}}' > {output}
        """

rule list_overlap_snps:
    """
    This command will:
    1. find the SNP intersection of two bedfiles,
    2. filter to SNPs which were tested against the same gene,
    3. filter to lines tagging a unique EB eQTL (so large LD blocks with multiple overlaps aren't multiply counted)
    """
    input:
        ebqtl="results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/yri_ldblock_tophits.sorted.bed",
        gtex="results/static_eqtl_followup/gtex/{tissue}.signif_variant_gene_pairs.sorted.bed"
    output:
        "results/static_eqtl_followup/{annotation}/pseudobulk_tmm/basic/{type}/{npcs}/eqtl_overlap.{tissue}.bed"
    shell:
        """
        module load bedtools
        bedtools intersect -a {input.ebqtl} -b {input.gtex} -wa -wb | awk '$4 == $10' | awk '!seen[$4,$6]++' > {output}
        """
