import glob

gtex_tissues = set(glob_wildcards("/project2/gilad/jpopp/GTEx_Analysis_v8_eQTL/{tissues}.v8.signif_variant_gene_pairs.txt").tissues)

# BEDTOOLS HELPERS      
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

# WRANGLE GTEX DATA
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

rule list_overlap_snps_sigtests_per_tissue:
    """
    This command will:
    1. find the SNP intersection of two bedfiles,
    2. filter to SNPs which were tested against the same gene
    """
    resources:
        mem_mb=10000
    input:
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
