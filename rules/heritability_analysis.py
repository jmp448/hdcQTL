import glob

baseline_cats=set(glob_wildcards("data/sldsc/baselineLD_v2.2_bedfiles/{categories}.bed").categories)
all_traits = set(glob_wildcards("TCSC/sumstats/{traits}.sumstats.gz").traits)
print(all_traits)

rule prune_ben_baseline:
    resources:
        mem_mb=10000,
        time="10:00"
    input:
        annot="data/ldsc/ben_hg38_baselineld_v2.2/baselineLD.{chrom}.annot.gz"
    output:
        baseline="data/ldsc/baselineLD_v2.2_annot/baselineLD.{chrom}.annot.gz"
    conda: 
        "../slurmy/tensorqtl.yml"
    script:
        "../code/complex_trait_analysis/subset_annot_to_baseline.py"

rule create_qtl_annots:
    resources:
        mem_mb=10000,
        time="10:00"
    input:
        baseline="data/ldsc/baselineLD_v2.2_annot/baselineLD.{chrom}.annot.gz",
        gtex_subset="results/subset_qtl_calling/significant_hits/8pcs/signif_variant_gene_pairs.tsv",
        eb="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/{variant_group}_variant_gene_pairs.tsv"
    output:
        baseline_qtl="data/ldsc/baseline_eqtl_{gtex_grouping}_{variant_group}_annot/baselineLD.{chrom}.annot.gz"
    conda: 
        "../slurmy/tensorqtl.yml"
    script:
        "../code/complex_trait_analysis/make_eqtl_annots.py"
    
rule create_qtl_annots_geneset_specific:
    resources:
        mem_mb=10000,
        time="10:00"
    input:
        baseline="data/ldsc/baselineLD_v2.2_annot/baselineLD.{chrom}.annot.gz",
        gtex_subset="results/subset_qtl_calling/significant_hits/8pcs/signif_variant_gene_pairs.tsv",
        eb="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/{variant_group}_variant_gene_pairs.tsv",
        geneset="data/ldsc/genesets/{geneset}.tsv"
    output:
        baseline_qtl="data/ldsc/baseline_geneset_{geneset}_{subset}_{gtex_grouping}_{variant_group}_annot/baselineLD.{chrom}.annot.gz"
    conda: 
        "../slurmy/tensorqtl.yml"
    script:
        "../code/complex_trait_analysis/make_geneset_annots.py"

rule estimate_ld_scores:
    resources:
        mem_mb=30000,
        time="1:00:00"
    input:
        annot_file="data/ldsc/{annot_group}/baselineLD.{chrom}.annot.gz",
        print_snps="data/ldsc/hg38_regression_snp_weight_files/w_hm3.noMHC.snplist"
    output:
        expand("data/ldsc/{{annot_group}}/baselineLD.{{chrom}}.{ext}", ext=["log", "l2.M", "l2.M_5_50", "l2.ldscore.gz"])
    params:
        plink_prefix="data/ldsc/hg38_plink_files/1000G.EUR.hg38"
    conda: 
        "../code/ldsc/environment.yml"
    shell:
        """
        python code/ldsc/ldsc.py \
            --l2 \
            --bfile {params.plink_prefix}.{wildcards.chrom} \
            --ld-wind-cm 1 \
            --annot {input.annot_file} \
            --out data/ldsc/{wildcards.annot_group}/baselineLD.{wildcards.chrom} \
            --print-snps {input.print_snps}
        """

rule partition_heritability:
    resources:
        mem_mb=30000,
        time="1:00:00"
    input:
        sumstats="TCSC/sumstats/{gwas}.sumstats.gz",
        ld_scores=expand("data/ldsc/{{annot_group}}/baselineLD.{chrom}.l2.ldscore.gz", chrom=range(1, 23)),
        weights_files=expand("data/ldsc/hg38_regression_snp_weight_files/weights.hm3_noMHC.{chrom}.l2.ldscore.gz", chrom=range(1, 23)),
        annot_file=expand("data/ldsc/{{annot_group}}/baselineLD.{chrom}.annot.gz", chrom=range(1, 23))
    output:
        log="results/ldsc/{annot_group}/{gwas}.log",
        results="results/ldsc/{annot_group}/{gwas}.results"
    params:
        plink_prefix="data/ldsc/hg38_plink_files/1000G.EUR.hg38",
        weights_prefix="data/ldsc/hg38_regression_snp_weight_files/weights.hm3_noMHC",
    conda: 
        "../code/ldsc/environment.yml"
    shell:
        """
        python code/ldsc/ldsc.py --h2 {input.sumstats} \
            --ref-ld-chr data/ldsc/{wildcards.annot_group}/baselineLD. \
            --w-ld-chr {params.weights_prefix}. \
            --overlap-annot \
            --frqfile-chr {params.plink_prefix}. \
            --print-coefficients \
            --out results/ldsc/{wildcards.annot_group}/{wildcards.gwas}
        """

rule sldsc_meta_analysis:
    input:
        expand("results/ldsc/{{annot_group}}/{gwas}.results", gwas=all_traits)
    output:
        "temp/{annot_group}.done.txt"
    shell:
        """
        echo did it > {output}
        """