import glob

baseline_cats=set(glob_wildcards("data/sldsc/baselineLD_v2.2_bedfiles/{categories}.bed").categories)
all_traits = set(glob_wildcards("TCSC/sumstats/{traits}.sumstats.gz").traits)

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

rule wrangle_geneset:
    resources:
        mem_mb=10000,
        time="10:00"
    input:
        gmt="data/gene_sets/c5.go.bp.v2023.1.Hs.symbols.gmt",
        gtf="data/gencode/gencode.hg38.filtered.gtf"
    output:
        gs="data/ldsc/genesets/{geneset}.tsv"
    conda: 
        "../slurmy/r-sva.yml"
    script:
        "../code/complex_trait_analysis/make_gs_list.R"

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

# HG19 to HG38 Liftover        
# rule sumstats2bed_sa:
#     resources:
#         mem_mb=50000,
#         time="30:00"
#     input:
#         sumstats="{sumstat_file}.sumstats.hg19.tsv.gz"
#     output:
#         bed="{sumstat_file}.sumstats.hg19.bed"
#     conda:
#         "../slurmy/r-sva.yml"
#     script:
#         "../code/complex_trait_analysis/sumstats2bed_sa.R"
        
rule sumstats2bed_ldl:
  # this works for all GWAS catalog summary stats
    resources:
        mem_mb=50000,
        time="30:00"
    input:
        sumstats="{sumstat_file}.sumstats.hg19.tsv.gz"
    output:
        bed="{sumstat_file}.sumstats.hg19.bed"
    params:
        chr_filter=7
    conda:
        "../slurmy/r-sva.yml"
    script:
        "../code/complex_trait_analysis/sumstats2bed_ldl.R"
        
rule sumstats2bed_pardinas:
  # this is for pardinas summary stats downloaded from https://walters.psycm.cf.ac.uk/ 
    resources:
        mem_mb=50000,
        time="30:00"
    input:
        sumstats="data/gwas/pardinas/clozuk_pgc2.meta.sumstats.txt"
    output:
        bed="data/gwas/pardinas/clozuk_pgc2.meta.sumstats.hg19.bed"
    conda:
        "../slurmy/r-sva.yml"
    script:
        "../code/complex_trait_analysis/sumstats2bed_pardinas.R"

    
rule crossmap_hg19_to_hg38:
    resources:
        mem_mb=50000,
        time="30:00"
    input:
        chain_file="data/genotypes/hg19ToHg38.over.chain.gz",
        bed_input="{sumstat_file}.sumstats.hg19.bed"
    output:
        bed_output="{sumstat_file}.sumstats.hg38.bed"
    conda: "../slurmy/genome-toolkit.yml"
    shell:
        """
        python $CONDA_PREFIX/bin/CrossMap.py bed {input.chain_file} {input.bed_input} {output.bed_output}
        """
   
# rule bed2sumstats_sa:
#     resources:
#         mem_mb=50000,
#         time="30:00"
#     input:
#         orig_sumstats="{sumstat_file}.sumstats.hg19.tsv.gz",
#         lifted_bed="{sumstat_file}.sumstats.hg38.bed",
#         annots=expand("data/ldsc/baselineLD_v2.2_annot/baselineLD.{chrom}.annot.gz", chrom=range(1, 23))
#     output:
#         lifted_sumstats="{sumstat_file}.sumstats.hg38.tsv.gz"
#     conda:
#         "../slurmy/r-sva.yml"
#     script:
#         "../code/complex_trait_analysis/bed2sumstats_sa.R"
        
# rule prep_allele_merge:
#     input:
#         snplist="data/ldsc/hg38_regression_snp_weight_files/w_hm3.noMHC.snplist",
#         annots=expand("data/ldsc/baseline_geneset_GOBP-TISSUE-DEV_full_gtexcomb_mash-signif_annot/baselineLD.{chrom}.annot.gz", chrom=range(1, 23))
#     output:
#         allele_list="data/ldsc/hg38_regression_snp_weight_files/w_hm3.noMHC.snplist_withalleles"
#     conda:
#         "../slurmy/r-sva.yml"
#     script:
#         "../code/complex_trait_analysis/prep_allele_merge.R"

rule munge_sumstats_sa:
    resources:
        mem_mb=50000,
        time="30:00"
    input: 
        sumstats="data/gwas/sinnott-armstrong/{trait}.sumstats.hg38.tsv.gz"
    output: 
        "TCSC/sumstats/{trait}.sumstats.gz"
    conda: 
        "../code/ldsc/environment.yml"
    shell:
        """
        python code/ldsc/munge_sumstats.py --sumstats {input.sumstats} \
            --out TCSC/sumstats/{wildcards.trait} 
        """
