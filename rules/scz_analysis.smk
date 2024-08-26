import glob

## PREP
### To get LD-tagged SNPs, need to make a BED file
rule make_genotype_bed_gwas_analysis:
    resources:
        mem_mb=100000
    input:
        genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
        inds="data/genotypes/all_individuals_53.tsv"
    output:
        expand("data/genotypes/yri_all.hg38.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/genotypes/yri_all.hg38"
    shell:
        "code/complex_trait_analysis/make_genotype_bed.sh {input.genotypes} {input.inds} {params.prefix}"

rule remove_duplicated_snps:
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
	      expand("data/genotypes/yri_all.hg38.{out}", out=['bed', 'bim', 'fam'])
    output:
    	  expand("data/genotypes/yri_all.hg38.deduplicated.{out}", out=['bed', 'bim', 'fam'])
    params:
        genotypes="data/genotypes/yri_all.hg38",
        temp_allsnp="temp/all_snps",
        temp_dupsnp="temp/duplicated_snps"
    shell:
	      """
	      module load plink
	      plink --bfile {params.genotypes} --write-snplist --out {params.temp_allsnp}
	      cat {params.temp_allsnp}.snplist | sort | uniq -d > {params.temp_dupsnp}.snplist
	      plink --bfile {params.genotypes} --exclude {params.temp_dupsnp}.snplist --make-bed --out {params.genotypes}.deduplicated
	      """


### SCZ ANALYSIS
# Wrangle SCZ summary stats
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
    log: "logs/crossmap_hg19_to_hg38_{sumstat_file}.log"
    shell:
        """
        python $CONDA_PREFIX/bin/CrossMap bed {input.chain_file} {input.bed_input} {output.bed_output} > {log} 2>&1
        """

## List eQTLs that tag a SCZ 
rule generate_scz_qtl_snplist:
    resources:
      mem_mb=100000
    input:
      eb_bed="results/static_eqtl_followup/qtl_sets/tensorqtl/original/signif_variant_gene_pairs.bed",
      gtex_bed="results/static_eqtl_followup/qtl_sets/tensorqtl/original/signif_variant_gene_pairs.all_tissue_overlap.bed",
      scz_sumstats="data/gwas/pardinas/clozuk_pgc2.meta.sumstats.hg38.bed"
    output:
      non_gtex_scz_eqtls="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap/signif_variant_gene_pairs.gtex_removed.snplist.txt",
      non_gtex_scz_eqtls_bed="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap/signif_variant_gene_pairs.gtex_removed.bed"
    conda:
        "../slurmy/r-sva.yml"
    script:
        "../code/complex_trait_analysis/generate_scz_qtl_snplist.R"

rule list_yri_tagged_scz_eqtls:
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
	      genotypes=expand("data/genotypes/yri_all.hg38.deduplicated.{out}", out=['bed', 'bim', 'fam']),
	      snps="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap/signif_variant_gene_pairs.gtex_removed.snplist.txt"
    output:
    	  expand("results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap/signif_variant_gene_pairs.gtex_removed.snplist.{ext}", ext=['tags.list', 'tags', 'nosex', 'log'])
    params:
        genotypes="data/genotypes/yri_all.hg38.deduplicated",
        outfiles="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap/signif_variant_gene_pairs.gtex_removed.snplist"
    shell:
	      """
	      module load plink
	      plink --bfile {params.genotypes} --show-tags {input.snps} --tag-r2 0.5 --tag-kb 1000 --list-all --out {params.outfiles}
	      """

rule expand_bed_to_tag_snps_scz:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        eqtls="results/static_eqtl_followup/qtl_sets/tensorqtl/original/signif_variant_gene_pairs.bed",
        tags="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap/signif_variant_gene_pairs.gtex_removed.snplist.tags.list",
        bim_file="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/all_celltypes_combined/genotypes_filtered_plink.bim"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap/tagging-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/complex_trait_analysis/expand_bed_to_tag_snps.R"
        
rule refine_qtl_set_scz:
    resources:
      mem_mb=100000
    input:
      non_gtex_gwas_eqtls="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap/signif_variant_gene_pairs.gtex_removed.bed",
      gtex_overlap_tag_variants="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap/tagging-signif_variant_gene_pairs.all_tissue_overlap.bed"
    output:
      filtered_bed="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap/signif_variant_gene_pairs.gtex_removed.tagged_gtex_removed.bed"
    conda:
        "../slurmy/r-sva.yml"
    script:
        "../code/complex_trait_analysis/refine_gwas_qtl_snplist.R"

## Trubetskoy SCZ analysis
rule sumstats2bed_trubetskoy:
    resources:
        mem_mb=50000,
        time="30:00"
    input:
        sumstats="data/gwas/trubetskoy_2022/trubetskoy-scz.sumstats.hg19.tsv.gz"
    output:
        bed="data/gwas/trubetskoy_2022/trubetskoy-scz.sumstats.hg19.bed"
    conda:
        "../slurmy/genome-toolkit.yml"
    script:
        "../code/complex_trait_analysis/sumstats2bed_trubetskoy.R"

rule generate_scz_qtl_snplist_trubetskoy:
    resources:
      mem_mb=100000
    input:
      eb_bed="results/static_eqtl_followup/qtl_sets/tensorqtl/original/signif_variant_gene_pairs.bed",
      gtex_bed="/project2/gilad/jpopp/ebQTL/results/static_eqtl_followup/qtl_sets/tensorqtl/original/signif_variant_gene_pairs.all_tissue_overlap.bed",
      scz_sumstats="data/gwas/trubetskoy_2022/trubetskoy-scz.sumstats.hg38.bed"
    output:
      non_gtex_scz_eqtls="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_trubetskoy/signif_variant_gene_pairs.gtex_removed.snplist.txt",
      non_gtex_scz_eqtls_bed="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_trubetskoy/signif_variant_gene_pairs.gtex_removed.bed"
    conda:
        "../slurmy/r-sva.yml"
    script:
        "../code/complex_trait_analysis/generate_scz_qtl_snplist.R"

rule list_yri_tagged_scz_eqtls_trubetskoy:
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
	      genotypes=expand("/project2/gilad/jpopp/ebQTL/data/genotypes/yri_all.hg38.deduplicated.{out}", out=['bed', 'bim', 'fam']),
	      snps="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_trubetskoy/signif_variant_gene_pairs.gtex_removed.snplist.txt"
    output:
    	  expand("results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_trubetskoy/signif_variant_gene_pairs.gtex_removed.snplist.{ext}", ext=['tags.list', 'tags', 'nosex', 'log'])
    params:
        genotypes="data/genotypes/yri_all.hg38.deduplicated",
        outfiles="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_trubetskoy/signif_variant_gene_pairs.gtex_removed.snplist"
    shell:
	      """
	      module load plink
	      plink --bfile {params.genotypes} --show-tags {input.snps} --tag-r2 0.5 --tag-kb 1000 --list-all --out {params.outfiles}
	      """

rule expand_bed_to_tag_snps_scz_trubetskoy:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        eqtls="results/static_eqtl_followup/qtl_sets/tensorqtl/original/signif_variant_gene_pairs.bed",
        tags="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_trubetskoy/signif_variant_gene_pairs.gtex_removed.snplist.tags.list",
        bim_file="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/all_celltypes_combined/genotypes_filtered_plink.bim"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_trubetskoy/tagging-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/complex_trait_analysis/expand_bed_to_tag_snps.R"
        
rule refine_qtl_set_scz_trubetskoy:
    resources:
      mem_mb=100000
    input:
      non_gtex_gwas_eqtls="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_trubetskoy/signif_variant_gene_pairs.gtex_removed.bed",
      gtex_overlap_tag_variants="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_trubetskoy/tagging-signif_variant_gene_pairs.all_tissue_overlap.bed"
    output:
      filtered_bed="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_trubetskoy/signif_variant_gene_pairs.gtex_removed.tagged_gtex_removed.bed"
    conda:
        "../slurmy/r-sva.yml"
    script:
        "../code/complex_trait_analysis/refine_gwas_qtl_snplist.R"

## Willer LDL analysis
rule sumstats2bed_willer:
    resources:
        mem_mb=50000,
        time="30:00"
    input:
        sumstats="data/gwas/willer_2022/willer-ldl.sumstats.hg19.tsv.gz"
    output:
        bed="data/gwas/willer_2022/willer-ldl.sumstats.hg19.bed"
    conda:
        "../slurmy/genome-toolkit.yml"
    script:
        "../code/complex_trait_analysis/sumstats2bed_willer.R"

rule generate_ldl_qtl_snplist_willer:
    resources:
      mem_mb=100000
    input:
      eb_bed="results/static_eqtl_followup/qtl_sets/tensorqtl/original/signif_variant_gene_pairs.bed",
      gtex_bed="/project2/gilad/jpopp/ebQTL/results/static_eqtl_followup/qtl_sets/tensorqtl/original/signif_variant_gene_pairs.all_tissue_overlap.bed",
      scz_sumstats="data/gwas/willer_2022/willer-ldl.sumstats.hg38.bed"
    output:
      non_gtex_scz_eqtls="results/static_eqtl_followup/qtl_sets/tensorqtl/ldl_overlap_willer/signif_variant_gene_pairs.gtex_removed.snplist.txt",
      non_gtex_scz_eqtls_bed="results/static_eqtl_followup/qtl_sets/tensorqtl/ldl_overlap_willer/signif_variant_gene_pairs.gtex_removed.bed"
    conda:
        "../slurmy/r-sva.yml"
    script:
        "../code/complex_trait_analysis/generate_scz_qtl_snplist.R"

rule list_yri_tagged_ldl_eqtls_willer:
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
	      genotypes=expand("/project2/gilad/jpopp/ebQTL/data/genotypes/yri_all.hg38.deduplicated.{out}", out=['bed', 'bim', 'fam']),
	      snps="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_willer/signif_variant_gene_pairs.gtex_removed.snplist.txt"
    output:
    	  expand("results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_willer/signif_variant_gene_pairs.gtex_removed.snplist.{ext}", ext=['tags.list', 'tags', 'nosex', 'log'])
    params:
        genotypes="data/genotypes/yri_all.hg38.deduplicated",
        outfiles="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_willer/signif_variant_gene_pairs.gtex_removed.snplist"
    shell:
	      """
	      module load plink
	      plink --bfile {params.genotypes} --show-tags {input.snps} --tag-r2 0.5 --tag-kb 1000 --list-all --out {params.outfiles}
	      """

rule expand_bed_to_tag_snps_ldl_willer:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        eqtls="results/static_eqtl_followup/qtl_sets/tensorqtl/original/signif_variant_gene_pairs.bed",
        tags="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_willer/signif_variant_gene_pairs.gtex_removed.snplist.tags.list",
        bim_file="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/all_celltypes_combined/genotypes_filtered_plink.bim"
    output:
        bedfile="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_willer/tagging-signif_variant_gene_pairs.bed"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/complex_trait_analysis/expand_bed_to_tag_snps.R"
        
rule refine_qtl_set_scz_willer:
    resources:
      mem_mb=100000
    input:
      non_gtex_gwas_eqtls="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_willer/signif_variant_gene_pairs.gtex_removed.bed",
      gtex_overlap_tag_variants="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_willer/tagging-signif_variant_gene_pairs.all_tissue_overlap.bed"
    output:
      filtered_bed="results/static_eqtl_followup/qtl_sets/tensorqtl/scz_overlap_willer/signif_variant_gene_pairs.gtex_removed.tagged_gtex_removed.bed"
    conda:
        "../slurmy/r-sva.yml"
    script:
        "../code/complex_trait_analysis/refine_gwas_qtl_snplist.R"

rule neur_pseudotime_pseudobulk:
    resources:
        mem_mb=500000,
        time="05:00:00"
    input:
        raw_expression="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/pseudobulk_raw.tsv",
        pseudobulk="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{trajectory}_{nbins}.pseudobulk_tmm.tsv",
        metadata="/project2/gilad/katie/ebQTL/CombinedFormationAndCollectionMetadata_102andPilot_SWAPSANDCONTAMINATIONADDED_012522.csv",
        sample_summary_manual="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/sample_summary_manual.tsv"
    output:
        raw_expression="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/pseudobulk_raw.tsv",
        norm_expression="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/expression.tsv",
        covariates="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/covariates.tsv",
        individuals="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/individuals.tsv",
        pseudotime="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/pseudotime.tsv"
    params:
        table_prefix = "data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/",
        fig_prefix = "figs/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/{pca}/"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/dynamic_qtl_calling/pseudobulk_tmm-agg-dynamic.R"

# baseline_cats=set(glob_wildcards("data/sldsc/baselineLD_v2.2_bedfiles/{categories}.bed").categories)
# all_traits = set(glob_wildcards("TCSC/sumstats/{traits}.sumstats.gz").traits)

# rule prune_ben_baseline:
#     resources:
#         mem_mb=10000,
#         time="10:00"
#     input:
#         annot="data/ldsc/ben_hg38_baselineld_v2.2/baselineLD.{chrom}.annot.gz"
#     output:
#         baseline="data/ldsc/baselineLD_v2.2_annot/baselineLD.{chrom}.annot.gz"
#     conda: 
#         "../slurmy/tensorqtl.yml"
#     script:
#         "../code/complex_trait_analysis/subset_annot_to_baseline.py"
# 
# rule create_qtl_annots:
#     resources:
#         mem_mb=10000,
#         time="10:00"
#     input:
#         baseline="data/ldsc/baselineLD_v2.2_annot/baselineLD.{chrom}.annot.gz",
#         gtex_subset="results/subset_qtl_calling/significant_hits/8pcs/signif_variant_gene_pairs.tsv",
#         eb="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/{variant_group}_variant_gene_pairs.tsv"
#     output:
#         baseline_qtl="data/ldsc/baseline_eqtl_{gtex_grouping}_{variant_group}_annot/baselineLD.{chrom}.annot.gz"
#     conda: 
#         "../slurmy/tensorqtl.yml"
#     script:
#         "../code/complex_trait_analysis/make_eqtl_annots.py"
# 
# rule wrangle_geneset:
#     resources:
#         mem_mb=10000,
#         time="10:00"
#     input:
#         gmt="data/gene_sets/c5.go.bp.v2023.1.Hs.symbols.gmt",
#         gtf="data/gencode/gencode.hg38.filtered.gtf"
#     output:
#         gs="data/ldsc/genesets/{geneset}.tsv"
#     conda: 
#         "../slurmy/r-sva.yml"
#     script:
#         "../code/complex_trait_analysis/make_gs_list.R"

# rule create_qtl_annots_geneset_specific:
#     resources:
#         mem_mb=10000,
#         time="10:00"
#     input:
#         baseline="data/ldsc/baselineLD_v2.2_annot/baselineLD.{chrom}.annot.gz",
#         gtex_subset="results/subset_qtl_calling/significant_hits/8pcs/signif_variant_gene_pairs.tsv",
#         eb="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/{variant_group}_variant_gene_pairs.tsv",
#         geneset="data/ldsc/genesets/{geneset}.tsv"
#     output:
#         baseline_qtl="data/ldsc/baseline_geneset_{geneset}_{subset}_{gtex_grouping}_{variant_group}_annot/baselineLD.{chrom}.annot.gz"
#     conda: 
#         "../slurmy/tensorqtl.yml"
#     script:
#         "../code/complex_trait_analysis/make_geneset_annots.py"
# 
# rule estimate_ld_scores:
#     resources:
#         mem_mb=30000,
#         time="1:00:00"
#     input:
#         annot_file="data/ldsc/{annot_group}/baselineLD.{chrom}.annot.gz",
#         print_snps="data/ldsc/hg38_regression_snp_weight_files/w_hm3.noMHC.snplist"
#     output:
#         expand("data/ldsc/{{annot_group}}/baselineLD.{{chrom}}.{ext}", ext=["log", "l2.M", "l2.M_5_50", "l2.ldscore.gz"])
#     params:
#         plink_prefix="data/ldsc/hg38_plink_files/1000G.EUR.hg38"
#     conda: 
#         "../code/ldsc/environment.yml"
#     shell:
#         """
#         python code/ldsc/ldsc.py \
#             --l2 \
#             --bfile {params.plink_prefix}.{wildcards.chrom} \
#             --ld-wind-cm 1 \
#             --annot {input.annot_file} \
#             --out data/ldsc/{wildcards.annot_group}/baselineLD.{wildcards.chrom} \
#             --print-snps {input.print_snps}
#         """
# 
# rule partition_heritability:
#     resources:
#         mem_mb=30000,
#         time="1:00:00"
#     input:
#         sumstats="TCSC/sumstats/{gwas}.sumstats.gz",
#         ld_scores=expand("data/ldsc/{{annot_group}}/baselineLD.{chrom}.l2.ldscore.gz", chrom=range(1, 23)),
#         weights_files=expand("data/ldsc/hg38_regression_snp_weight_files/weights.hm3_noMHC.{chrom}.l2.ldscore.gz", chrom=range(1, 23)),
#         annot_file=expand("data/ldsc/{{annot_group}}/baselineLD.{chrom}.annot.gz", chrom=range(1, 23))
#     output:
#         log="results/ldsc/{annot_group}/{gwas}.log",
#         results="results/ldsc/{annot_group}/{gwas}.results"
#     params:
#         plink_prefix="data/ldsc/hg38_plink_files/1000G.EUR.hg38",
#         weights_prefix="data/ldsc/hg38_regression_snp_weight_files/weights.hm3_noMHC",
#     conda: 
#         "../code/ldsc/environment.yml"
#     shell:
#         """
#         python code/ldsc/ldsc.py --h2 {input.sumstats} \
#             --ref-ld-chr data/ldsc/{wildcards.annot_group}/baselineLD. \
#             --w-ld-chr {params.weights_prefix}. \
#             --overlap-annot \
#             --frqfile-chr {params.plink_prefix}. \
#             --print-coefficients \
#             --out results/ldsc/{wildcards.annot_group}/{wildcards.gwas}
#         """
# 
# rule sldsc_meta_analysis:
#     input:
#         expand("results/ldsc/{{annot_group}}/{gwas}.results", gwas=all_traits)
#     output:
#         "temp/{annot_group}.done.txt"
#     shell:
#         """
#         echo did it > {output}
#         """

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
        
# rule sumstats2bed_ldl:
#   # this works for all GWAS catalog summary stats
#     resources:
#         mem_mb=50000,
#         time="30:00"
#     input:
#         sumstats="{sumstat_file}.sumstats.hg19.tsv.gz"
#     output:
#         bed="{sumstat_file}.sumstats.hg19.bed"
#     params:
#         chr_filter=7
#     conda:
#         "../slurmy/r-sva.yml"
#     script:
#         "../code/complex_trait_analysis/sumstats2bed_ldl.R"
        
   
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

# rule munge_sumstats_sa:
#     resources:
#         mem_mb=50000,
#         time="30:00"
#     input: 
#         sumstats="data/gwas/sinnott-armstrong/{trait}.sumstats.hg38.tsv.gz"
#     output: 
#         "TCSC/sumstats/{trait}.sumstats.gz"
#     conda: 
#         "../code/ldsc/environment.yml"
#     shell:
#         """
#         python code/ldsc/munge_sumstats.py --sumstats {input.sumstats} \
#             --out TCSC/sumstats/{wildcards.trait} 
#         """
