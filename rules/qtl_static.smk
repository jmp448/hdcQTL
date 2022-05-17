"""
rule mashr_pseudobulk_tmm:
    resources:
        mem_mb=250000
    input:
        anndata="data/single_cell_objects/{annotation}.csc.h5ad",
        pseudobulk="data/single_cell_objects/{annotation}.pseudobulk.tsv",
        counts="data/single_cell_objects/{annotation}.pseudobulk.cell_counts.tsv"
    output:
        expand("data/static/{{aggregation}}/type/{{annotation}}/{type}/{out}.tsv", type=CELLTYPES, out=['expression', 'covariates', 'individuals']),
        sample_summary="data/static/{aggregation}/type/{annotation}/sample_summary.tsv",
        celltype_summary="data/static/{aggregation}/type/{annotation}/samples_per_celltype.tsv"
    conda: "slurmy/r-pseudobulk.yml"
    script:
        "code/mashr/{wildcards.aggregation}_static_tmm.R"
"""

rule mashr_pseudobulk_scran:
    resources:
        mem_mb=250000
    input:
        pseudobulk="data/single_cell_objects/{annotation}.pseudobulk.tsv",
        sample_summary="data/static/{aggregation}/type/{annotation}/sample_summary.tsv",
        celltype_summary="data/static/{aggregation}/type/{annotation}/samples_per_celltype.tsv"
    output:
        expand("data/static/{{aggregation}}/type/{{annotation}}/{type}/{out}.tsv", type=CELLTYPES, out=['expression', 'pve', 'covariates', 'individuals'])
    conda: "slurmy/r-pseudobulk-scran.yml"
    script:
        "code/mashr/{wildcards.aggregation}_static_scran.R"


rule mashr_genotype_filter:
    input:
	      genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
	      inds="data/static/{aggregation}/type/{annotation}/{type}/individuals.tsv"
    output:
    	  "data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.recode.vcf"
    shell:
	      "code/mashr/genotype_filter.sh {input.genotypes} {wildcards.aggregation} {wildcards.annotation} {wildcards.type} {input.inds}"

rule mashr_genotype_012:
    input:
	      "data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.recode.vcf"
    output:
    	  expand("data/static/{{aggregation}}/type/{{annotation}}/{{type}}/genotypes_filtered.{out}", out=['012', '012.indv', '012.pos'])
    shell:
	      "code/mashr/genotype_012.sh {input} {wildcards.aggregation} {wildcards.annotation} {wildcards.type}"

rule mashr_genotype_transpose:
    resources:
        mem_mb=50000
    input:
	      "data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.012"
    output:
	      "data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.012.transpose"
    shell:
	      "code/mashr/genotype_transpose.sh {input} {output}"

rule mashr_genotype_reformat:
    resources:
        mem_mb=50000
    input:
        genotypes="data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.012.transpose",
        individuals="data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.012.indv",
        snp_locs="data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.012.pos"
    output:
        snp_locs="data/static/{aggregation}/type/{annotation}/{type}/snp_locs.tsv",
        genotypes="data/static/{aggregation}/type/{annotation}/{type}/genotypes.tsv"
    shell:
        "code/mashr/genotype_reformat.sh {input.genotypes} {input.individuals} {input.snp_locs} {output.snp_locs} {output.genotypes}" 

rule matrix_eqtl:
    resources:
        mem_mb=75000,
        time="00:30:00"
    input:
        genotypes="data/static/pseudobulk/type/{annotation}/{type}/genotypes.tsv",
        snp_locs="data/static/pseudobulk/type/{annotation}/{type}/snp_locs.tsv",
        expression="data/static/pseudobulk/type/{annotation}/{type}/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv",
        covariates="data/static/pseudobulk/type/{annotation}/{type}/covariates.tsv"
    output:
        eqtls="results/mashr/{aggregation}/type/{annotation}/eqtls/{type}/eqtls.tsv",
        df="results/mashr/{aggregation}/type/{annotation}/eqtls/{type}/df.tsv"
    conda: "slurmy/r-matrixEQTL.yml"
    script:
        "code/mashr/matrixEQTL.R"

rule mtc_static:
    resources:
        mem_mb=75000,
        time="00:15:00"
    input:
        eqtls="results/mashr/{aggregation}/type/{annotation}/eqtls/{type}/eqtls.tsv",
        df="results/mashr/{aggregation}/type/{annotation}/eqtls/{type}/df.tsv"
    output:
        all_tests="results/mashr/{aggregation}/type/{annotation}/eqtls/{type}/eqtls.mtc.tsv",
        top_tests="results/mashr/{aggregation}/type/{annotation}/eqtls/{type}/eqtls.tophits.tsv"
    conda: "slurmy/r-mashr.yml"
    script:
        "code/mashr/mtc.R"

rule mashr:
    resources:
        mem_mb=75000,
        time="03:00:00"
    input:
        qtls=expand("results/mashr/{{aggregation}}/type/{{annotation}}/eqtls/{type}/eqtls.mtc.tsv", type=CELLTYPES)
    output:
        trained="results/mashr/{aggregation}/type/{annotation}/eqtls/mashr.trained.rds",
        tophits="results/mashr/{aggregation}/type/{annotation}/eqtls/mashr.tophits.rds"
    conda: "slurmy/r-mashr.yml"
    script:
        "code/mashr/mashr.R"