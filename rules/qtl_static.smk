def get_prefix(path, c=1):
    return '/'.join(path.split('/')[:-c])

rule pseudobulk_qc:
    resources:
        mem_mb=250000
    input:
        pseudobulk="data/single_cell_objects/{annotation}.{aggregation}.tsv",
        sample_summary="data/static/{annotation}/{aggregation}/sample_summary.tsv",
        celltype_summary="data/static/{annotation}/{aggregation}/samples_per_celltype.tsv"
    output:
        tables=expand("data/static/{{annotation}}/{{aggregation}}/{{type}}/{out}.tsv", out=['expression', 'covariates', 'individuals']),
        plots=expand("figs/static/{{annotation}}/{{aggregation}}/{{type}}/{out}.png", out=['scree', 'pcs'])
    params:
        table_prefix = lambda wildcards, output: get_prefix(output.tables[0], c=2),
        fig_prefix = lambda wildcards, output: get_prefix(output.plots[0], c=2)
    conda: "../slurmy/r-pseudobulk-scran.yml"
    script:
        "../code/mashr/{wildcards.aggregation}-qc.R"

rule genotype_filter:
    input:
	      genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
	      inds="data/static/{annotation}/{aggregation}/{type}/individuals.tsv"
    output:
    	  "data/static/{annotation}/{aggregation}/{type}/genotypes_filtered.recode.vcf"
    shell:
	      "code/mashr/genotype_filter.sh {input.genotypes} {wildcards.aggregation} {wildcards.annotation} {wildcards.type} {input.inds}"

rule genotype_012:
    input:
	      "data/static/{annotation}/{aggregation}/{type}/genotypes_filtered.recode.vcf"
    output:
    	  expand("data/static/{{annotation}}/{{aggregation}}/{{type}}/genotypes_filtered.{out}", out=['012', '012.indv', '012.pos'])
    shell:
	      "code/mashr/genotype_012.sh {input} {wildcards.aggregation} {wildcards.annotation} {wildcards.type}"

rule genotype_transpose:
    resources:
        mem_mb=50000
    input:
	      "data/static/{annotation}/{aggregation}/{type}/genotypes_filtered.012"
    output:
	      "data/static/{annotation}/{aggregation}/{type}/genotypes_filtered.012.transpose"
    shell:
	      "code/mashr/genotype_transpose.sh {input} {output}"

rule genotype_reformat:
    resources:
        mem_mb=50000
    input:
        genotypes="data/static/{annotation}/{aggregation}/{type}/genotypes_filtered.012.transpose",
        individuals="data/static/{annotation}/{aggregation}/{type}/genotypes_filtered.012.indv",
        snp_locs="data/static/{annotation}/{aggregation}/{type}/genotypes_filtered.012.pos"
    output:
        snp_locs="data/static/{annotation}/{aggregation}/{type}/snp_locs.tsv",
        genotypes="data/static/{annotation}/{aggregation}/{type}/genotypes.tsv"
    shell:
        "code/mashr/genotype_reformat.sh {input.genotypes} {input.individuals} {input.snp_locs} {output.snp_locs} {output.genotypes}" 

rule matrix_eqtl:
    resources:
        mem_mb=75000,
        time="00:30:00"
    input:
        genotypes="data/static/{annotation}/{aggregation}/{type}/genotypes.tsv",
        snp_locs="data/static/{annotation}/{aggregation}/{type}/snp_locs.tsv",
        expression="data/static/{annotation}/{aggregation}/{type}/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv",
        covariates="data/static/{annotation}/{aggregation}/{type}/covariates.tsv"
    output:
        eqtls="results/static/{annotation}/{aggregation}/{type}/eqtls.tsv",
        df="results/static/{annotation}/{aggregation}/{type}/df.tsv"
    conda: "slurmy/r-matrixEQTL.yml"
    script:
        "code/mashr/matrixEQTL.R"

rule mtc_static:
    resources:
        mem_mb=75000,
        time="00:15:00"
    input:
        eqtls="results/static/{annotation}/{aggregation}/{type}/eqtls.tsv",
        df="results/static/{annotation}/{aggregation}/{type}/df.tsv"
    output:
        all_tests="results/static/{annotation}/{aggregation}/{type}/eqtls.mtc.tsv",
        top_tests="results/static/{annotation}/{aggregation}/{type}/eqtls.tophits.tsv"
    conda: "slurmy/r-mashr.yml"
    script:
        "code/mashr/mtc.R"

rule mashr:
    resources:
        mem_mb=75000,
        time="03:00:00"
    input:
        qtls=expand("results/static/{{annotation}}/{{aggregation}}/{{type}}/eqtls.mtc.tsv")
    output:
        trained="results/static/{annotation}/{aggregation}/mashr.trained.rds",
        tophits="results/static/{annotation}/{aggregation}/mashr.tophits.rds"
    conda: "slurmy/r-mashr.yml"
    script:
        "code/mashr/mashr.R"