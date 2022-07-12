import numpy as np
import pandas as pd

def get_prefix(path, c=1):
    return '/'.join(path.split('/')[:-c])

def list_celltypes(wildcards):
    pb_file = f'data/single_cell_objects/{wildcards.annotation}.{wildcards.aggregation}.tsv'
    samples = list(pd.read_csv(pb_file, sep='\t', nrows=0).columns)[1:]
    pb_clusters = list(np.unique([s.split("_")[1] for s in samples]))
    return pb_clusters

def list_eqtl_files(wildcards):
    pb_file = f'data/single_cell_objects/{wildcards.annotation}.{wildcards.aggregation}.tsv'
    samples = list(pd.read_csv(pb_file, sep='\t', nrows=0).columns)[1:]
    pb_clusters = list(np.unique([s.split("_")[1] for s in samples]))
    return [f"results/static/{wildcards.annotation}/{wildcards.aggregation}/{c}/eqtls.mtc.tsv"
            for c in pb_clusters]

rule pseudobulk_qc:
    resources:
        mem_mb=250000,
        time="15:00"
    input:
        pseudobulk="data/single_cell_objects/{annotation}.{aggregation}.tsv",
        sample_summary="data/static/{annotation}/{aggregation}/sample_summary.tsv"
    output:
        sample_summary_manual="data/static/{annotation}/{aggregation}/{decomp}/sample_summary_manual.tsv"
    params:
        table_prefix = "data/static/{annotation}/{aggregation}/{decomp}",
        fig_prefix = "figs/static/{annotation}/{aggregation}/{decomp}"
    conda: "../slurmy/r-pseudobulk-scran.yml"
    script:
        "../code/mashr/{wildcards.aggregation}-{wildcards.decomp}-qc.R"
        
rule pseudobulk_agg:
    resources:
        mem_mb=250000
    input:
        pseudobulk="data/single_cell_objects/{annotation}.{aggregation}.tsv",
        sample_summary_manual="data/static/{annotation}/{aggregation}/{decomp}/sample_summary_manual.tsv",
        celltypes="data/single_cell_objects/{annotation}.{aggregation}.tsv"
    output:
        all_expression="data/static/{annotation}/{aggregation}/{decomp}/pseudobulk_all.tsv"
    params:
        table_prefix = "data/static/{annotation}/{aggregation}/{decomp}",
        fig_prefix = "figs/static/{annotation}/{aggregation}/{decomp}"
    conda: "../slurmy/r-pseudobulk-scran.yml"
    script:
        "../code/mashr/{wildcards.aggregation}-{wildcards.decomp}-agg.R"

rule genotype_filter:
    input:
	      genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
	      inds="data/static/{annotation}/{aggregation}/{type}/individuals.tsv"
    output:
    	  "data/static/{annotation}/{aggregation}/{type}/genotypes_filtered.recode.vcf"
    params:
        prefix="data/static/{annotation}/{aggregation}/{type}/genotypes_filtered"
    shell:
	      "code/mashr/genotype_filter.sh {input.genotypes} {input.inds} {params.prefix}"

rule genotype_012:
    input:
	      "data/static/{annotation}/{aggregation}/{type}/genotypes_filtered.recode.vcf"
    output:
    	  expand("data/static/{{annotation}}/{{aggregation}}/{{type}}/genotypes_filtered.{out}", out=['012', '012.indv', '012.pos'])
    params:
        prefix="data/static/{annotation}/{aggregation}/{type}/genotypes_filtered"
    shell:
	      "code/mashr/genotype_012.sh {input} {params.prefix}"

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
    conda: "../slurmy/r-matrixEQTL.yml"
    script:
        "../code/mashr/matrixeqtl_nominal.R"

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
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/mashr/matrixeqtl_mtc.R"

rule mashr:
    resources:
        mem_mb=75000,
        time="03:00:00"
    input:
        unpack(list_eqtl_files)
    output:
        trained="results/static/{annotation}/{aggregation}/mashr.trained.rds",
        tophits="results/static/{annotation}/{aggregation}/mashr.tophits.rds",
        random="results/static/{annotation}/{aggregation}/mashr.random_tests.tsv",
        top="results/static/{annotation}/{aggregation}/mashr.top_hits.tsv",
        datasets="results/static/{annotation}/{aggregation}/mashr.training_data.RData"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/mashr/mashr.R"
        
rule ashr:
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
        eqtls="results/static/{annotation}/{aggregation}/Shared/eqtls.mtc.tsv",
        random="results/static/{annotation}/{aggregation}/mashr.random_tests.tsv",
        top="results/static/{annotation}/{aggregation}/mashr.top_hits.tsv",
        df="results/static/{annotation}/{aggregation}/Shared/df.tsv"
    output:
        trained="results/static/{annotation}/{aggregation}/ashr.trained.rds",
        tophits="results/static/{annotation}/{aggregation}/ashr.tophits.rds"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/mashr/ashr.R"