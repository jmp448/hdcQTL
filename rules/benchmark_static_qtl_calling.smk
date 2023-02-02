import numpy as np
import pandas as pd
from itertools import product

def get_prefix(path, c=1):
    return '/'.join(path.split('/')[:-c])
    
def get_decomp_type(dc):
    return dc.split('_')[0]

def list_matrixeqtl_outputs(wildcards, npc=None):
    pseudobulk_file = f'data/single_cell_objects/{wildcards.annotation}.{wildcards.aggregation}.tsv'
    samples = list(pd.read_csv(pseudobulk_file, sep='\t', nrows=0).columns)[1:]
    celltypes = list(np.unique([s.split("_")[1] for s in samples]))
    combos = [celltypes, [str(i)+"pcs" for i in range(1, 11)]] # all combos of cell type / n pcs
    celltypes_pcs = ['/'.join(s) for s in product(*combos)]
    file_list = [f"results/static/{wildcards.annotation}/{wildcards.aggregation}/{wildcards.decomp}/{ctpc}/matrixeqtl.cis_qtl_pairs.all.mtc.tsv"
          for ctpc in celltypes_pcs]
    return file_list

def list_matrixeqtl_outputs_npcs(wildcards):
    pseudobulk_file = f'data/single_cell_objects/{wildcards.annotation}.{wildcards.aggregation}.tsv'
    samples = list(pd.read_csv(pseudobulk_file, sep='\t', nrows=0).columns)[1:]
    celltypes = list(np.unique([s.split("_")[1] for s in samples]))
    file_list = [f"results/static/{wildcards.annotation}/{wildcards.aggregation}/{wildcards.decomp}/{c}/{wildcards.npcs}/matrixeqtl.cis_qtl_pairs.all.mtc.tsv"
          for c in celltypes]
    return file_list
    
rule list_all_individuals:
    input:
	      "data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz"
    output:
    	  "data/genotypes/all_individuals.tsv"
    shell:
	      "code/static_qtl_methods_comp/list_all_individuals.sh {input} {output}"

rule pseudobulk_qc:
    resources:
        mem_mb=250000,
        time="1:30:00"
    input:
        pseudobulk="data/benchmark_static_qtl_calling/{annotation}/{annotation}.{aggregation}.tsv",
        sample_summary="data/benchmark_static_qtl_calling/{annotation}/{aggregation}/sample_summary.tsv"
    output:
        sample_summary_manual="data/benchmark_static_qtl_calling/{annotation}/{aggregation}/sample_summary_manual.tsv"
    params:
        table_prefix = "data/static/{annotation}/{aggregation}/{decomp}",
        fig_prefix = "figs/static/{annotation}/{aggregation}/{decomp}",
        decomp_type= get_decomp_type("{decomp}")
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/static_qtl_methods_comp/{wildcards.aggregation}-{wildcards.decomp}-qc.R"
        
rule pseudobulk_agg:
    resources:
        mem_mb=250000,
        time="1:30:00"
    input:
        pseudobulk="data/single_cell_objects/{annotation}.{aggregation}.tsv",
        sample_summary_manual="data/static/{annotation}/{aggregation}/{decomp}/sample_summary_manual.tsv",
        celltypes="data/single_cell_objects/{annotation}.{aggregation}.tsv"
    output:
        all_expression="data/static/{annotation}/{aggregation}/{decomp}/pseudobulk_all.tsv"
    params:
        table_prefix = "data/static/{annotation}/{aggregation}/{decomp}",
        fig_prefix = "figs/static/{annotation}/{aggregation}/{decomp}"
    conda: "../slurmy/r-pseudobulk.yml"
    script:
        "../code/static_qtl_methods_comp/{wildcards.aggregation}-{wildcards.decomp}-agg.R"

rule genotype_filter:
    input:
	      genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
	      inds="data/static/{annotation}/{aggregation}/{decomp}/{type}/individuals.tsv"
    output:
    	  "data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes_filtered.recode.vcf"
    params:
        prefix="data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes_filtered"
    shell:
	      "code/static_qtl_methods_comp/genotype_filter.sh {input.genotypes} {input.inds} {params.prefix}"

rule genotype_012:
    input:
	      "data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes_filtered.recode.vcf"
    output:
    	  expand("data/static/{{annotation}}/{{aggregation}}/{{decomp}}/{{type}}/genotypes_filtered.{out}", out=['012', '012.indv', '012.pos'])
    params:
        prefix="data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes_filtered"
    shell:
	      "code/static_qtl_methods_comp/genotype_012.sh {input} {params.prefix}"

rule genotype_transpose:
    resources:
        mem_mb=50000
    input:
	      "data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes_filtered.012"
    output:
	      "data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes_filtered.012.transpose"
    shell:
	      "code/static_qtl_methods_comp/genotype_transpose.sh {input} {output}"

rule genotype_reformat:
    resources:
        mem_mb=50000
    input:
        genotypes="data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes_filtered.012.transpose",
        individuals="data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes_filtered.012.indv",
        snp_locs="data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes_filtered.012.pos"
    output:
        snp_locs="data/static/{annotation}/{aggregation}/{decomp}/{type}/snp_locs.tsv",
        genotypes="data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes.tsv"
    params:
        temp_loc="temp/genotype_reformat.{annotation}.{aggregation}.{decomp}.{type}"
    shell:
        "code/static_qtl_methods_comp/genotype_reformat.sh {input.genotypes} {input.individuals} {input.snp_locs} {params.temp_loc} {output.snp_locs} {output.genotypes}" 

rule matrixeqtl_nominal:
    resources:
        mem_mb=75000,
        time="00:30:00"
    input:
        genotypes="data/static/{annotation}/{aggregation}/{decomp}/{type}/genotypes.tsv",
        snp_locs="data/static/{annotation}/{aggregation}/{decomp}/{type}/snp_locs.tsv",
        expression="data/static/{annotation}/{aggregation}/{decomp}/{type}/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv",
        covariates="data/static/{annotation}/{aggregation}/{decomp}/{type}/covariates.tsv"
    output:
        eqtls="results/static/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.all.tsv",
        df="results/static/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.df.tsv"
    conda: "../slurmy/r-matrixEQTL.yml"
    script:
        "../code/static_qtl_methods_comp/matrixeqtl_nominal.R"

rule matrixeqtl_mtc:
    resources:
        mem_mb=75000,
        time="00:15:00"
    input:
        eqtls="results/static/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.all.tsv",
        df="results/static/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.df.tsv"
    output:
        all_tests="results/static/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.all.mtc.tsv",
        top_tests="results/static/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.tophits.tsv",
        n_hits="results/static/{annotation}/{aggregation}/{decomp}/{type}/{npcs}/matrixeqtl.cis_qtl_pairs.nhits.tsv"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_qtl_methods_comp/matrixeqtl_mtc.R"

rule run_matrixeqtl_allpcs:
    input:
        unpack(list_matrixeqtl_outputs)
    output:
        "temp/{annotation}.{aggregation}.{decomp}.matrixeqtl.done"
    shell:
        "echo booyah > {output}"
        
rule run_matrixeqtl_npcs:
    input:
        unpack(list_matrixeqtl_outputs_npcs)
    output:
        "temp/{annotation}.{aggregation}.{decomp}.matrixeqtl.{npcs}.done"
    shell:
        "echo booyah > {output}"
  
rule ashr:
    resources:
        mem_mb=50000,
        time="03:00:00"
    input:
        eqtls="results/static/{annotation}/{aggregation}/{decomp}/Shared/matrixeqtl.cis_qtl_pairs.all.mtc.tsv",
        random="results/static/{annotation}/{aggregation}/{decomp}/mashr.random_tests.tsv",
        top="results/static/{annotation}/{aggregation}/{decomp}/mashr.top_hits.tsv",
        df="results/static/{annotation}/{aggregation}/{decomp}/Shared/matrixeqtl.df.tsv"
    output:
        trained="results/static/{annotation}/{aggregation}/{decomp}/ashr.trained.rds",
        tophits="results/static/{annotation}/{aggregation}/{decomp}/ashr.tophits.rds"
    conda: "../slurmy/r-mashr.yml"
    script:
        "../code/static_qtl_methods_comp/ashr.R"