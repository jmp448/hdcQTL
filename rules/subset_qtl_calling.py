import glob

gtex_tissues = set(glob_wildcards("data/trans_qtl_calling/gtex/expression/{tissues}.v8.normalized_expression.bed.gz").tissues)

rule plink_genotype_reformat_gtex:
    # This command requires GTEx data access, which is protected
    resources:
        mem_mb=100000,
        time="02:00:00"
    input:
        genotypes="/home/jpopp/scratch16-abattle4/lab_data/GTEx_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz"
    output:
        expand("data/subset_qtl_calling/genotypes/plink.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/subset_qtl_calling/genotypes/plink"
    shell:
        "code/subset_qtl_calling/plink_genotype_reformat_gtex.sh {input.genotypes} {params.prefix}"

rule tensorqtl_permutations_subset:
    resources:
        mem_mb=50000,
        time="30:00"
    input:
        exp="/home/jpopp/scratch16-abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_expression_matrices/{tissue}.v8.normalized_expression.bed.gz",
        genotypes=expand("data/subset_qtl_calling/genotypes/plink.{out}", out=['bed', 'bim', 'fam']),
        all_tests="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/eb_gtex_harmonized_tests.txt",
        cov="/home/jpopp/scratch16-abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_covariates/{tissue}.v8.covariates.txt"
    output:
        donor_sample="results/subset_qtl_calling/permutations/{npcs}pcs/{tissue}.donor_sample.tsv",
        cis_df="results/subset_qtl_calling/permutations/{npcs}pcs/{tissue}.tsv"
    params:
        plink_prefix="data/subset_qtl_calling/genotypes/plink"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/subset_qtl_calling/tensorqtl_permutations_subset.py"

rule tensorqtl_merge_subset:
    input:
        expand("results/subset_qtl_calling/permutations/{{npcs}}pcs/{tissue}.tsv", tissue = gtex_tissues)
    output:
        all_qtls="results/subset_qtl_calling/permutations/{npcs}pcs/merged/tensorqtl_permutations.all.tsv",
        top_qtls="results/subset_qtl_calling/permutations/{npcs}pcs/merged/tensorqtl_permutations.sighits.tsv"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/subset_qtl_calling/tensorqtl_merge_subset.R"

rule tensorqtl_nominal_subset:
    resources:
        mem_mb=50000,
        time="01:00:00"
    input:
        genotypes=expand("data/subset_qtl_calling/genotypes/plink.{out}", out=['bed', 'bim', 'fam']),
        all_tests="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/eb_gtex_harmonized_tests.txt",
        exp="/home/jpopp/scratch16-abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_expression_matrices/{tissue}.v8.normalized_expression.bed.gz",
        cov="/home/jpopp/scratch16-abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_covariates/{tissue}.v8.covariates.txt",
        donor_sample="results/subset_qtl_calling/permutations/{npcs}pcs/{tissue}.donor_sample.tsv"
    output:
        expand("results/subset_qtl_calling/nominal/{{npcs}}pcs/{{tissue}}.cis_qtl_pairs.chr{i}.parquet", i=range(1, 23))
    params:
        plink_prefix="data/subset_qtl_calling/genotypes/plink",
        output_prefix="results/subset_qtl_calling/nominal/{npcs}pcs/{tissue}"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/subset_qtl_calling/tensorqtl_nominal_subset.py"

        
rule tensorqtl_merge_nominal_subset:
    resources:
        mem_mb=25000,
        time="30:00"
    input:
        expand("results/subset_qtl_calling/nominal/{{npcs}}pcs/{tissue}.cis_qtl_pairs.chr{i}.parquet", tissue = gtex_tissues,
              i=range(1, 23))
    output:
        merged_df="results/subset_qtl_calling/nominal/{npcs}pcs/merged/tensorqtl_nominal.all.tsv"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/subset_qtl_calling/tensorqtl_merge_nominal_subset.py"