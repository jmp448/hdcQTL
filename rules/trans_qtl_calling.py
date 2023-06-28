#TODO update listing of candidate variants to be based on GTEx MAF, not just YRI
#TODO update the snakemake environment so that it has the right conda version without manually installing 4.11.0
#TODO note bandaid - manually changed read_phenotype_bed in tensorqtl

import pandas as pd

def map_tissue_to_str(t):
    tissue_to_str_map = pd.read_csv("data/trans_qtl_calling/gtex/tissue_to_filename_map.csv")
    s = tissue_to_str_map.loc[tissue_to_str_map['filenaming_str']==t]['tissue_str'].values[0]
    return s

## LIST VARIANTS AND GENES FOR TRANS EQTL CALLING
rule list_donors_gtex_tissue:
    input:
        "/home/jpopp/scratch16-abattle4/lab_data/GTEx_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"
    output:
        "data/trans_qtl_calling/gtex/donors_per_tissue/donors-{tissue}.txt"
    params:
        tissue_string=lambda wildcards: map_tissue_to_str(wildcards.tissue)
    shell:
        """
        cut -f1,14 {input} | grep '{params.tissue_string}' | cut -d: -f2 | cut -d- -f1,2 | sort -u | tail -n +2 > {output}
        """
        # cut -f1,14 {input} | grep Nerve - Tibial | cut -d: -f2 | cut -d- -f1,2 | sort -u | tail -n +2 > /scratch16/abattle4/jpopp/GTEx_misc/{tissue}_donors.txt

rule get_tissue_maf:
    resources:
        mem_mb=10000,
        time="01:30:00"
    input:
        genotypes="/home/jpopp/scratch16-abattle4/lab_data/GTEx_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz",
        inds="data/trans_qtl_calling/gtex/donors_per_tissue/donors-{tissue}.txt",
        snp_ids="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/tested_variants.txt"
    output:
        "data/trans_qtl_calling/gtex/maf/{tissue}.frq"
    params:
        prefix="data/trans_qtl_calling/gtex/maf/"
    shell:
        "code/trans_qtl_calling/compute_af_gtex.sh {input.genotypes} {input.inds} {params.prefix}/{wildcards.tissue} {input.snp_ids}"

rule list_samples_gtex_tissue:
    input:
        "/home/jpopp/scratch16-abattle4/lab_data/GTEx_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"
    output:
        "data/trans_qtl_calling/gtex/samples_per_tissue/samples-{tissue}.txt"
    params:
        tissue_string=lambda wildcards: map_tissue_to_str(wildcards.tissue)
    shell:
        """
        cut -f1,14 {input} | grep '{params.tissue_string}' | cut -d: -f2 | cut -f1 | sort -u | tail -n +2 > {output}
        """
        
rule list_trans_qtl_candidate_variants:
    resources:
        mem_mb=50000
    input:
        tests_list="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/eb_gtex_harmonized_tests.txt",
        afs="data/trans_qtl_calling/gtex/maf/{tissue}.frq",
        eb_hits="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif_variant_gene_pairs.tsv",
        eb_gtex_overlap="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/signif_variant_gene_pairs.full_gtex_overlap.bed",
        #gtf="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
        gtf="data/gencode/gencode.hg38.filtered.gtf",
        gmt="data/gene_sets/c5.go.bp.v2022.1.Hs.symbols.gmt"
    output:
        candidate_info="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/trans_eqtl_variant_candidate_info.{gs}.{tissue}.tsv",
        candidates="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/trans_eqtl_variant_candidates.{gs}.{tissue}.txt"
        #match_details="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/trans_eqtl_variant_matchers.{gs}.{tissue}.Rdata"
    conda: 
        "../slurmy/r-mashr.yml"
    script:
        "../code/trans_qtl_calling/list_trans_qtl_candidate_variants.R"

rule list_trans_qtl_candidate_genes:
    resources:
        mem_mb=50000
    input:
        gtf="data/gencode/gencode.hg38.filtered.gtf",
        gmt="data/gene_sets/c5.go.bp.v2022.1.Hs.symbols.gmt"
    output:
        gene_set="data/gene_sets/{gs}.ensg.tsv"
    conda:
        "../slurmy/r-mashr.yml"
    script:
        "../code/trans_qtl_calling/list_trans_qtl_candidate_genes.R"

rule plink_genotype_reformat_trans:
    # This command requires GTEx data access, which is protected
    resources:
        mem_mb=100000,
        time="02:00:00"
    input:
        genotypes="/home/jpopp/scratch16-abattle4/lab_data/GTEx_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz",
        inds="data/trans_qtl_calling/gtex/donors_per_tissue/donors-{tissue}.txt",
        trans_candidates="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/trans_eqtl_variant_candidates.{gs}.{tissue}.txt"
    output:
        expand("data/trans_qtl_calling/gtex/genotypes_filtered/plink.{{tissue}}.{{gs}}.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/trans_qtl_calling/gtex/genotypes_filtered/plink"
    shell:
        "code/static_eqtl_followup/plink_genotype_reformat_trans.sh {input.genotypes} {input.inds} {input.trans_candidates} {params.prefix}.{wildcards.tissue}.{wildcards.gs}"

rule tensorqtl_trans_allgenes:
    resources:
        mem_mb=100000,
        time="01:00:00"
    input:
        genotypes=expand("data/trans_qtl_calling/gtex/genotypes_filtered/plink.{{tissue}}.{{candidate_gs}}.{out}", out=['bed', 'bim', 'fam']),
        exp="data/trans_qtl_calling/gtex/expression/{tissue}.v8.normalized_expression.bed.gz",
        cov="data/trans_qtl_calling/gtex/covariates/{tissue}.v8.covariates.txt"
    output:
        output_loc="results/trans_qtl_calling/{annotation}/pseudobulk_tmm/basic/{tissue}.{candidate_gs}-variants.all-genes.tsv"
    params:
        plink_prefix="data/trans_qtl_calling/{annotation}/pseudobulk_tmm/basic/genotypes_filtered/plink.{tissue}.{gs}"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/trans_qtl_calling/trans_qtl_calling_allgenes.py"
        
rule tensorqtl_trans_pathwaygenes:
    resources:
        mem_mb=100000,
        time="01:00:00"
    input:
        genotypes=expand("data/trans_qtl_calling/gtex/genotypes_filtered/plink.{{tissue}}.{{candidate_gs}}.{out}", out=['bed', 'bim', 'fam']),
        exp="data/trans_qtl_calling/gtex/expression/{tissue}.v8.normalized_expression.bed.gz",
        cov="data/trans_qtl_calling/gtex/covariates/{tissue}.v8.covariates.txt",
        gs_genes="data/gene_sets/{affected_gs}.ensg.tsv"
    output:
        output_loc="results/trans_qtl_calling/{annotation}/pseudobulk_tmm/basic/{tissue}.{candidate_gs}-variants.{affected_gs}-genes.tsv"
    params:
        plink_prefix="data/trans_qtl_calling/gtex/genotypes_filtered/plink.{tissue}.{candidate_gs}"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/trans_qtl_calling/trans_qtl_calling_pathwaygenes.py"

## CELL TYPE PROPORTION QTL CALLING
rule wrangle_ctprops_per_tissue:
    input:
        tissue_samples="data/trans_qtl_calling/gtex/samples_per_tissue/samples-{tissue}.txt",
        ctprops="data/trans_qtl_calling/gtex/celltype_proportions/celltype_proportions.v8.xCell.7celltypes.txt"
    output:
        tissue_proportions="data/trans_qtl_calling/gtex/celltype_proportions/proportions-{tissue}.txt"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/trans_qtl_calling/wrangle_ctprops_per_tissue.py"

rule tensorqtl_ctprops:
    resources:
        mem_mb=100000,
        time="01:00:00"
    input:
        genotypes=expand("data/trans_qtl_calling/gtex/genotypes_filtered/plink.{{tissue}}.{{candidate_gs}}.{out}", out=['bed', 'bim', 'fam']),
        ctprops="data/trans_qtl_calling/gtex/celltype_proportions/proportions-{tissue}.txt",
        cov="data/trans_qtl_calling/gtex/covariates/{tissue}.v8.covariates.txt"
    output:
        output_loc="results/trans_qtl_calling/{annotation}/pseudobulk_tmm/basic/{tissue}.{candidate_gs}-variants.celltype-proportions.tsv"
    params:
        plink_prefix="data/trans_qtl_calling/gtex/genotypes_filtered/plink.{tissue}.{candidate_gs}"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/trans_qtl_calling/ctprop_qtl_calling.py"
