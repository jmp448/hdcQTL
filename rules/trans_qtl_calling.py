#TODO update listing of candidate variants to be based on GTEx MAF, not just YRI
#TODO update the snakemake environment so that it has the right conda version without manually installing 4.11.0
#TODO note bandaid - manually changed read_phenotype_bed in tensorqtl
rule list_donors_gtex_tissue:
    input:
        "/home/jpopp/scratch16-abattle4/lab_data/GTEx_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"
    output:
        "data/trans_qtl_calling/gtex/{tissue}_donors.txt"
    params:
        tissue_string="Nerve - Tibial"
    shell:
        """
        cut -f1,14 {input} | grep '{params.tissue_string}' | cut -d: -f2 | cut -d- -f1,2 | sort -u | tail -n +2 > {output}
        """
        # cut -f1,14 {input} | grep Nerve - Tibial | cut -d: -f2 | cut -d- -f1,2 | sort -u | tail -n +2 > /scratch16/abattle4/jpopp/GTEx_misc/{tissue}_donors.txt

rule list_samples_gtex_tissue:
    input:
        "/home/jpopp/scratch16-abattle4/lab_data/GTEx_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"
    output:
        "data/trans_qtl_calling/gtex/{tissue}_samples.txt"
    params:
        tissue_string="Nerve - Tibial"
    shell:
        """
        cut -f1,14 {input} | grep '{params.tissue_string}' | cut -d: -f2 | sort -u | tail -n +2 > {output}
        """

rule get_tissue_maf:
    resources:
        mem_mb=10000,
        time="30:00"
    input:
        genotypes="/home/jpopp/scratch16-abattle4/lab_data/GTEx_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz",
        inds="data/trans_qtl_calling/gtex/{tissue}_donors.txt",
        positions="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/eb_gtex_harmonized_tests.txt"
	  output:
	      "data/genotypes/gtex_maf_{tissue}.frq"
	  params:
	      prefix="data/genotypes/gtex_maf"
	  shell:
	      "code/trans_qtl_calling/compute_af_gtex.sh {input.genotypes} {input.inds} {params.prefix}_{tissue} {input.positions}"

  
rule plink_genotype_reformat_trans:
    # This command requires GTEx data access, which is protected
    resources:
        mem_mb=100000
    input:
        genotypes="/home/jpopp/scratch16-abattle4/lab_data/GTEx_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz",
        inds="data/trans_qtl_calling/gtex/{tissue}_donors.txt",
        trans_candidates="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/trans_eqtl_variant_candidates.{gs}.txt"
    output:
        expand("data/trans_qtl_calling/gtex/genotypes_filtered_plink.{{tissue}}.{{gs}}.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/trans_qtl_calling/gtex/genotypes_filtered_plink"
    shell:
        "code/static_eqtl_followup/plink_genotype_reformat_trans.sh {input.genotypes} {input.inds} {input.trans_candidates} {params.prefix}.{wildcards.tissue}.{wildcards.gs}"

rule tensorqtl_trans:
    resources:
        mem_mb=100000,
        partition="gpu2",
        gres="gpu:1",
        nodes=1
    input:
        genotypes=expand("data/trans_qtl_calling/{{annotation}}/pseudobulk_tmm/basic/{{type}}/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam']),
        exp="data/trans_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/expression.bed.gz",
        cov="data/trans_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/covariates.tsv"
    output:
        expand("results/trans_qtl_calling/{{annotation}}/pseudobulk_tmm/basic/{{type}}/tensorqtl.cis_qtl_pairs.chr{i}.parquet", i=range(1, 23))
    params:
        plink_prefix="data/trans_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/genotypes_filtered_plink",
        output_prefix="results/trans_qtl_calling/{annotation}/pseudobulk_tmm/basic/{type}/tensorqtl"
    conda:
        "../slurmy/tensorqtl.yml"
    script:
        "../code/static_qtl_calling/tensorqtl_nominal.py"
