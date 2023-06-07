def strip_tissue:
  

rule list_samples_gtex_tissue:
    input:
        "~/scratch16-abattle4/lab_data/GTEx_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"
    output:
        "~/data-abattle4/jpopp/GTEx_misc/{tissue}_donors.txt"
    params:
        tissue_string="Nerve - Tibial"
    shell:
        """
        cut -f1,14 {input} | grep {params.tissue_string} | cut -d: -f2 | cut -d- -f1,2 | sort -u | tail -n +2 > {output}
        """
        # cut -f1,14 {input} | grep Nerve - Tibial | cut -d: -f2 | cut -d- -f1,2 | sort -u | tail -n +2 > /scratch16/abattle4/jpopp/GTEx_misc/{tissue}_donors.txt
  
rule plink_genotype_reformat_trans:
    # This command requires GTEx data access, which is protected
    resources:
        mem_mb=100000
    input:
        genotypes="~/scratch16-abattle4/lab_data/GTEx_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz",
        inds="~/data-abattle4/jpopp/GTEx_misc/{tissue}_donors.txt",
        trans_candidates="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/{npcs}/trans_eqtl_variant_candidates.{gs}.txt"
    output:
        expand("data/dynamic_qtl_calling/{{trajectory}}_{{nbins}}/pseudobulk_tmm/genotypes_filtered_plink.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/dynamic_qtl_calling/{trajectory}_{nbins}/pseudobulk_tmm/genotypes_filtered_plink"
    shell:
        "code/static_eqtl_followup/plink_genotype_reformat_trans.sh {input.genotypes} {input.inds} {input.trans_candidates} {params.prefix}"
