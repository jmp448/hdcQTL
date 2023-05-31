rule make_kinship:
    resources:
        mem_mb=30000,
        time="1:00:00"
    input:
        genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
        inds="data/genotypes/all_individuals_53.tsv"
	  output:
	      expand("data/genotypes/yri_kinship.{ext}", ext=['rel', 'rel.id', 'log', 'nosex'])#,
	      #pruned_bed=expand("data/genotypes/pruned_bed.{ext}", ext=['bed', 'log'])
	  params:
	      kinship_prefix="data/genotypes/yri_kinship",
	      bed_prefix="data/genotypes/pruned_bed"
	  shell:
	      "code/cellregmap_eqtl_calling/make_kinship.sh {input.genotypes} {input.inds} {params.kinship_prefix} {params.bed_prefix}"

# rule run_interaction_test:
  
