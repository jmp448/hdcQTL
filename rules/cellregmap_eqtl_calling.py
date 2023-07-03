rule make_genotype_bed:
    resources:
        mem_mb=100000
    input:
        genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
        inds="data/genotypes/all_individuals_53.tsv"
    output:
        expand("data/genotypes/yri_maf0.1_all.hg38.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="data/genotypes/yri_maf0.1_all.hg38"
    shell:
        "code/cellregmap_eqtl_calling/make_genotype_bed.sh {input.genotypes} {input.inds} {params.prefix}"

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

rule wrangle_kinship:
    resources:
        mem_mb=10000,
        time="1:00:00"
    input:
        kinship="data/genotypes/yri_kinship.rel",
        kinship_id="data/genotypes/yri_kinship.rel.id"
    output:
        kinship_tsv="data/genotypes/yri_kinship.tsv"
    conda: 
        "../slurmy/cellregmap.yml"
    script:
        "../code/cellregmap_eqtl_calling/wrangle_kinship.py"

rule compress_expression:
    resources:
        mem_mb=30000,
        time="1:00:00"
    input:
        pseudocell_adata="data/single_cell_objects/eb_pseudocells_raw.h5ad"
    output:
        exp="data/single_cell_objects/eb_pseudocells_normalized.nc"
    conda: 
        "../slurmy/cellregmap.yml"
    script:
        "../code/cellregmap_eqtl_calling/compress_expression.py"

rule make_pseudocell_metadata:
    resources:
        mem_mb=30000,
        time="1:00:00"
    input:
        pseudocell_map="data/fast_topics/cell_pseudocell_mapping.tsv",
        metadata="/project2/gilad/katie/ebQTL/CombinedFormationAndCollectionMetadata_102andPilot_SWAPSANDCONTAMINATIONADDED_012522.csv"
    output:
        pseudocell_metadata="data/cellregmap/pseudocell_metadata.tsv"
    conda: 
        "../slurmy/cellregmap.yml"
    script:
        "../code/cellregmap_eqtl_calling/make_pseudocell_metadata.py"


# rule run_interaction_test:
  
