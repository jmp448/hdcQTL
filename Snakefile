import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()
ZHANG_FILES, = glob_wildcards("data/gwas/zhang/sumstats/{fname}.sumstats")
TRAITS=['TST', 'SHBG']
#CELLTYPES=["IPSC", "MESODERM", "ENDODERM", "EPITHELIAL", "ENDOTHELIAL", "EARLYECTO", "NEUR1", "NEUR2", "NEUR3"]
#CELLTYPES=["EarlyEctoderm","endoderm", "MesodermEndothelialHematopoetic", "NeuralCrest", "PluripotentUndifferentiated", "PostMitoticNeurons", "developingEye", "DevelopingSensoryGanglia", "neuralTube"]
# CELLTYPES=['AFP-ALB-positive-cells', 'Acinar-cells', 'Amacrine-cells',
       # 'Antigen-presenting-cells', 'Astrocytes', 'Bipolar-cells',
       # 'CCL19-CCL21-positive-cells', 'CLC-IL5RA-positive-cells',
       # 'CSH1-CSH2-positive-cells', 'Cardiomyocytes', 'Chromaffin-cells',
       # 'Ciliated-epithelial-cells',
       # 'Corneal-and-conjunctival-epithelial-cells', 'Ductal-cells',
       # 'ELF3-AGBL2-positive-cells', 'ENS-glia', 'ENS-neurons', 'Ectoderm',
       # 'Endocardial-cells', 'Endoderm', 'Epicardial-fat-cells',
       # 'Erythroblasts', 'Excitatory-neurons', 'Extravillous-trophoblasts',
       # 'Ganglion-cells', 'Granule-neurons', 'Hematopoietic-stem-cells',
       # 'Hepatoblasts', 'Horizontal-cells', 'IGFBP1-DKK1-positive-cells',
       # 'Inhibitory-interneurons', 'Inhibitory-neurons',
       # 'Intestinal-epithelial-cells', 'Islet-endocrine-cells',
       # 'Lens-fibre-cells', 'Limbic-system-neurons',
       # 'Lymphatic-endothelial-cells', 'Lymphoid-cells',
       # 'MUC13-DMBT1-positive-cells', 'Megakaryocytes', 'Mesangial-cells',
       # 'Mesoderm', 'Mesothelial-cells', 'Metanephric-cells', 'Microglia',
       # 'Myeloid-cells', 'Neural-crest', 'Neuroendocrine-cells',
       # 'Oligodendrocytes', 'PAEP-MECOM-positive-cells',
       # 'PDE1C-ACSM3-positive-cells', 'PDE11A-FAM19A2-positive-cells',
       # 'Parietal-and-chief-cells', 'Photoreceptor-cells',
       # 'PluripotentUndifferentiated', 'Purkinje-neurons',
       # 'Retinal-pigment-cells', 'Retinal-progenitors-and-Muller-glia',
       # 'SATB2-LRRC7-positive-cells', 'SKOR2-NPSR10-positive-cells',
       # 'SLC24A4-PEX5L-positive-cells', 'Satellite-cells', 'Schwann-cells',
       # 'Skeletal-muscle-cells', 'Smooth-muscle-cells',
       # 'Squamous-epithelial-cells', 'Stellate-cells', 'Stromal-cells',
       # 'Sympathoblasts', 'Syncytiotrophoblasts-and-villous-cytotrophoblasts',
       # 'Thymic-epithelial-cells', 'Thymocytes', 'Trophoblast-giant-cells',
       # 'Unipolar-brush-cells', 'Ureteric-bud-cells',
       # 'Vascular-endothelial-cells', 'Visceral-neurons']

CELLTYPES=['Astrocytes', 'Visceral-neurons', 'Amacrine-cells',
'Oligodendrocytes', 'Ganglion-cells', 'Smooth-muscle-cells', 
'Limbic-system-neurons', 'Vascular-endothelial-cells', 
'Photoreceptor-cells', 'CCL19-CCL21-positive-cells', 
'CLC-IL5RA-positive-cells', 'Mesangial-cells', 'AFP-ALB-positive-cells', 
'Bipolar-cells', 'Satellite-cells', 'ENS-neurons', 'ENS-glia', 
'Granule-neurons', 'Hepatoblasts', 'Microglia', 'IGFBP1-DKK1-positive-cells', 
'Purkinje-neurons', 'Epicardial-fat-cells', 'Cardiomyocytes']

#rule download_gencode_data:
#    input:
#        HTTP.remote("ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.basic.annotation.gff3.gz"", keep_local=True)
#    output:
#        "data/gencode/gencode.v39.basic.annotation.gff3"
#    run:
#      shell("wget -O data/gencode/gencode.v39.basic.annotation.gff3.gz {input}")
#      shell("gzip -d data/gencode/gencode.v39.basic.annotation.gff3.gz")
#
#rule process_gencode_data:
#    input:
#      "data/gencode/gencode.v39.basic.annotation.gff3"
#    output:
#      "data/gencode/gencode.v39.basic.annotation.filtered.gff3",
#      "data/gencode/gencode.v39.basic.annotation.filtered.bed",
#      "data/gencode/gencode.v39.basic.annotation.filtered.tss.bed"
#    script:
#      "code/mashr/process_gff.R"

rule process_gtf:
    input:
        gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf",
        tss_loc="data/gencode/gencode.hg38.filtered.tss.tsv",
        bed_loc="data/gencode/gencode.hg38.filtered.tss.bed"
    script:
        "code/mashr/gene_locs.R"

rule list_snps_zhang:
    input:
        expand("data/gwas/zhang/sumstats/{fname}.sumstats", fname=ZHANG_FILES)
    output:
        "data/gwas/{study}/all_snps.txt"
    run:
        shell("cut -f 1 {input}*.sumstats | sort | uniq | sed -e s/SNP//g > {output}")

rule download_fca_sc:
    input:
        counts=HTTP.remote("atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/gene_count_sampled.RDS", keep_local=True),
        cell_metadata=HTTP.remote("atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/df_cell.RDS", keep_local=True),
        gene_metadata=HTTP.remote("atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/df_gene.RDS", keep_local=True)
    output:
        counts="data/fca/counts.sampled.rds",
        cell_metadata="data/fca/cell_metadata.rds",
        gene_metadata="data/fca/gene_metadata.rds"
    run:
        shell("wget -O {output.counts} {input.counts}")
        shell("wget -O {output.cell_metadata} {input.cell_metadata}")
        shell("wget -O {output.gene_metadata} {input.gene_metadata}")

rule download_fca_celltype:
    input:
        expression=HTTP.remote("https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/FCA_RNA_supp_files/gene_expression_celltype.txt"),
        de_genes=HTTP.remote("https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/DE_gene_77_main_cell_type.csv")
    output:
        expression="data/fca/expression.celltype.csv",
        de_genes="data/fca/de_genes.csv"
    run:
        shell("wget -O {output.expression} {input.expression}")
        shell("wget -O {output.de_genes} {input.de_genes}")

rule create_h5ad_fca:
    input:
        "data/fca/counts.sampled.rds",
        "data/fca/cell_metadata.rds",
        "data/fca/gene_metadata.rds"
    output:
        "data/single_cell_objects/fca.sampled.h5ad",
        "data/fca/organ_celltype.tsv"
    script:
        "code/fca/create_h5ad_fca.R"

rule scdrs_score:
    input:
        "data/single_cell_objects/Lowpass.3seqbatches.merged.scvi_processed_and_scaled.h5ad",
        "data/scDRS/gs_files/{geneset}.gs",
        "data/scDRS/gs_files/{geneset}.map"
    output:
        "results/scDRS/{geneset}/{trait}/scores.tsv",
        "results/scDRS/{geneset}/{trait}/scores.ctrl.tsv"
    script:
        "code/scdrs/scDRS.py"

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
        
rule mashr_no_endothelial:
    resources:
        mem_mb=75000,
        time="01:00:00"
    input:
        qtls=expand("results/mashr/{{aggregation}}/type/{{annotation}}/eqtls/{type}/eqtls.mtc.tsv", type=[t for t in CELLTYPES if t != "ENDOTHELIAL"])
    output:
        trained="results/mashr/{aggregation}/type/{annotation}/eqtls/mashr.trained.noendo.rds",
        tophits="results/mashr/{aggregation}/type/{annotation}/eqtls/mashr.tophits.noendo.rds"
    conda: "slurmy/r-mashr.yml"
    script:
        "code/mashr/mashr.R"
        
rule make_conda: 
    input:
        "test.tmp"
    output:
        "test.tmp2"
    conda: "slurmy/tensorqtl.yml"
    shell:
        "echo booyah > {output}"
        
###
#rule download_gs_files:
#    input:
#        HTTP.remote("ndownloader.figshare.com/files/30853708", keep_local=True)
#    output:
#        gsfile="data/scDRS/gs_files/gs_file.zip"
#    run:
#        shell("wget -O data/scDRS/gs_files/gs_file.zip {input}")
###
