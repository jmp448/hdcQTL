import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

include: "rules/qtl_static.smk"

HTTP = HTTPRemoteProvider()

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


rule process_gtf:
    input:
        gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf",
        tss_loc="data/gencode/gencode.hg38.filtered.tss.tsv",
        bed_loc="data/gencode/gencode.hg38.filtered.tss.bed"
    script:
        "code/mashr/gene_locs.R"

        
rule make_conda: 
    input:
        "test.tmp"
    output:
        "test.tmp2"
    conda: "slurmy/tensorqtl.yml"
    shell:
        "echo booyah > {output}"

