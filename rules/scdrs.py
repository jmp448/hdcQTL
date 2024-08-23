all_traits_zhang=['ADHD', 'AD', 'AF', 'BIP', 'CD', 'CEL', 'CAD', 'DPW', 'FG', 'RTOL', 'IBD', 'ISI', 'INT', 'SLE', 'MDD', 'MS', 'CIR', 'RCTT', 'RA', 'SWB', 'SD', 'T1D', 'T2D', 'UC', 'VNR', 'WRY', 'AA', 'AP', 'CHOL', 'GL', 'HDL', 'HBA1C', 'LDL', 'SHBG', 'TST', 'TB', 'TP', 'TG', 'EC', 'LYM', 'MCH', 'MONO', 'PLT', 'RDW', 'RBCC', 'WBCC', 'BMD-HT', 'BALD', 'BMI', 'HGHT', 'WHR', 'DBP', 'SBP', 'BC', 'ECOL', 'EYR', 'SMOK', 'AIT', 'ECZ', 'ASM', 'CVD', 'HTN', 'HT', 'RR-ENT', 'BMR', 'FEV', 'FVC', 'NTC', 'MORN', 'HAIR', 'AAM', 'AAMP', 'CHIL']
oconnor_traits=['CHIL', 'SMOK', 'ECOL', 'MORN', 'NTC', 'FVC', 'CVD', 'SBP', 'BMI', 'IBD', 'HGHT', 'WHR', 'AAM', 'FEV', 'WBCC', 'ECZ', 'ASM', 'EC', 'RA', 'PLT', 'AIT', 'BMD-HT', 'AD', 'T2D', 'RDW', 'RBCC', 'BALD', 'AAMP']
subset_traits=["ADHD", "AD", "AF", "BIP", "CAD", "FG", "ISI", "MDD", "CIR", "RCTT", "SCZ", "SD", "T2D", "AA", "AP", "CHOL", "GL", "HDL", "HBA1C", "LDL", "SHBG", "TST", "TB", "TP", "TG", "BMD-HT", "BMI", "HGHT", "WHR", "DBP", "SBP", "BC", "SMOK", "CVD", "HTN", "HT", "RR-ENT", "BMR", "FEV", "FVC"]

rule scdrs:
    resources:
        mem_mb=200000,
        time="12:00:00"
    input:
        adata="data/trajectory_inference/eb_scdrs_preprocessed.adata",
        gs="data/gene_sets/zhang.magma_10kb_1000.74_traits.gs",
        gs_map="data/gene_sets/zhang.magma_10kb_1000.74_traits.map"
    output:
        "results/scDRS/zhang.magma_10kb_1000.74_traits/{trait}/scores.tsv"
    conda:
        "../slurmy/scvi.yml"
    script:
        "../code/complex_trait_analysis/scDRS.py"

rule scdrs_all:
    input:
        expand("results/scDRS/zhang.magma_10kb_1000.74_traits/{trait}/scores.tsv", trait=subset_traits)
    output:
        "temp/scdrs.done"
    shell:
        "echo done > {output}"

rule transfer_connectivities:
    resources:
        mem_mb=200000,
        time="2:00:00"
    input:
        adata_scdrs="data/trajectory_inference/eb_scdrs_preprocessed.adata",
        adata_conn="data/single_cell_objects/eb_pflog1ppfnorm.hvg.umap_embedding.h5ad"
    output:
        scdrs_conn="data/single_cell_objects/eb_scdrs_preprocessed.connectivities_included.h5ad"
    conda:
        "../slurmy/scvi.yml"
    script:
        "../code/complex_trait_analysis/sc_connectivities_liftover.py"

rule scdrs_celltype_analysis:
    resources:
        mem_mb=200000,
        time="2:00:00"
    input:
        adata="data/single_cell_objects/eb_scdrs_preprocessed.connectivities_included.h5ad",
        celltypes="data/fca/eb_cellid_labels.tsv",
        scores="results/scDRS/zhang.magma_10kb_1000.74_traits/{trait}/scores.tsv"
    output:
        celltype_scores="results/scDRS/zhang.magma_10kb_1000.74_traits/{trait}/celltype_scores.tsv"
    conda:
        "../slurmy/scvi.yml"
    script:
        "../code/complex_trait_analysis/scDRS_celltypes.py"

rule scdrs_celltype_analysis_combined:
    input:
        expand("results/scDRS/zhang.magma_10kb_1000.74_traits/{trait}/celltype_scores.tsv", trait=subset_traits)
    output:
        "results/scDRS/zhang.magma_10kb_1000.74_traits/combined_celltype_scores.tsv"
    shell:
        "echo done > {output}"
