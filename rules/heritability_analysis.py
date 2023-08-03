import glob

baseline_cats=set(glob_wildcards("data/sldsc/baselineLD_v2.2_bedfiles/{categories}.bed").categories)

rule prune_ben_baseline:
    resources:
      mem_mb=10000,
      time="10:00"
    input:
        annot="data/ldsc/ben_hg38_baselineld_v2.2/baselineLD.{chr}.annot.gz"
    output:
        baseline="data/ldsc/baselineLD_v2.2_annot/baselineLD.{chr}.annot.gz"
    conda: 
        "../slurmy/tensorqtl.yml"
    script:
        "../code/complex_trait_analysis/subset_annot_to_baseline.py"

rule estimate_ld_scores:
    resources:
        mem_mb=30000,
        time="1:00:00"
    input:
        annot_file="data/ldsc/{annot_group}/baselineLD.{chr}.annot.gz",
        print_snps="data/ldsc/"
    output:
        expand("data/ldsc/{{annot_group}}/baselineLD.{{chr}}.{ext}", ext=["l2.ldscore", "l2.M_5_50"])
    params:
        plink_prefix="data/ldsc/hg38_plink_files/1000G.EUR.hg38"
    conda: 
        "../code/ldsc/environment.yml"
    shell:
        """
        python code/ldsc/ldsc.py \
            --l2 \
            --bfile {params.plink_prefix}.{wildcards.chr} \
            --ld-wind-cm 1 \
            --annot {input.annot_file} \
            --out data/ldsc/{wildcards.annot_group}/baselineLD.{chr} \
            --print-snps temp/ldsc_wtf.txt
        """
