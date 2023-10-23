import pandas as pd
import glob

rule get_roadmap_promoters:
    input:
        "data/roadmap_epigenomes/full_chromhmm/{epigenome}_15_coreMarks_hg38lift_mnemonics.bed.gz"
    output:
        "data/roadmap_epigenomes/promoters/{epigenome}_promoters.hg38.bed"
    shell:
      """
      zcat {input} | awk '$4 == "1_TssA" || $4 == "2_TssAFlnk" {{print}}' > {output}
      """

rule get_roadmap_enhancers:
    input:
        "data/roadmap_epigenomes/full_chromhmm/{epigenome}_15_coreMarks_hg38lift_mnemonics.bed.gz"
    output:
        "data/roadmap_epigenomes/enhancers/{epigenome}_enhancers.hg38.bed"
    shell:
      """
      zcat {input} | awk '$4 == "6_EnhG" || $4 == "7_Enh" || $4 == "12_EnhBiv" {{print}}' > {output}
      """

rule combine_roadmap_promoters_ipsc:
    input:
        expand("data/roadmap_epigenomes/promoters/{epigenome}_promoters.hg38.bed", epigenome=['E018', 'E019', 'E020', 'E021', 'E022'])
    output:
        "data/roadmap_epigenomes/promoters/promoters.ipsc_combined.hg38.bed"
    # params:
    #     agg_tempfile="temp/agg_tempfile.txt"
    shell:
      """
      module load bedtools
      cat {input} | bedtools sort -i - | bedtools merge -i - > {output}
      """

rule combine_roadmap_enhancers_ipsc:
    input:
        expand("data/roadmap_epigenomes/enhancers/{epigenome}_enhancers.hg38.bed", epigenome=['E018', 'E019', 'E020', 'E021', 'E022'])
    output:
        "data/roadmap_epigenomes/enhancers/enhancers.ipsc_combined.hg38.bed"
    # params:
    #     agg_tempfile="temp/agg_tempfile.txt"
    shell:
      """
      module load bedtools
      cat {input} | bedtools sort -i - | bedtools merge -i - > {output}
      """
      
rule combine_roadmap_promoters_brain:
    input:
        expand("data/roadmap_epigenomes/promoters/{epigenome}_promoters.hg38.bed", epigenome=['E071', 'E074', 'E068', 'E069', 'E072', 'E067', 'E073', 'E070', 'E082', 'E081'])
    output:
        "data/roadmap_epigenomes/promoters/promoters.brain_combined.hg38.bed"
    # params:
    #     agg_tempfile="temp/agg_tempfile.txt"
    shell:
      """
      module load bedtools
      cat {input} | bedtools sort -i - | bedtools merge -i - > {output}
      """

rule combine_roadmap_enhancers_brain:
    input:
        expand("data/roadmap_epigenomes/enhancers/{epigenome}_enhancers.hg38.bed", epigenome=['E071', 'E074', 'E068', 'E069', 'E072', 'E067', 'E073', 'E070', 'E082', 'E081'])
    output:
        "data/roadmap_epigenomes/enhancers/enhancers.brain_combined.hg38.bed"
    shell:
      """
      module load bedtools
      cat {input} | bedtools sort -i - | bedtools merge -i - > {output}
      """

rule combine_roadmap_promoters_heart:
    input:
        expand("data/roadmap_epigenomes/promoters/{epigenome}_promoters.hg38.bed", epigenome=['E065', 'E083', 'E095', 'E104', 'E105'])
    output:
        "data/roadmap_epigenomes/promoters/promoters.heart_combined.hg38.bed"
    shell:
      """
      module load bedtools
      cat {input} | bedtools sort -i - | bedtools merge -i - > {output}
      """

rule combine_roadmap_enhancers_heart:
    input:
        expand("data/roadmap_epigenomes/enhancers/{epigenome}_enhancers.hg38.bed", epigenome=['E065', 'E083', 'E095', 'E104', 'E105'])
    output:
        "data/roadmap_epigenomes/enhancers/enhancers.heart_combined.hg38.bed"
    shell:
      """
      module load bedtools
      cat {input} | bedtools sort -i - | bedtools merge -i - > {output}
      """

rule combine_roadmap_promoters_all:
    input:
        expand("data/roadmap_epigenomes/promoters/{epigenome}_promoters.hg38.bed", epigenome=set(glob_wildcards("data/roadmap_epigenomes/full_chromhmm/{samples}_15_coreMarks_hg38lift_mnemonics.bed.gz").samples))
    output:
        "data/roadmap_epigenomes/promoters/promoters.all_combined.hg38.bed"
    shell:
      """
      module load bedtools
      cat {input} | bedtools sort -i - | bedtools merge -i - > {output}
      """

rule combine_roadmap_enhancers_all:
    input:
        expand("data/roadmap_epigenomes/enhancers/{epigenome}_enhancers.hg38.bed", epigenome=set(glob_wildcards("data/roadmap_epigenomes/full_chromhmm/{samples}_15_coreMarks_hg38lift_mnemonics.bed.gz").samples))
    output:
        "data/roadmap_epigenomes/enhancers/enhancers.all_combined.hg38.bed"
    shell:
      """
      module load bedtools
      cat {input} | bedtools sort -i - | bedtools merge -i - > {output}
      """
