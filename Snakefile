import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

include: "rules/sc_preprocessing.py"
include: "rules/fca_annotation.py"
include: "rules/benchmark_static_qtl_calling.py"
include: "rules/benchmark_mash.py"
include: "rules/static_qtl_calling.py"
include: "rules/static_eqtl_followup.py"
include: "rules/trajectory_inference.py"
include: "rules/dynamic_qtl_calling.py"
include: "rules/fast_topics.py"
include: "rules/cellregmap_eqtl_calling.py"
include: "rules/trans_qtl_calling.py"

HTTP = HTTPRemoteProvider()

rule make_conda: 
    input:
        "test.tmp"
    output:
        "test.tmp2"
    conda: "slurmy/tensorqtl.yml"
    shell:
        "echo booyah > {output}"
