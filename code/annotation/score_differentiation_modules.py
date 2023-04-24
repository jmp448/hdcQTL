import numpy as np
import pandas as pd
import scipy
import scanpy as sc
import scdrs
import matplotlib.pyplot as plt
import seaborn as sns

diff_gs_loc = snakemake.input['diff_gs']
fca_loc = snakemake.input['fca']
eb_loc = snakemake.input['eb']

table_prefix = snakemake.params['table_prefix']
fig_prefix = snakemake.params['fig_prefix']

differentiation_gs = scdrs.util.load_gs(diff_gs_loc)
fca = sc.read_h5ad(fca_loc)
eb = sc.read_h5ad(eb_loc)

# create comparable color palettes for the two datasets
fca_types = list(fca.obs.celltype.cat.categories)
fca_palette = sns.color_palette(n_colors=len(fca_types))
eb_types = list(eb.obs.celltype.cat.categories)
eb_palette = sns.color_palette([fca_palette[fca_types.index(t)] for t in eb_types])

for g in list(differentiation_gs.keys()):
    geneset = differentiation_gs[g][0]
    
    fca_scores = scdrs.score_cell(fca, geneset, n_ctrl=100)
    fca_scores['celltype'] = fca.obs['celltype']
    fca_scores.to_csv(table_prefix + "fca/" + g + ".csv", sep="\t")
    
    eb_scores = scdrs.score_cell(eb, geneset, n_ctrl=100)
    eb_scores['celltype'] = eb.obs['celltype']
    eb_scores.to_csv(table_prefix + "eb_subsampled/" + g + ".csv", sep="\t")
    
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(15, 10))
    sns.violinplot(data=fca_scores, x="celltype", y="norm_score", palette=fca_palette, inner=None, ax=ax[0])
    sns.violinplot(data=eb_scores, x="celltype", y="norm_score", palette=eb_palette, inner=None, ax=ax[1])

    ax[0].title.set_text(g + " (FCA)")
    ax[0].tick_params(rotation=90)
    ax[1].title.set_text(g + " (EB)")
    ax[1].tick_params(rotation=90)

    fig.tight_layout(pad=2.0)
    plt.savefig(fig_prefix + g + ".png")
