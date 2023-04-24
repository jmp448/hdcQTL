import pyarrow.parquet as pq
import pyarrow as pa
import pandas as pd

floc_list = list(snakemake.input)
merged_df_loc = snakemake.output['merged_df']
beta_loc = snakemake.output['beta_df']
se_loc = snakemake.output['se_df']

tables = []

for floc in floc_list:
    table = pq.read_table(floc)
    table = table.append_column('celltype', pa.array([floc.split("/")[5]] * table.num_rows))
    table = table.select(['phenotype_id', 'variant_id', 'celltype', 'slope', 'slope_se'])
    tables.append(table)

# Concatenate
merged_table = pa.concat_tables(tables)

# Convert to pandas data frame
merged_df = merged_table.to_pandas()
merged_df.to_csv(merged_df_loc, sep="\t", index=False)

# Wrangle to matrices of effect sizes and standard errors for mash
# I'll drop any tests that weren't assessed in all cell types to make this more efficient
merged_df['gv'] = merged_df['phenotype_id'].astype(str) + '_' + merged_df['variant_id'].astype(str)

beta_ests = merged_df[['gv', 'celltype', 'slope']].pivot(index='gv', columns='celltype', values='slope')
beta_ests.to_csv(beta_loc, sep="\t", index=True, index_label='gv')

standard_errors = merged_df[['gv', 'celltype', 'slope_se']].pivot(index='gv', columns='celltype', values='slope_se')
standard_errors.to_csv(se_loc, sep="\t", index=True, index_label='gv')
