import pyarrow.parquet as pq
import pyarrow as pa
import pandas as pd

floc_list = list(snakemake.input)
merged_df_loc = snakemake.output['merged_df']

tables = []

for floc in floc_list:
    table = pd.read_csv(floc, usecols=['phenotype_id', 'variant_id', 'b', 'b_se', 'pval'])
    table['tissue'] = floc.split("/")[5]
    tables.append(table)

merged_table = pd.concat(tables)

merged_table.to_csv(merged_df_loc, sep="\t", index=False)