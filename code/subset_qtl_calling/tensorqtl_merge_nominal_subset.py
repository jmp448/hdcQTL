import pyarrow.parquet as pq
import pyarrow as pa
import pandas as pd
import re

floc_list = list(snakemake.input)
merged_df_loc = snakemake.output['merged_df']

tables = []

for floc in floc_list:
    table = pq.read_table(floc, columns=['phenotype_id', 'variant_id', 'slope', 'slope_se', 'pval_nominal'])
    table = table.append_column('tissue', pa.array([re.split(r'[./]', floc)[4]] * table.num_rows))
    tables.append(table)

# Concatenate
merged_table = pa.concat_tables(tables)

# Convert to pandas data frame
merged_df = merged_table.to_pandas()
merged_df.to_csv(merged_df_loc, sep="\t", index=False)