import pandas as pd
from sklearn.utils import resample
from itertools import chain
import random

random.seed(42)

ld_block_map_loc = snakemake.input['ld_block_map']
resampled_snp_bedfiles = snakemake.output

ld_block_map = pd.read_csv(ld_block_map_loc, sep="\t", names=['ldblock_chr', 'ldblock_start', 'ldblock_end', 'ldblock_liftoverqual', 'qtl_chr', 'qtl_start', 'qtl_end', 'eb_ensg', 'eb_hgnc', 'eb_rs', 'eb_context'])

ld_block_bed = ld_block_map[['ldblock_chr', 'ldblock_start', 'ldblock_end']]
snp_bed = ld_block_map[['qtl_chr', 'qtl_start', 'qtl_end', 'eb_ensg', 'eb_hgnc', 'eb_rs', 'eb_context']]
snp_bed = snp_bed.applymap(lambda x: x.strip() if isinstance(x, str) else x)

ld_blocks = pd.DataFrame(ld_block_bed.apply(lambda row: '-'.join(row.astype(str)), axis=1), columns=['block']).reset_index()
ld_blocks = ld_blocks.groupby('block').agg(list)

for fname in resampled_snp_bedfiles:
    block_sample = resample(ld_blocks.index.to_list())
    ld_blocks_resampled = ld_blocks.loc[block_sample]
    
    snp_sample = list(chain.from_iterable(ld_blocks_resampled['index']))
    snps_resampled = snp_bed.loc[snp_sample]
    
    snps_resampled.to_csv(fname, sep="\t", header=False, index=False)
