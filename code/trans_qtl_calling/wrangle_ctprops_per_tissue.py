import pandas as pd

# define paths to data
samples_loc = '../../data/trans_qtl_calling/gtex/samples_per_tissue/samples-Nerve_Tibial.txt'
ctprops_loc = '../../data/trans_qtl_calling/gtex/celltype_proportions/celltype_proportions.v8.xCell.7celltypes.txt'
celltype = 'Nerve_Tibial'
outfile = '../../data/trans_qtl_calling/gtex/celltype_proportions/proportions-Nerve_Tibial.txt'

samples_loc = snakemake.input['tissue_samples']
ctprops_loc = snakemake.input['ctprops']
outfile = snakemake.output['tissue_proportions']

# list samples in this tissue
samples = list(pd.read_csv(samples_loc, header=None, names=['sample'])['sample'])

# subset cell type proportions for this tissue
ctprops = pd.read_csv(ctprops_loc, sep="\t")
ctprops = ctprops[ctprops.columns.intersection(samples).insert(0, 'cell_type')].set_index('cell_type')

# wrangle donor IDs
donor_ids = ctprops.columns.str.extract(r'([^-]+-[^-]+)').iloc[:, 0].tolist()
ctprops.columns = donor_ids

# save 
ctprops.to_csv(outfile, sep="\t")