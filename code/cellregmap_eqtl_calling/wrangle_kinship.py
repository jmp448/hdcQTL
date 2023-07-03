import pandas as pd

kinship_loc = snakemake.input['kinship']
kinship_id_loc = snakemake.input['kinship_id']

kinship_tsv = snakemake.output['kinship_tsv']

K_id = pd.read_csv(kinship_id_loc, sep="\t", index_col=0, header = None)[1].tolist()
K = pd.read_csv(kinship_loc, sep="\t", names = K_id)
K.index = K_id
assert all(K.columns == K.index) #symmetric matrix, donors x donors

K.to_csv(kinship_tsv, sep = "\t")