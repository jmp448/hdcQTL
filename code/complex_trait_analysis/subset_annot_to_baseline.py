import pandas as pd

annot_loc = snakemake.input['annot']
baseline_loc = snakemake.output['baseline']

annots = pd.read_csv(annot_loc, sep="\t")

# remove the customized added annotations
ben_addons = [s for s in list(annots.columns) if s[-6:]=="MaxCPP"]
annots = annots.drop(columns=ben_addons)

annots.to_csv(baseline_loc, sep="\t", index=False)