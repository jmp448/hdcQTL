import pandas as pd

annot_file = snakemake.input['baseline']
gtex_subsampled_eqtls_loc = snakemake.input['gtex_subset']
eb_eqtls_loc = snakemake.input['eb']
annot_out_qtl_loc = snakemake.output['baseline_qtl']
variant_group = snakemake.wildcards['variant_group']
gtex_grouping=snakemake.wildcards['gtex_grouping']

assert gtex_grouping in ['gtex', 'gtexLOO', 'gtexcomb', 'gtexpilotcomb']
# Load data
baseline = pd.read_csv(annot_file, sep="\t")
gtex_eqtls = pd.read_csv(gtex_subsampled_eqtls_loc, sep="\t")
eb_eqtls = pd.read_csv(eb_eqtls_loc, sep="\t")

# Wrangle EB QTL data
if variant_group == "mash-signif":
    eb_eqtls_bed = eb_eqtls[['EB_VARIANT_ID']].drop_duplicates().rename(columns={'EB_VARIANT_ID': 'SNP'})
    eb_eqtls_bed['EB'] = 1
else:
    eb_eqtls_bed = eb_eqtls[['variant_id']].drop_duplicates().rename(columns={'variant_id': 'SNP'})
    eb_eqtls_bed['EB'] = 1

# Wrangle further so each tissue has its own column and QTL significance is encoded 0/1
if gtex_grouping == 'gtexpilotcomb':
    gtex_pilot_tissues = ['Adipose_Subcutaneous', 'Artery_Tibial', 'Heart_Left_Ventricle',
                          'Lung', 'Muscle_Skeletal', 'Nerve_Tibial', 
                          'Skin_Sun_Exposed_Lower_leg', 'Thyroid', 'Whole_Blood']
    gtex_eqtls = gtex_eqtls[['variant_id', 'tissue']]
    gtex_eqtls = gtex_eqtls.loc[gtex_eqtls['tissue'].isin(gtex_pilot_tissues)]
    gtex_eqtls_bed = gtex_eqtls.pivot_table(index='variant_id', columns='tissue', aggfunc=lambda x: 1, fill_value=0).reset_index()
else:
    gtex_eqtls_bed = gtex_eqtls[['variant_id', 'tissue']].pivot_table(index='variant_id', columns='tissue', aggfunc=lambda x: 1, fill_value=0).reset_index()

# Convert to bed format
gtex_eqtls_bed['CHR'] = [int(s.split('_')[0][3:]) for s in gtex_eqtls_bed['variant_id']]
gtex_eqtls_bed['BP'] = [int(s.split('_')[1]) for s in gtex_eqtls_bed['variant_id']]
gtex_eqtls_bed = gtex_eqtls_bed.drop(columns=['variant_id'])

# Combine w baseline
if gtex_grouping in ['gtex', 'gtexLOO']:
    baseline_gtex = baseline.merge(gtex_eqtls_bed, on=['CHR', 'BP'], how='left').fillna(0)
    baseline_gtex_eb = baseline_gtex.merge(eb_eqtls_bed, on=['SNP'], how='left').fillna(0)
    baseline_gtex_eb.iloc[:, -50:] = baseline_gtex_eb.iloc[:, -50:].astype(int)
    # Tweak so that indicator is significance accounting for all but the present tissue ('leave one out')
    if gtex_grouping == 'gtexLOO':
        qtls_combined = baseline_gtex_eb.iloc[:, -50:]
        qtls_combined = qtls_combined.loc[(qtls_combined!=0).any(axis=1)]
        qtls_omitted = 1 - qtls_combined.rename(columns={t: f'Omit_{t}' for t in qtls_combined.columns})
        baseline_gtex_eb = baseline.join(qtls_omitted).fillna(0)
        baseline_gtex_eb.iloc[:, -50:] = baseline_gtex_eb.iloc[:, -50:].astype(int)

elif gtex_grouping in ['gtexcomb', 'gtexpilotcomb']:
    gtexcomb_eqtls_bed = gtex_eqtls_bed[['CHR', 'BP']]
    if gtex_grouping == 'gtexcomb':
        gtexcomb_eqtls_bed['GTEx_eQTL'] = 1
    elif gtex_grouping == 'gtexpilotcomb':
        gtexcomb_eqtls_bed['GTEx_Pilot_eQTL'] = 1
    baseline_gtex = baseline.merge(gtexcomb_eqtls_bed, on=['CHR', 'BP'], how='left').fillna(0)
    baseline_gtex_eb = baseline_gtex.merge(eb_eqtls_bed, on=['SNP'], how='left').fillna(0)
    baseline_gtex_eb.iloc[:, -2:] = baseline_gtex_eb.iloc[:, -2:].astype(int)

# Save annotation file
baseline_gtex_eb.to_csv(annot_out_qtl_loc, sep="\t", index=False)
