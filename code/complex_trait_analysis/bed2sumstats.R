library(data.table)

orig_sumstats_loc <- snakemake@input[['orig_sumstats']]
lifted_bed_loc <- snakemake@input[['lifted_bed']]
annot_locs <- snakemake@input[['annots']]
output_sumstats_loc <- snakemake@output[['lifted_sumstats']]

# Load original sumstats for allele info, P vals, etc
orig_sumstats <- fread(orig_sumstats_loc)[CHROM %in% paste0(seq(22)),]
orig_colnames <- colnames(orig_sumstats)
orig_sumstats <- unique(orig_sumstats)[, HG19_POS := paste0(CHROM, "_", POS)]
keeper_cols <- c("HG19_POS", setdiff(orig_colnames, c("CHROM", "POS", "ID")))
orig_sumstats <- orig_sumstats[, ..keeper_cols]

# Load the lifted positions
lifted_bed <- fread(lifted_bed_loc, header=F, col.names = c("CHROM", "START", "POS", "HG19_POS"))
lifted_bed <- unique(lifted_bed)[CHROM %in% paste0(seq(22)), .(CHROM, POS, HG19_POS)]

# Load rsIDs from annotations
rsids <- data.table()
for (annot_loc in annot_locs) {
  print(annot_loc)
  chr_rsids <- fread(annot_loc, select = c("CHR", "BP", "SNP"))
  rsids <- rbind(rsids, chr_rsids)
}

print('arrived')
lifted_sumstats <- lifted_bed[orig_sumstats, on = "HG19_POS", nomatch = 0L]
lifted_sumstats <- lifted_sumstats[, ID := paste0(CHROM, ":", POS, "_", REF, "_", ALT1)]
lifted_sumstats <- lifted_sumstats[,..orig_colnames]
setnames(lifted_sumstats, c("OBS_CT", "REF", "ALT1", "A1_FREQ", "A1"), c("N", "A1", "A2", "MAF", "MA"))
lifted_sumstats <- lifted_sumstats[, CHROM := as.integer(CHROM)]
lifted_sumstats <- lifted_sumstats[, P := as.numeric(P)]
lifted_sumstats <- lifted_sumstats[rsids, on = .(CHROM=CHR, POS=BP), nomatch = 0L]
fwrite(lifted_sumstats, output_sumstats_loc, sep = "\t")
