library(data.table)
library(vroom)

sumstats_loc <- snakemake@input[['sumstats']]
bed_loc <- snakemake@output[['bed']]
chr_filter <- as.numeric(snakemake@params[['chr_filter']])

if (!is.na(chr_filter)) {
  sumstats <- fread(sumstats_loc)[chromosome == chr_filter, .(chromosome, base_pair_location, p_value)]
} else {
  sumstats <- fread(sumstats_loc)[, .(chromosome, base_pair_location, p_value)]
}
sumstats <- unique(sumstats)[, HG19_POS := paste0(chromosome, "_", base_pair_location)]

bed <- copy(sumstats)[, `:=`(START = base_pair_location - 1)]
setnames(bed, old = "chromosome", new = "#CHR")
setnames(bed, old = "base_pair_location", new = "END")
setnames(bed, old = "p_value", new = "P")
setcolorder(bed, c("#CHR", "START", "END", "HG19_POS", "P"))
fwrite(bed, bed_loc, sep = "\t")
