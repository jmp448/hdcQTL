library(data.table)
library(vroom)

sumstats_loc <- snakemake@input[['sumstats']]
bed_loc <- snakemake@output[['bed']]
  
sumstats <- fread(sumstats_loc)[, .(CHROM, POS, PVAL)]
sumstats <- unique(sumstats)[, HG19_POS := paste0(CHROM, "_", POS)]

bed <- copy(sumstats)[, `:=`(START = POS - 1)]
setnames(bed, old = "CHROM", new = "#CHR")
setnames(bed, old = "POS", new = "END")
setnames(bed, old = "PVAL", new = "P")
setcolorder(bed, c("#CHR", "START", "END", "HG19_POS", "P"))
fwrite(bed, bed_loc, sep = "\t")
