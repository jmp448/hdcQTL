library(data.table)
library(vroom)

sumstats_loc <- snakemake@input[['sumstats']]
bed_loc <- snakemake@output[['bed']]
  
sumstats <- fread(sumstats_loc)[, .(CHROM, POS_b37, pvalue)]
sumstats <- unique(sumstats)[, HG19_POS := paste0(CHROM, "_", POS_b37)]

bed <- copy(sumstats)[, `:=`(START = POS_b37 - 1)]
setnames(bed, old = "CHROM", new = "#CHR")
setnames(bed, old = "POS_b37", new = "END")
setnames(bed, old = "pvalue", new = "P")
setcolorder(bed, c("#CHR", "START", "END", "HG19_POS", "P"))
fwrite(bed, bed_loc, sep = "\t")
