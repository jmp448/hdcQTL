library(data.table)
library(vroom)

sumstats_loc <- snakemake@input[['sumstats']]
bed_loc <- snakemake@output[['bed']]
  
sumstats <- fread(sumstats_loc)[, .(CHR, BP, P)]
sumstats <- unique(sumstats)[, HG19_POS := paste0(CHR, "_", BP)]

bed <- copy(sumstats)[, `:=`(START = BP - 1)]
setnames(bed, old = "CHR", new = "#CHR")
setnames(bed, old = "BP", new = "END")
setcolorder(bed, c("#CHR", "START", "END", "HG19_POS", "P"))
fwrite(bed, bed_loc, sep = "\t")
