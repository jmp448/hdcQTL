library(MatrixEQTL)
library(tidyverse)

snp_file <- snakemake@input[["genotypes"]]
snp_loc_file <- snakemake@input[["snp_locs"]]
expression_file <- snakemake@input[["expression"]]
gene_loc_file <- snakemake@input[["gene_locs"]]
covariate_file <- snakemake@input[["covariates"]]
npcs <- as.numeric(str_replace(snakemake@wildcards[["npcs"]], "pcs", ""))
sample_summary_file <- snakemake@input[["sample_summary"]]
eqtl_file <- snakemake@output[["eqtls"]]
df_file <- snakemake@output[["df"]]

print(npcs)
output.dir <- eqtl_file %>% str_match(".*/") %>% as.character
print(output.dir)
dir.create(output.dir, recursive=T)

useModel = modelLINEAR

# Output file name
output_file_name_cis = eqtl_file
output_file_name_tra = tempfile()

# Include all cis eQTLs, no trans
pvOutputThreshold_cis = 1
pvOutputThreshold_tra = 0

# Distance for local gene-SNP pairs
cisDist = 50000

## Load genotype data
snps = SlicedData$new()
snps$fileDelimiter = "\t"      # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 50000      # read file in slices of 50,000 rows
snps$LoadFile(snp_file)

## Load gene expression data
gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 10000      # read file in slices of 10,000 rows
gene$LoadFile(expression_file)

## Load covariates
all_pcs = read_tsv(covariate_file)
used_pcs = filter(all_pcs, covariate %in% paste0("PC", seq(1, npcs))) %>%
  column_to_rownames("covariate") %>%
  as.matrix
cvrt = SlicedData$new(used_pcs)

# Error covariance matrix
sample_summary = read_tsv(sample_summary_file) %>%
  dplyr::arrange(ind_type) %>%
  filter(!dropped)
errorCovariance = diag(as.numeric(1/sample_summary$n_cells_filtered))

## Run the analysis
snpspos = read.table(snp_loc_file, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_loc_file, header = TRUE, stringsAsFactors = FALSE)

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = TRUE)

unlink(output_file_name_tra)

## Results:
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')

write(me$param$dfFull, file=df_file)
