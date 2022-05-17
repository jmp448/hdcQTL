library(MatrixEQTL)
library(tidyverse)

snp_file <- snakemake@input[["genotypes"]]
snp_loc_file <- snakemake@input[["snp_locs"]]
expression_file <- snakemake@input[["expression"]]
gene_loc_file <- snakemake@input[["gene_locs"]]
covariate_file <- snakemake@input[["covariates"]]
eqtl_file <- snakemake@output[["eqtls"]]
df_file <- snakemake@output[["df"]]

output.dir <- eqtl_file %>% str_match(".*/") %>% as.character
dir.create(output.dir, showWarnings = FALSE)

useModel = modelLINEAR

# Output file name
output_file_name_cis = eqtl_file
output_file_name_tra = tempfile()

# Include all cis eQTLs, no trans
pvOutputThreshold_cis = 1
pvOutputThreshold_tra = 0

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric()

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
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
if(!is.null(covariate_file)) {
  cvrt$LoadFile(covariate_file)
}

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
