library(scran)
library(scater)
library(BiocParallel)

eb_lowpass_sce_dir <- "/project2/gilad/jpopp/ebQTL/data/single_cell_objects/Lowpass.3seqbatches.merged.TEMP.sce"
eb_lowpass_libsizeprocessed_dir <- "/project2/gilad/jpopp/ebQTL/data/single_cell_objects/Lowpass.3seqbatches.merged.TEMP.libsize.normalized.sce"
eb_lowpass_libsizedec_dir <- "/project2/gilad/jpopp/ebQTL/data/single_cell_objects/Lowpass.3seqbatches.merged.TEMP.libsize.decomposed.rds"
eb_lowpass_processed_dir <- "/project2/gilad/jpopp/ebQTL/data/single_cell_objects/Lowpass.3seqbatches.merged.TEMP.scran.normalized.sce"
eb_lowpass_dec_dir <- "/project2/gilad/jpopp/ebQTL/data/single_cell_objects/Lowpass.3seqbatches.merged.TEMP.decomposed.rds"

eb_lowpass <- readRDS(eb_lowpass_sce_dir)
# eb_lowpass_sub <- eb_lowpass[,sample(seq(1, ncol(eb_lowpass)), 10000)]

# just normalize by library size
lib.factors <- librarySizeFactors(eb_lowpass)
eb_lowpass_libsize <- logNormCounts(eb_lowpass, size.factors=lib.factors)
saveRDS(eb_lowpass_libsize, eb_lowpass_libsizeprocessed_dir)

# variance decomposition w library size normalization
system.time(eb_lowpass_libsize_dec <- modelGeneVar(eb_lowpass_libsize,
                                           BPPARAM=MulticoreParam(multicoreWorkers())))
saveRDS(eb_lowpass_libsize_dec, eb_lowpass_libsizedec_dir)


# normalize by deconvolution
## "quick and dirty" clustering
print(paste0("Num workers: ", multicoreWorkers(), "...\n"))
print("Quick cluster...\n")
system.time(clusts <- quickCluster(eb_lowpass,
                                   min.size=50,
                                   BPPARAM=MulticoreParam(multicoreWorkers())))

# compute sum factors by deconvolution
print("\nComputing sum factors...\n")
system.time(eb_lowpass <- computeSumFactors(eb_lowpass,
                                                clusters=clusts,
                                                BPPARAM=MulticoreParam(multicoreWorkers())))
eb_lowpass <- logNormCounts(eb_lowpass)
saveRDS(eb_lowpass, eb_lowpass_processed_dir)

# scran variance decomposition (to get hvf)
print("\nModeling mean-variance...\n")
system.time(eb_lowpass_dec <- modelGeneVar(eb_lowpass,
                                           BPPARAM=MulticoreParam(multicoreWorkers())))
saveRDS(eb_lowpass_dec, eb_lowpass_dec_dir)