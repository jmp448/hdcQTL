suppressPackageStartupMessages({
  library("splatter")
  library("scater")
  library("scran")
  library("VariantAnnotation")
  library("ggplot2")
  library("readr")
})
set.seed(42)

# Estimate some parameters from real data
real_counts <- as(t(readMM("/project2/gilad/jpopp/ebQTL/temp/meso_5kgenes_oneind_counts.mtx")), "matrix")
real_params <- splatPopEstimate(counts=real_counts)

# Get details about # cells
sample_summary <- read_tsv("/project2/gilad/jpopp/ebQTL/data/static/highpass_cellid_all/pseudobulk-scran/sample_summary.tsv")
ncells_max <- max(sample_summary$n_cells)
n_ind <- length(unique(sample_summary$individual))

# Load the relevant VCF and GFF
vcf <- readVcf("/project2/gilad/jpopp/ebQTL/data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz")


vcf <- mockVCF(n.samples=6)
gff <- mockGFF(n.genes=5000)

all_params <- setParams(real_params, batchCells=rep(ncells_max, n_ind), 
                        batch.facLoc=0, batch.facScale=1e-10, # negligible batch effect
                        batch.size=1)
sim <- splatPopSimulate(params=all_params, vcf = vcf, gff = gff, sparsify = FALSE)

sim <- logNormCounts(sim)
sim <- runPCA(sim, ncomponents = 5)
plotPCA(sim, colour_by = "Sample")