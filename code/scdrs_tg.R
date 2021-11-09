library(tidyverse)
library(scran)
library(scater)

# load gene set
magma_genesets <- read_tsv("/project2/gilad/jpopp/ebQTL/data/scDRS/gs_file/magma_10kb_1000.74_traits.gs")
tg_geneset <- magma_genesets %>%
  filter(TRAIT=="UKB_460K.biochemistry_Triglycerides") %>%
  .$GENESET %>%
  str_split(",") %>%
  .[[1]]

# load EB data - iPSC to endoderm
eb_endo <- readRDS("/project2/gilad/jpopp/ebQTL/data/single_cell_objects/Lowpass.3seqbatches.merged.endoderm.raw.sce")
lib.factors <- librarySizeFactors(eb_endo)
eb_endo_libsize <- logNormCounts(eb_endo, size_factors=lib.factors)

# model variance
eb_endo_libsize_dec <- modelGeneVar(eb_endo_libsize,
                                    BPPARAM=MulticoreParam(multicoreWorkers()-2))

# visualize the iPSC -> endoderm transition
hvg.endo <- getTopHVGs(eb_endo_libsize_dec, fdr.threshold=0.05)
eb_endo_pcs <- runPCA(eb_endo_libsize, subset_row = hvg.endo)
pcviz <- as_tibble(reducedDim(eb_endo_pcs)) %>%
  `colnames<-`(paste0("PC", seq(1, 50))) %>%
  mutate(ALB=logcounts(eb_endo_libsize)["ALB",]) %>%
  mutate(POU5F1=logcounts(eb_endo_libsize)["POU5F1",]) %>%
  mutate(size_factor=sizeFactors(eb_endo_libsize)) %>%
  mutate(leiden=colData(eb_endo_libsize)$leiden)

ggplot(pcviz, aes(x=PC1, y=PC2, color=leiden)) +
  geom_point(alpha=0.5, size=0.2) # doesn't look great, but not going to worry about it

# get a control set of genes
gene_data <- tibble("gene"=rownames(eb_endo_libsize_dec),
                    "mean"=unname(eb_endo_libsize_dec$mean),
                    "total.variance"=unname(eb_endo_libsize_dec$total),
                    "tech.variance"=unname(eb_endo_libsize_dec$tech)) %>%
  mutate(mean_bin=ntile(mean, 15)) %>%
  mutate(var_bin=ntile(total.variance, 15)) %>%
  unite("mean_var_bin", mean_bin, var_bin) %>%
  mutate(is.gwas.gene=gene %in% tg_geneset)

n_per_bin <- gene_data %>%
  group_by(mean_var_bin) %>%
  summarize(n_genes=sum(is.gwas.gene))

control_genes <- gene_data %>%
  filter(!is.gwas.gene) %>%
  select(mean_var_bin, gene) %>%
  group_by(mean_var_bin) %>%
  nest %>%
  ungroup %>%
  left_join(n_per_bin, by="mean_var_bin") %>%
  mutate(samp=map2(data, n_genes, sample_n)) %>%
  select(-data) %>%
  unnest(samp)

# compute 