library(tidyverse)
library(vroom)
library(RColorBrewer)

gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
eb_medians_loc="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/pseudobulk_tmm_median.tsv"
gtex_loc="results/static_eqtl_followup/gtex/GTEx_gene_median_tpm.gct"
# Get a mapping from HGNC to ENSG
# get a mapping from HGNC to ENSG
pull_gene_type <- function(attr) {
  str_split(attr, "\"")[[1]][6]
}

pull_gene_name <- function(attr) {
  str_split(attr, "\"")[[1]][8]
}

pull_gene_ensg <- function(attr) {
  str_split(attr, "\"")[[1]][2]
}

gencode <- vroom(gtf_loc, col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), skip=5) %>%
  mutate(hgnc=map_chr(attribute, pull_gene_name)) %>%
  mutate(ensg=map_chr(attribute, pull_gene_ensg))

gencode <- vroom(gtf_loc, col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), skip=5) %>%
  filter(seqname %in% paste0("chr", seq(1, 22))) %>%
  filter(feature == "gene") %>%
  mutate(type=map_chr(attribute, pull_gene_type)) %>%
  mutate(hgnc=map_chr(attribute, pull_gene_name)) %>%
  mutate(ensg=map_chr(attribute, pull_gene_ensg)) %>%
  filter(type=="protein_coding")

# 4 HGNC symbols are duplicated here - these have overlapping loci, and will be removed from our analysis
# TBCE, MATR3, HSPA14, GGT1
hgnc.dup <- filter(gencode, hgnc %in% names(table(gencode$hgnc)[table(gencode$hgnc)>1]))
gencode <- gencode %>% filter(!hgnc %in% hgnc.dup$hgnc)
hgnc_ensg_dict <- select(gencode, c(hgnc, ensg))

# Load the EB expression data, and convert gene names
clean_celltype_names <- function(s) {
  paste0("EB_", s)
}
eb_medians <- vroom(eb_medians_loc) %>%
  left_join(hgnc_ensg_dict, by=c("gene"="hgnc")) %>%
  select(-c(gene)) %>%
  drop_na(ensg) %>% # drops 19 genes that didn't map to an ENSG ID
  rename_with(clean_celltype_names, !ensg)

# Load the GTEx data
clean_tissue_names <- function(s) {
  s_cleaned <- str_replace(s, " - ", "_") %>% 
    str_extract("[^(]+") %>%
    str_replace(" ", "_")
  paste0("GTEX_", s_cleaned)
}
gtex_medians <- vroom(gtex_loc, skip=2) %>%
  mutate(ensg=str_extract(Name, "[^.]+"), .keep="unused", .before=1) %>%
  distinct(ensg, .keep_all = T) %>% # drop duplicated ENSG's
  select(-c(Description)) %>%
  rename_with(clean_tissue_names, !ensg)

common_genes <- intersect(eb_medians$ensg, gtex_medians$ensg)

eb_medians_mat <- filter(eb_medians, ensg %in% common_genes) %>%
  arrange(ensg) %>%
  column_to_rownames("ensg") %>% as.matrix()
gtex_medians_mat <- filter(gtex_medians, ensg %in% common_genes) %>%
  arrange(ensg) %>%
  column_to_rownames("ensg") %>% as.matrix()

nonzero_cor <- function(x, y) {
  nonzeros <- which(x + y > 0)
  cor(x[nonzeros], y[nonzeros], method="spearman", use="pairwise.complete.obs")
}

# First, evaluate the correlation within EBs
cor_eb <- apply(eb_medians_mat, 2, function(x) {
  apply(eb_medians_mat, 2, function(y) {
    nonzero_cor(x, y)
  })
})

corrplot(cor_eb,
         order="hclust",
         hclust.method="ward.D2",
         method="color",
         is.corr=F, tl.cex=0.75, cl.cex=0.75)

# Next, evaluate correlation in GTEx
cor_gtex <- apply(gtex_medians_mat, 2, function(x) {
  apply(gtex_medians_mat, 2, function(y) {
    nonzero_cor(x, y)
  })
})

corrplot(cor_gtex,
         order="hclust",
         hclust.method="ward.D2",
         method="color",
         is.corr=F, tl.cex=0.6, cl.cex=0.6)

cor_gtex_eb <- apply(gtex_medians_mat, 2, function(x) {
  apply(eb_medians_mat, 2, function(y) {
    nonzero_cor(x, y)
  })
})

hc_eb <- hclust(as.dist(1 - cor(pairwise_correlations)), method="ward.D2")$order
hc_gtex <- hclust(as.dist(1 - cor(t(pairwise_correlations))), method="ward.D2")$order

heatmap.2(t(cor_gtex_eb),
          density.info=NULL,
          trace="none",
          cexRow = 0.65,
          cexCol = 0.65,
          margins = c(15, 12),
          col=colorRampPalette(brewer.pal(9, "YlOrBr"))(100))

corrplot(cor_gtex_eb[c("EB_Hepatoblasts","EB_Cardiomyocytes", "EB_CNS-neurons"),
                      c("GTEX_Liver", "GTEX_Heart_Left_Ventricle", "GTEX_Brain_Cortex")],
         addCoef.col = "black",
         method="color", is.corr=F, tl.cex=1.2, cl.cex=1.2)

corrplot(cor_gtex_eb[c("EB_Hepatoblasts","EB_Cardiomyocytes", "EB_CNS-neurons"),
                     c("GTEX_Liver", "GTEX_Heart_Left_Ventricle", "GTEX_Brain_Cortex",
                       "GTEX_Cells_Cultured_fibroblasts")],
         addCoef.col = "black",
         method="color", is.corr=F, tl.cex=1.2, cl.cex=1.2)

