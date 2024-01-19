library(httr)
library(vroom)
library(tidyverse)

# Filter a list of associations to genome-wide significant hits
filter_associations <- function(assoc_list) {
  signif_assoc <- lapply(assoc_list, function(assoc) {
    if (assoc$pval < 5e-8) {
      return(assoc)
    }
  })
  signif_assoc <- Filter(Negate(is.null), signif_assoc) # Remove NULL entries
  return(signif_assoc)
}

# Pull PheWAS associations from the OpenTargets GraphQL
pull_associations <- function(variant_id) {
  query_string = "
    query mapGenetoPhenotype($myVariantId: String!){
    pheWAS(variantId: $myVariantId) {
      associations {
        studyId
        pval
      }
    }
  }
  "
  base_url <- "https://api.genetics.opentargets.org/graphql"
  variables <- list("myVariantId" = variant_id)
  post_body <- list(query = query_string, variables = variables)
  r <- POST(url=base_url, body=post_body, encode='json')
  hits <- filter_associations(content(r)$data$pheWAS$associations)
  hits_df <- map_df(hits, ~ tibble::tibble(study = .x$studyId, pval = .x$pval))
  hits_df$variant <- variant_id
  hits_df
}

# Load mash QTLs
mash_eqtls <- vroom("results/static_eqtl_followup/qtl_sets/mash/original/mash-signif_variant_gene_pairs.bed")

# Load a BIM file for conversion
bim <- vroom("data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/all_celltypes_combined/genotypes_filtered_plink.bim", col_names=c("CHR", "RS", "CM", "POS", "ALT", "REF")) %>% distinct()
mash_translator <- filter(bim, RS %in% mash_eqtls$EB_VARIANT_ID) %>%
  distinct() %>%
  mutate(variant_id = paste0(str_replace(CHR, "chr", ""), "_", POS, "_", REF, "_", ALT))

# Write to a list of rsIDs
mash_gwas_hits <-  bind_rows(map_df(mash_translator$variant_id, pull_associations))
mash_gwas_variants <- filter(mash_translator, variant_id %in% mash_gwas_hits$variant)
mash_gwas <- filter(mash_eqtls, EB_VARIANT_ID %in% mash_gwas_variants$RS)

mash_snplist <- dplyr::select(mash_gwas, EB_VARIANT_ID) %>%
  distinct() %>%
  write_tsv("results/cellregmap_eqtl_calling/eb_cellid/pseudobulk_tmm/basic/mash_gwas_snplist.txt", col_names=F)

# # Filter once the above was already run
# mash_snplist <- read_tsv("results/cellregmap_eqtl_calling/eb_cellid/pseudobulk_tmm/basic/novel_mash_gwas_snplist.txt", col_names="EB_VARIANT_ID")
# mash_gwas <- filter(novel_mash_eqtls, EB_VARIANT_ID %in% mash_snplist$EB_VARIANT_ID)
# 
# ### PARSE RESULTS
# # Filter out any hits that tag a GTEx eQTL
# mash_tagged_overlap <- read.table("results/static_eqtl_followup/qtl_sets/mash/original/mash-signif-tags_variant_gene_pairs.all_tissue_overlap.bed",
#                                 col.names=c("EB_CHR", "EB_START", "EB_END", "EB_ENSG", "EB_HGNC", "EB_VARIANT_ID", "EB_TAGGED_SNP",
#                                             "GTEX_CHR", "GTEX_START", "GTEX_END", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT")) %>%
#   as_tibble()
# mash_fully_novel <- filter(mash_gwas, !EB_VARIANT_ID %in% mash_tagged_overlap$EB_TAGGED_SNP)
# write_tsv(mash_fully_novel, "results/gwas_overlap/mash_fully_novel.tsv")
# 
# #### Just looking through results
# mash_novel_hits <- mash_fully_novel %>% unite(test, EB_HGNC, EB_VARIANT_ID, sep="_")
# novel_posteriors <- get_pm(m)[mash_novel_hits$test,]

