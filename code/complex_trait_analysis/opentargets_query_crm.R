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

# Load CRM QTLs
crm_eqtls <- vroom("results/static_eqtl_followup/qtl_sets/dynamic-eqtls/crm-signif_variant_gene_pairs.bed")

# Load a BIM file for conversion -- dynamic eqtls
bim <- vroom("data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/all_celltypes_combined/genotypes_filtered_plink.bim", col_names=c("CHR", "RS", "CM", "POS", "ALT", "REF")) %>% distinct()
crm_translator <- filter(bim, RS %in% crm_eqtls$EB_VARIANT_ID) %>%
  distinct() %>%
  mutate(variant_id = paste0(str_replace(CHR, "chr", ""), "_", POS, "_", REF, "_", ALT))

crm_gwas_hits <-  bind_rows(map_df(crm_translator$variant_id, pull_associations))
crm_gwas_variants <- filter(crm_translator, variant_id %in% crm_gwas_hits$variant)
crm_gwas <- filter(crm_eqtls, EB_VARIANT_ID %in% crm_gwas_variants$RS)

# Write to a list of rsIDs
crm_snplist <- dplyr::select(crm_gwas, EB_VARIANT_ID) %>%
  distinct() %>%
  write_tsv("results/cellregmap_eqtl_calling/eb_cellid/pseudobulk_tmm/basic/crm_gwas_snplist.txt", col_names=F)