library(httr)
library(vroom)
library(tidyverse)

#inputs
dynamic_qtl_bed_loc <- "results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.bed"
bim_loc <- "data/dynamic_qtl_calling/all_trajectories_combined/genotypes_filtered_plink.bim"

#outputs
snplist_loc <- "results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.snplist.txt"
dynamic_qtl_gwas_bed_loc <- "results/static_eqtl_followup/qtl_sets/dynamic-eqtls/opentargets_overlap/dynamic-signif_variant_gene_pairs.gtex_removed.opentargets_overlap.bed"

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

dynamic_eqtls_gtex_removed <- vroom(dynamic_qtl_bed_loc)
bim <- vroom(bim_loc, col_names=c("CHR", "RS", "POS_CM", "POS_BP", "ALLELE_1", "ALLELE_2")) 
opentargets_translator <- filter(bim, RS %in% dynamic_eqtls_gtex_removed$EB_VARIANT_ID) %>%
  mutate(variant_id=paste0(str_replace(CHR, "chr", ""), "_", POS_BP, "_", ALLELE_2, "_", ALLELE_1))

dynamic_eqtls_opentargets_overlap <- bind_rows(map_df(opentargets_translator$variant_id, pull_associations))
dynamic_eqtl_gwas_variants <- opentargets_translator %>%
  filter(variant_id %in% dynamic_eqtls_opentargets_overlap$variant) %>%
  pull(RS)
write(dynamic_eqtl_gwas_variants, snplist_loc)

dynamic_eqtls_gtex_removed_gwas_bed <- dynamic_eqtls_gtex_removed %>%
  filter(EB_VARIANT_ID %in% dynamic_eqtl_gwas_variants) %>%
  write_tsv(dynamic_qtl_gwas_bed_loc)
