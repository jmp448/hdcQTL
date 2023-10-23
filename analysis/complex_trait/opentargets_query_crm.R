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
crm_eqtls <- vroom("results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/crm-signif_variant_gene_pairs.bed")

# Load the GTEx overlapping hits
overlap_crm <- vroom("results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/crm-signif_variant_gene_pairs.full_gtex_overlap.bed",
                              col_names=c("CHR", "START", "STOP", "EB_ENSG", "EB_HGNC", 
                                          "RSID", "EB_CONTEXT", "CHR1", "START1", 
                                          "STOP1", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT"))
novel_crm_eqtls <- anti_join(crm_eqtls, overlap_crm, by=c("EB_ENSG"="EB_ENSG", "EB_VARIANT_ID"="RSID"))

# Load a BIM file for conversion -- dynamic eqtls
bim <- vroom("data/genotypes/yri_maf0.1_all.hg38.bim", col_names=c("CHR", "RS", "CM", "POS", "REF", "ALT")) %>% distinct()
crm_translator <- filter(bim, RS %in% novel_crm_eqtls$EB_VARIANT_ID) %>%
  distinct() %>%
  mutate(variant_id = paste0(str_replace(CHR, "chr", ""), "_", POS, "_", REF, "_", ALT)) %>%
  mutate(variant_id_swapped = paste0(str_replace(CHR, "chr", ""), "_", POS, "_", ALT, "_", REF))

crm_gwas_hits <-  bind_rows(map_df(crm_translator$variant_id_swapped, pull_associations))
crm_gwas_variants <- filter(crm_translator, variant_id_swapped %in% crm_gwas_hits$variant)
crm_gwas <- filter(novel_crm_eqtls, EB_VARIANT_ID %in% crm_gwas_variants$RS)

# Write to a list of rsIDs
crm_snplist <- dplyr::select(crm_gwas, EB_VARIANT_ID) %>%
  distinct() %>%
  write_tsv("results/cellregmap_eqtl_calling/eb_cellid/pseudobulk_tmm/basic/novel_crm_gwas_snplist.txt", col_names=F)

### PARSE RESULTS
# Filter out any CM hits that tag a GTEx eQTL
crm_tagged_overlap <- read.table("results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/crm-signif-tags_variant_gene_pairs.full_gtex_overlap.bed",
                                col.names=c("EB_CHR", "EB_START", "EB_END", "EB_ENSG", "EB_HGNC", "EB_VARIANT_ID", "EB_TAGGED_SNP",
                                            "GTEX_CHR", "GTEX_START", "GTEX_END", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT")) %>%
  as_tibble()
crm_fully_novel <- filter(crm_gwas, !EB_VARIANT_ID %in% crm_tagged_overlap$EB_TAGGED_SNP)
write_tsv(crm_fully_novel, "results/gwas_overlap/crm_fully_novel.tsv")
