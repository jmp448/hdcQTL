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

# Load the CM/ HEP dynamic hits
cm_dynamic_eqtls <- vroom("results/dynamic_qtl_calling/eb-cm_15binstrimmed/dynamic-eb-cm-signif_variant_gene_pairs.bed")
hep_dynamic_eqtls <- vroom("results/dynamic_qtl_calling/eb-hep_15binstrimmed/dynamic-eb-hep-signif_variant_gene_pairs.bed")
neur_dynamic_eqtls <- vroom("results/dynamic_qtl_calling/eb-neur_15binstrimmed/dynamic-eb-neur-signif_variant_gene_pairs.bed")

# Load the GTEx overlapping hits
overlap_dynamic_qtls <- vroom("results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/dynamic_variant_gene_pairs.full_gtex_overlap.bed",
                              col_names=c("CHR", "START", "STOP", "EB_ENSG", "EB_HGNC", 
                                          "RSID", "TRAJECTORY", "CHR1", "START1", 
                                          "STOP1", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT"))

novel_cm_dynamic_eqtls <- anti_join(cm_dynamic_eqtls, overlap_dynamic_qtls, by=c("EB_ENSG"="EB_ENSG", "EB_VARIANT_ID"="RSID"))
novel_hep_dynamic_eqtls <- anti_join(hep_dynamic_eqtls, overlap_dynamic_qtls, by=c("EB_ENSG"="EB_ENSG", "EB_VARIANT_ID"="RSID"))
novel_neur_dynamic_eqtls <- anti_join(neur_dynamic_eqtls, overlap_dynamic_qtls, by=c("EB_ENSG"="EB_ENSG", "EB_VARIANT_ID"="RSID"))

# Load a BIM file for conversion
bim <- bind_rows(
  vroom("data/dynamic_qtl_calling/eb-cm_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim", col_names=c("CHR", "RS", "CM", "POS", "REF", "ALT")),
  vroom("data/dynamic_qtl_calling/eb-hep_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim", col_names=c("CHR", "RS", "CM", "POS", "REF", "ALT")),
  vroom("data/dynamic_qtl_calling/eb-neur_15binstrimmed/pseudobulk_tmm/nipals/genotypes_filtered_plink.bim", col_names=c("CHR", "RS", "CM", "POS", "REF", "ALT"))
) %>% distinct()
cm_translator <- filter(bim, RS %in% novel_cm_dynamic_eqtls$EB_VARIANT_ID) %>%
  distinct() %>%
  mutate(variant_id = paste0(str_replace(CHR, "chr", ""), "_", POS, "_", REF, "_", ALT))%>%
  mutate(variant_id_swapped = paste0(str_replace(CHR, "chr", ""), "_", POS, "_", ALT, "_", REF))
hep_translator <- filter(bim, RS %in% novel_hep_dynamic_eqtls$EB_VARIANT_ID) %>%
  distinct() %>%
  mutate(variant_id = paste0(str_replace(CHR, "chr", ""), "_", POS, "_", REF, "_", ALT))%>%
  mutate(variant_id_swapped = paste0(str_replace(CHR, "chr", ""), "_", POS, "_", ALT, "_", REF))
neur_translator <- filter(bim, RS %in% novel_neur_dynamic_eqtls$EB_VARIANT_ID) %>%
  distinct() %>%
  mutate(variant_id = paste0(str_replace(CHR, "chr", ""), "_", POS, "_", REF, "_", ALT))%>%
  mutate(variant_id_swapped = paste0(str_replace(CHR, "chr", ""), "_", POS, "_", ALT, "_", REF))

cm_dynamic_gwas_hits <- bind_rows(map_df(cm_translator$variant_id_swapped, pull_associations))
cm_dynamic_eqtls_gwas_variants <- filter(cm_translator, variant_id_swapped %in% cm_dynamic_gwas_hits$variant)
cm_dynamic_eqtls_gwas <- filter(novel_cm_dynamic_eqtls, EB_VARIANT_ID %in% cm_dynamic_eqtls_gwas_variants$RS)

hep_dynamic_gwas_hits <- bind_rows(map_df(hep_translator$variant_id_swapped, pull_associations))
hep_dynamic_eqtls_gwas_variants <- filter(hep_translator, variant_id_swapped %in% hep_dynamic_gwas_hits$variant)
hep_dynamic_eqtls_gwas <- filter(novel_hep_dynamic_eqtls, EB_VARIANT_ID %in% hep_dynamic_eqtls_gwas_variants$RS)

neur_dynamic_gwas_hits <- bind_rows(map_df(neur_translator$variant_id_swapped, pull_associations))
neur_dynamic_eqtls_gwas_variants <- filter(neur_translator, variant_id_swapped %in% neur_dynamic_gwas_hits$variant)
neur_dynamic_eqtls_gwas <- filter(novel_neur_dynamic_eqtls, EB_VARIANT_ID %in% neur_dynamic_eqtls_gwas_variants$RS)

# Write to a list of rsIDs
cm_snplist <- dplyr::select(cm_dynamic_eqtls_gwas, EB_VARIANT_ID) %>%
  distinct() %>%
  write_tsv("results/dynamic_qtl_calling/eb-cm_15binstrimmed/pseudobulk_tmm/nipals/novel_dynamic_gwas_snplist.txt", col_names=F)

hep_snplist <- dplyr::select(hep_dynamic_eqtls_gwas, EB_VARIANT_ID) %>%
  distinct() %>%
  write_tsv("results/dynamic_qtl_calling/eb-hep_15binstrimmed/pseudobulk_tmm/nipals/novel_dynamic_gwas_snplist.txt", col_names=F)

neur_snplist <- dplyr::select(neur_dynamic_eqtls_gwas, EB_VARIANT_ID) %>%
  distinct() %>%
  write_tsv("results/dynamic_qtl_calling/eb-neur_15binstrimmed/pseudobulk_tmm/nipals/novel_dynamic_gwas_snplist.txt", col_names=F)

### PARSE RESULTS
# Filter out any CM hits that tag a GTEx eQTL
cm_tagged_overlap <- read.table("results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/dynamic-eb-cm-signif-tags_variant_gene_pairs.full_gtex_overlap.bed",
                                col.names=c("EB_CHR", "EB_START", "EB_END", "EB_ENSG", "EB_HGNC", "EB_VARIANT_ID", "EB_TAGGED_SNP",
                                            "GTEX_CHR", "GTEX_START", "GTEX_END", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT")) %>%
  as_tibble()
cm_fully_novel <- filter(cm_dynamic_eqtls_gwas, !EB_VARIANT_ID %in% cm_tagged_overlap$EB_TAGGED_SNP)
write_tsv(cm_fully_novel, "results/gwas_overlap/cm_fully_novel.tsv")

hep_tagged_overlap <- read.table("results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/dynamic-eb-hep-signif-tags_variant_gene_pairs.full_gtex_overlap.bed",
                                col.names=c("EB_CHR", "EB_START", "EB_END", "EB_ENSG", "EB_HGNC", "EB_VARIANT_ID", "EB_TAGGED_SNP",
                                            "GTEX_CHR", "GTEX_START", "GTEX_END", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT")) %>%
  as_tibble()
hep_fully_novel <- filter(hep_dynamic_eqtls_gwas, !EB_VARIANT_ID %in% hep_tagged_overlap$EB_TAGGED_SNP)
write_tsv(hep_fully_novel, "results/gwas_overlap/hep_fully_novel.tsv")

neur_tagged_overlap <- read.table("results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/dynamic-eb-neur-signif-tags_variant_gene_pairs.full_gtex_overlap.bed",
                                 col.names=c("EB_CHR", "EB_START", "EB_END", "EB_ENSG", "EB_HGNC", "EB_VARIANT_ID", "EB_TAGGED_SNP",
                                             "GTEX_CHR", "GTEX_START", "GTEX_END", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT")) %>%
  as_tibble()
neur_fully_novel <- filter(neur_dynamic_eqtls_gwas, !EB_VARIANT_ID %in% neur_tagged_overlap$EB_TAGGED_SNP)
write_tsv(neur_fully_novel, "results/gwas_overlap/neur_fully_novel.tsv")

