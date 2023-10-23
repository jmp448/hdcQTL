# Wrangle GO BP gene sets
gmt_file <- "/project2/gilad/jpopp/ebQTL/data/gene_sets/c5.go.bp.v2022.1.Hs.symbols.gmt"
gs_file <- "/project2/gilad/jpopp/ebQTL/data/gene_sets/c5.go.bp.Hs.symbols.gs"

# Redo for 2023 version
gmt_file <- "/project2/gilad/jpopp/ebQTL/data/gene_sets/c5.go.bp.v2023.1.Hs.symbols.gmt"
gs_file <- "/project2/gilad/jpopp/ebQTL/data/gene_sets/c5.go.bp.Hs.2023.symbols.gs"

# Wrangle GO MF gene sets
gmt_file <- "/project2/gilad/jpopp/ebQTL/data/gene_sets/c5.go.mf.v2022.1.Hs.symbols.gmt"
gs_file <- "/project2/gilad/jpopp/ebQTL/data/gene_sets/c5.go.mf.Hs.symbols.gs"

# Read in the gmt file using readLines()
gmt_lines <- readLines(gmt_file)

# Split each line into a list of elements using strsplit()
gmt_list <- lapply(gmt_lines, function(x) unlist(strsplit(x, "\t")))
gmt_list2 <- lapply(gmt_list, function(x) list(x[1], x[2], x[3:length(x)]))

# Convert the list to a tibble
gmt_df <- as_tibble(do.call("rbind", gmt_list2)) %>%
  mutate(across(V1:V2, unlist)) %>%
  dplyr::rename(geneset=V1, link=V2, genes=V3)

# Convert to .gs format for scDRS
gmt_gs <- gmt_df %>%
  dplyr::select(geneset, genes) %>%
  dplyr::rename(TRAIT=geneset) %>%
  mutate(GENESET=sapply(genes, paste, collapse=","), .keep="unused") %>%
  write_tsv(gs_file)