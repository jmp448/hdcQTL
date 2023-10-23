library(tidyverse)

center.scale <- function(g) {
  # center and scale within a gene's expression vector
  if (sd(g) == 0) {
    g - g
  } else {
    (g-mean(g))/sd(g)
  }
}

cell.line.pca <- function(x, npc=5) {
  if (!is_tibble(x)) {
    x <- as_tibble(x)
  }
  # reshape so columns are (gene, pseudobulk bin, and) individuals
  # filter out the days that weren't fully observed
  x <- x %>% gather(!gene, key="sample", value="counts") %>% 
    mutate(day=str_extract(sample, "[^_]+$")) %>%
    mutate(ind=str_extract(sample, "[^_]+")) %>%
    dplyr::select(!sample) %>% spread(ind, counts) %>%
    drop_na()
  x <- x %>% dplyr::select(!c(gene, day)) %>% apply(1, center.scale) %>% 
    t %>% as_tibble %>% 
    mutate(gene_day=paste(x$gene, x$day, sep="_")) %>%
    relocate(gene_day)
  novargenes <- x %>% dplyr::select(!gene_day) %>% apply(1, var) %>% as_tibble %>%
    mutate(gene_day=x$gene_day) %>% filter(value==0) %>% .$gene_day
  x <- x %>% filter(!gene_day %in% novargenes)
  x.svd <- x %>% dplyr::select(!gene_day) %>% t %>% svd(nu=npc, nv=npc)
  pcgenes <- x.svd %>% .$v %>% as_tibble %>% mutate(gene_day=x$gene_day)
  cell.line.pcs <- x.svd %>% .$u %>% as_tibble %>% 
    `colnames<-`(paste0("PC_", seq(1, npc))) %>%
    mutate(ind=colnames(x)[-c(1)]) %>% relocate(ind)
  var.exp <- x.svd %>% .$d %>% `^`(2) %>% `/`(sum((x.svd$d)^2))
  list("cell.line.pcs"=cell.line.pcs, "pc.genes"=pcgenes, "pve"=var.exp)
}

regular.pca <- function(x, npc=10) {
  if (!is_tibble(x)) {
    x <- as_tibble(x)
  }
  x <- x %>% dplyr::select(!c(gene)) %>% apply(1, center.scale) %>% 
    t %>% as_tibble %>% mutate(gene=x$gene) %>% relocate(gene)
  novargenes <- x %>% dplyr::select(!gene) %>% apply(1, var) %>% as_tibble %>%
    mutate(gene=x$gene) %>% filter(value==0) %>% .$gene
  x <- x %>% filter(!gene %in% novargenes)
  x.svd <- x %>% dplyr::select(!gene) %>% t %>% svd(nu=npc, nv=npc)
  x.svd$u <- x.svd$u %>% as_tibble() %>% `colnames<-`(paste0("PC", seq(1,npc))) %>%
    mutate(sample=colnames(x)[-c(1)]) %>% relocate(sample)
  x.svd$v <- x.svd$v %>% as_tibble() %>% `colnames<-`(paste0("PC", seq(1,npc))) %>%
    mutate(gene=x$gene) %>% relocate(gene)
  x.svd
}

prob.cell.line.pca <- function(x, npc=5, pca_method) {
  if (!is_tibble(x)) {
    x <- as_tibble(x)
  }
  # reshape so columns are (gene, pseudobulk bin, and) individuals
  x <- x %>% gather(!gene, key="sample", value="counts") %>% 
    mutate(day=str_extract(sample, "[^_]+$")) %>%
    mutate(ind=str_extract(sample, "[^_]+")) %>%
    dplyr::select(!sample) %>% spread(ind, counts)
  tailcol <- tail(colnames(x),-2)
  x <- x %>% dplyr::select(!c(gene, day)) %>% apply(1, scale) %>% 
    t %>% as_tibble %>% 
    mutate(gene_day=paste(x$gene, x$day, sep="_")) %>%
    relocate(gene_day)
  colnames(x) <- c('gene_day',tailcol)
  novargenes <- x %>% dplyr::select(!gene_day) %>% apply(1, var) %>% as_tibble %>%
    mutate(gene_day=x$gene_day) %>% filter(value==0) %>% .$gene_day
  x <- x %>% filter(!gene_day %in% novargenes)
  
  x.svd <- x %>% dplyr::select(!gene_day) %>% t
  x.svd <- pcaMethods::pca(x.svd,method=pca_method,nPcs = npc, seed=2021)
  pcgenes <- x.svd@loadings %>% as_tibble %>% mutate(gene_day=x$gene_day)
  prob.cell.line.pcs <- x.svd@scores %>% as_tibble %>% 
    `colnames<-`(paste0("PC_", seq(1, npc))) %>%
    mutate(ind=colnames(x)[-c(1)]) %>% relocate(ind)
  var.exp <- x.svd@sDev %>% `^`(2) %>% `/`(sum((x.svd@sDev)^2))
  list("cell.line.pcs"=prob.cell.line.pcs, "pc.genes"=pcgenes, "pve"=var.exp)
}