
# Ultilize of FANTOM 5 data and tumor purity data to adjust ICP expression in TCGA ---------------------------------

library(magrittr)
library(org.Hs.eg.db)
library(clusterProfiler)

# data path ---------------------------------------------------------------

immune_path <- "/project/huff/huff/immune_checkpoint"
tcga_path <- file.path(immune_path,"/data/TCGA_data")


# load data ---------------------------------------------------------------

ICP_fantom.gene_exp.cell_line.Immune_cell.combine <- readr::read_rds(file.path(immune_path,"genelist_data","ICP_fantom.gene_exp.cell_line.Immune_cell.combine.rds.gz"))

TCGA_tissue <- readr::read_tsv("/project/huff/huff/data/TCGA/TCGA_cancer_tissue_classification.txt")
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)

filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

gene_list_expr <- readr::read_rds(file.path(tcga_path,"pancan33_expr.rds.gz")) %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-expr) 


