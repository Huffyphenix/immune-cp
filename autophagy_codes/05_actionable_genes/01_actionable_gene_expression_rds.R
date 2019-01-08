library(magrittr)
maker_gene_path <- "/extraspace/yye1/share_data/markerGenes"
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")


expr <- readr::read_rds(file.path(tcga_path, "pancan_expr_20160513.rds.gz"))
gene_list <- readr::read_tsv(file.path(maker_gene_path, "actionable.genes.txt")) %>% 
  dplyr::rename(symbol = Gene)

filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

expr %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-expr) -> gene_list_expr

readr::write_rds(gene_list_expr, path = file.path(tcga_path, "actionable_gene_pancan_expr.rds.gz"), compress = "gz")
