library(ggplot2)
`%>%` <- magrittr::`%>%`


# Path
out_path <- "/project/huff/huff/immune_checkpoint/result_20171025"
tcga_path <- "/project/huff/huff/immune_checkpoint/data/TCGA_data"
expr_path <-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")
mirna_path <- file.path(out_path, "t_4_mirna")


# load methylation and gene list
mirna_target <- readr::read_rds(file.path(tcga_path, "mirna_gene_target.rds.gz"))

gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)
gene_list$col<-as.character(gene_list$col)
gene_list$size<-as.character(gene_list$size)

mirna_expr <- readr::read_rds(file.path(tcga_path, "pancan33_mirna_expr.rds.gz"))

gene_list %>% 
  dplyr::left_join(mirna_target, by = "symbol") -> gene_list_mirna
gene_list_mirna %>% readr::write_tsv(path = file.path(mirna_path, "02_a_gene_list_mirna.tsv"))

fn_filter_gene_list <- function(.x, gene_list) {
  gene_list_mirna %>%
    dplyr::rename(name = mirna) %>% 
    tidyr::drop_na() %>% 
    dplyr::inner_join(.x, by = "name")
}

mirna_expr %>% 
  dplyr::mutate(mirna = purrr::map(.x = mirna, .f = fn_filter_gene_list)) -> gene_list_mirna_expr 

gene_list_mirna_expr %>% 
  tidyr::unnest()



save.image(file = file.path(mirna_path, ".rda_02_mirna_a_gene_list_target.rda"))
load(file = file.path(mirna_path, ".rda_02_mirna_a_gene_list_target.rda"))
