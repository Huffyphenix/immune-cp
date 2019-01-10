###########################################################################
# expression profile of ICP genes in cancers ------------------------------

library(magrittr)
library(ggplot2)

# data path config
immune_path <- "/project/huff/huff/immune_checkpoint"
tcga_path <- file.path(immune_path,"/data/TCGA_data")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")
expr_path <-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")

# load data
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))

# calculation of average expression of genes in cancers
# 3 groups by function: activate, inhibit and two sides
fn_average <- function(.data){
  .data %>%
    tidyr::gather(-symbol,-entrez_id,key="barcode",value="expr") %>%
    dplyr::left_join(gene_list,by="symbol") %>%
    dplyr::filter(!is.na(expr)) %>%
    dplyr::group_by(barcode,functionWithImmune) %>%
    dplyr::mutate(average_exp = mean(expr)) %>%
    dplyr::select(barcode,functionWithImmune,average_exp) %>%
    unique() %>%
    dplyr::ungroup()
}

gene_list_expr %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(ICP_mean_expr = purrr::map(filter_expr,fn_average)) %>%
  dplyr::select(-filter_expr) %>%
  dplyr::ungroup() -> ICP_mean_expr_in_cancers

ICP_mean_expr_in_cancers %>%
  tidyr::unnest() %>%
  dplyr::mutate(average_exp = log2(average_exp)) %>%
  ggplot(aes(x=cancer_types,y=average_exp)) +
  geom_violin(aes(fill=functionWithImmune)) +
  facet_wrap(~functionWithImmune) +
  coord_flip() +
  theme_bw() +
  xlab("Cancer Types") +
  ylab("log2 (Expression)") +
  theme(
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(size = 12),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black")
  )
  