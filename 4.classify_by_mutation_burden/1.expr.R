####################################################################
# expression difference between high and low mutation burden samples

library(magrittr)
library(ggplot2)
# data path ---------------------------------------------------------------

immune_path <- "/project/huff/huff/immune_checkpoint"
burden_path <- "/project/huff/huff/data/TCGA"
tcga_path <- file.path(immune_path,"/data/TCGA_data")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint") 

result_path <- "/project/huff/huff/immune_checkpoint/result/mutation_load_class"
# load data ---------------------------------------------------------------

mutation_burden_class <- readr::read_rds(file.path(burden_path,"classfication_of_26_cancers_by_mutation_burden192.rds.gz")) %>%
  tidyr::unnest()

expr_data <- readr::read_rds(file.path(tcga_path,"pancan33_expr.rds.gz"))

gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)

# sample infomation
mutation_burden_class %>% 
  dplyr::filter(mutation_status == "high_muation_burden") %>%
  .$cancer_types %>%
  table() %>%  
  sort(decreasing = TRUE) %>%
  barplot(ylab="Number of samples",
          main="Distribution of samples with high mutation load (>192) in cancers",
          cex.main = 2, cex.axis = 2)

# filter data -------------------------------------------------------------

filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

expr_data %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-expr) %>%
  dplyr::filter(cancer_types %in% mutation_burden_class$cancer_types) -> gene_list_expr

gene_list_expr %>%
  tidyr::unnest() %>%
  tidyr::gather(-cancer_types,-symbol,-entrez_id,key="barcode",value="expr") %>%
  dplyr::filter(! substr(barcode,14,14) == 1) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::select(-cancer_types) %>%
  dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  dplyr::filter(!is.na(expr)) %>%
  tidyr::nest(-symbol) -> gene_list_expr_with_mutation_load

# Differential expression analysis ----------------------------------------------------------------

# function to calculate expression difference -----

calculate_fc_pvalue <- function(.x) {
  
  # pvalue & fdr
  broom::tidy(t.test(expr ~ mutation_status, data = .x)) %>%
    dplyr::mutate(`log2FC(h/l)`=log2(estimate1/estimate2)) %>%
    dplyr::mutate(high_mutation_burden_expr = estimate1,low_mutation_burden_expr = estimate2) %>%
    dplyr::select(high_mutation_burden_expr,low_mutation_burden_expr,`log2FC(h/l)`,p.value) -> df_pvalue
  return(df_pvalue)
}

# calculate the fc and p value for expression change from high to low mutation burden
gene_list_expr_with_mutation_load %>%
  dplyr::mutate(diff = purrr::map(data,calculate_fc_pvalue)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> gene_list_fc

gene_list_fc %>%
  dplyr::filter(p.value <= 0.05) %>%
  dplyr::filter(abs(`log2FC(h/l)`) >= 0.585)

gene_list_fc %>%
  dplyr::mutate(DE=ifelse(p.value<=0.05 & `log2FC(h/l)`>0.585, "Up", "None")) %>%
  dplyr::mutate(DE=ifelse(p.value<=0.05 & `log2FC(h/l)`<(-0.585), "Down", DE)) %>%
  dplyr::mutate(high_mutation_burden_expr=round(high_mutation_burden_expr,2),low_mutation_burden_expr=round(low_mutation_burden_expr,2),`log2FC(h/l)`=round(`log2FC(h/l)`,2),p.value=signif(p.value,2)) %>%
  readr::write_tsv(file.path(result_path,"1.mRNA_DE","DE_of_IMC_between_mutation_groups.tsv"))

## DE information
gene_list_fc %>%
  dplyr::mutate(DE=ifelse(p.value<=0.05 & `log2FC(h/l)`>0.585, "Up", "None")) %>%
  dplyr::mutate(DE=ifelse(p.value<=0.05 & `log2FC(h/l)`<(-0.585), "Down", DE)) %>%
  dplyr::select(DE) %>%
  table() %>%
  sort(decreasing = T) -> de_bar_data

de_bar_data %>% 
  barplot(ylab="Number of genes",
          main="DE Immune Checkpoint between high and low mutation load group",
          cex.main = 1, cex.axis = 2, cex.names = 2, ylim = c(0,50))
text(x=c(1,2,3)-0.2,y=as.numeric(de_bar_data)+1.5, as.character(de_bar_data))

## plot 
gene_list_expr_with_mutation_load %>%
  dplyr::filter(symbol %in% c("CD274","PDCD1","CTLA4","CD80")) %>%
  tidyr::unnest() %>%
  dplyr::mutate(`Log2(expr)`=log2(expr)) %>%
  dplyr:::mutate(mutation_status=ifelse(mutation_status=="low_mutation_burden","Low","High")) %>%
  dplyr::rename("Mutation Load" = "mutation_status") %>%
  ggplot(aes(x= `Mutation Load`,y=`Log2(expr)`)) +
  geom_violin(trim=F,aes(fill= `Mutation Load`)) +
  geom_boxplot(width=0.1) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey")
  ) +
  facet_grid(~symbol) 
  
                  

# cluster analysis --------------------------------------------------------

## data preparation
gene_list_expr_with_mutation_load %>%
  tidyr::unnest() %>%
  dplyr::select(-entrez_id, -sm_count) %>%
  dplyr::group_by(barcode,symbol) %>%
  dplyr::mutate(expr=mean(expr)) %>%          # use mean value to represent expression of the samples from same patient
  dplyr::ungroup() %>%
  unique() -> gene_list_expr_simplify

gene_list_expr_simplify %>%
  dplyr::select(-cancer_types,-mutation_status) %>%
  tidyr::spread(key=barcode, value=expr)  %>%
  as.data.frame() -> gene_list_expr_scale


rownames(gene_list_expr_scale) <- gene_list_expr_scale$symbol
gene_list_expr_scale[,-1] -> gene_list_expr_scale
gene_list_expr_scale %>% scale() -> gene_list_expr_scale

## side bar preparation

library(ComplexHeatmap)

### Up annotation
gene_list_expr_simplify %>%
  dplyr::select(barcode, mutation_status,cancer_types) %>%
  unique() %>%
  tidyr::spread(key=barcode) %>%
  as.data.frame() -> mutation_status_anno

rownames(mutation_status_anno) <- "group"

up_anno_1 <- HeatmapAnnotation(df = mutation_status_anno,
                                  col = list(group=c("T" = "#8C8C8C", "N" = "#FFC1C1")),
                                  width = unit(0.5, "cm"),
                                  name = "mutation freq group")
draw(up_anno_1,1:10)

### cancer annotation
gene_list_expr_simplify %>%
  dplyr::select(barcode, cancer_types) %>%
  unique() %>%
  tidyr::spread(key=barcode,value=cancer_types) %>%
  as.data.frame() -> cancer_anno