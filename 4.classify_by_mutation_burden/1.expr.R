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
gene_list_type <- read.table(file.path(gene_list_path, "checkpoint.type"),header=T)

## sample infomation
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
  dplyr::mutate(T_N = ifelse(substr(barcode,14,14) == 1,"Normal","Tumor")) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::select(-cancer_types) %>%
  dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  dplyr::filter(!is.na(expr)) %>%
  dplyr::mutate(mutation_status = ifelse(T_N=="Normal","Normal",mutation_status)) %>%
  tidyr::nest(-symbol) -> gene_list_expr_with_mutation_load

# Differential expression analysis ----------------------------------------------------------------

# calculate the fc and p value for expression change from high to low mutation burden
## High, Low mutation burden
### function to calculate expression difference -----

calculate_fc_pvalue <- function(.x) {
    # pvalue & fdr
  .x %>%
    dplyr::filter(! mutation_status=="Normal") -> .x
  broom::tidy(t.test(expr ~ mutation_status, data = .x)) %>%
    dplyr::mutate(`log2FC(h/l)`=log2(estimate1/estimate2)) %>%
    dplyr::mutate(high_mutation_burden_expr = estimate1,low_mutation_burden_expr = estimate2) %>%
    dplyr::select(high_mutation_burden_expr,low_mutation_burden_expr,`log2FC(h/l)`,p.value) -> df_pvalue
  return(df_pvalue)
}


gene_list_expr_with_mutation_load %>%
  dplyr::mutate(diff = purrr::map(data,calculate_fc_pvalue)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> gene_list_fc

gene_list_fc %>%
  dplyr::filter(p.value <= 0.05) %>%
  dplyr::filter(abs(`log2FC(h/l)`) >= 0.585) -> gene_list_fc_sig

gene_list_fc %>%
  dplyr::mutate(DE=ifelse(p.value<=0.05 & `log2FC(h/l)`>0.585, "Up", "None")) %>%
  dplyr::mutate(DE=ifelse(p.value<=0.05 & `log2FC(h/l)`<(-0.585), "Down", DE)) %>%
  dplyr::mutate(high_mutation_burden_expr=round(high_mutation_burden_expr,2),low_mutation_burden_expr=round(low_mutation_burden_expr,2),`log2FC(h/l)`=round(`log2FC(h/l)`,2),p.value=signif(p.value,2)) %>%
  readr::write_tsv(file.path(result_path,"1.mRNA_DE","DE_of_IMC_between_mutation_groups.tsv"))

## Normal, High, and Low mutation burden in each cancer types
### function to calculate
cancers_with_enough_high_mutation <- c("LUAD","SKCM","BLCA","COAD","LUSC","STAD","HNSC","UCEC","BRCA","LIHC","READ","ACC")
fn_get_FC <- function(.x){
  .x %>%
    dplyr::filter(! mutation_status %in% "Normal") %>%
    dplyr::filter(cancer_types %in% cancers_with_enough_high_mutation)-> .x
  # get pvalue and FC
  .x %>%
    dplyr::group_by(cancer_types)  %>%
    dplyr::do(broom::tidy(t.test(expr ~ mutation_status, data = .))) %>%
    dplyr::mutate(`log2FC(h/l)`=log2(estimate1/estimate2)) %>%
    dplyr::mutate(high_mutation_burden_expr = estimate1,low_mutation_burden_expr = estimate2) %>%
    dplyr::select(high_mutation_burden_expr,low_mutation_burden_expr,`log2FC(h/l)`,p.value) -> df_pvalue
  return(df_pvalue)
}

gene_list_expr_with_mutation_load %>%
  dplyr::mutate(diff = purrr::map(data,fn_get_FC)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> gene_list_fc_in_cancers

gene_list_fc_in_cancers %>%
  dplyr::filter(p.value <= 0.05) %>%
  dplyr::filter(abs(`log2FC(h/l)`) > 0.585) -> gene_list_fc_in_cancers_sig

## DE information plot
### barplot for overall De results between high and low mutation load groups
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

### buble plot for DE results in each cancers
gene_list_fc_in_cancers_sig %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::filter(! `log2FC(h/l)` == "-Inf") %>%
  dplyr::mutate(sum=sum(`log2FC(h/l)`)) %>%
  dplyr::select(cancer_types,sum) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::arrange(sum) -> cancer_rank
gene_list_fc_in_cancers_sig %>%
  dplyr::group_by(symbol) %>%
  dplyr::filter(! `log2FC(h/l)` == "-Inf") %>%
  dplyr::mutate(sum=sum(`log2FC(h/l)`)) %>%
  dplyr::select(symbol,sum) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::arrange(sum) -> symbol_rank

gene_list_fc_in_cancers_sig %>%
  dplyr::mutate(log10P=-log10(p.value)) %>%
  dplyr::filter(! `log2FC(h/l)` == "-Inf") %>%
  ggplot(aes(x=cancer_types,y=symbol)) +
  geom_point(aes(size=log10P,color=`log2FC(h/l)`)) +
  scale_color_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) +
  scale_y_discrete(limit = symbol_rank$symbol) +
  scale_x_discrete(limit = cancer_rank$cancer_types) 


## specific gene violin plot 
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
  # dplyr::filter(symbol %in% gene_list_fc_sig$symbol) %>%
  dplyr::filter(cancer_types %in% c("LUAD","SKCM","BLCA","COAD","LUSC","STAD","HNSC","UCEC")) %>%  # cancers with more samples with high mutation load
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
gene_list_expr_scale %>% t() %>% scale() %>% t() -> gene_list_expr_scale

## side bar preparation

library(ComplexHeatmap)

### Up annotation

gene_list_expr_simplify %>%
  dplyr::select(barcode,mutation_status,cancer_types) %>%
  unique() -> sample_anno

RColorBrewer::brewer.pal(8,"Set1") -> cancer_color
gene_list_expr_simplify$cancer_types %>% unique() -> cancer_tyes
data.frame(cancer_color=cancer_color,cancer_types=cancer_tyes) -> cancer_anno
fn_get_cancer_color <- function(i){
  cancer_anno[i,2] -> x
  cancer_anno[i,1] -> y
  print(paste(x, "=", y))
}

up_anno <- HeatmapAnnotation(df = data.frame(mutation_status=sample_anno$mutation_status,
                                             cancer_types=sample_anno$cancer_types),
                                  col = list(mutation_status = c("low_mutation_burden" = "#8C8C8C",
                                                               "high_muation_burden" = "#FFC1C1"),
                                             cancer_types = c("BLCA" = "#E41A1C",
                                                            "HNSC" = "#377EB8",
                                                            "LUAD" = "#4DAF4A",
                                                            "SKCM" = "#984EA3",
                                                            "LUSC" = "#FF7F00",
                                                            "UCEC" = "#FFFF33",
                                                            "STAD" = "#A65628",
                                                            "COAD" = "#F781BF")
                                             ),
                                  width = unit(0.5, "cm"))
draw(up_anno, 1:10)


### side annotation
gene_list_expr_simplify %>%
  dplyr::inner_join(gene_list_type,by="symbol") %>%
  dplyr::select(symbol,type,functionWithImmune) %>%
  unique() -> gene_anno
side_anno <- rowAnnotation(df=data.frame(type=gene_anno$type,function_type=gene_anno$functionWithImmune),
                           width = unit(0.5, "cm")
                           )
grid.newpage()
draw(side_anno,1:10)

### determining the optimal number of clusters
set.seed(123)
library(NbClust)
gene_list_expr_scale %>%
  NbClust(distance = "euclidean",
          min.nc = 2, max.nc = 10, 
          method = "complete", index ="all") -> res.nbclust

### get heatmap
library(circlize)
library(dendextend)
row_dend = hclust(dist(gene_list_expr_scale)) # row clustering
col_dend = hclust(dist(t(gene_list_expr_scale))) # column clustering
Heatmap(gene_list_expr_scale, name = "expression",  km = 5, 
        col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
        top_annotation = up_anno, 
        show_row_names = T, show_column_names = FALSE,
        column_dend_height = unit(30, "mm")
        #cluster_rows = color_branches(row_dend, k = 4)     # add color on the row tree branches
        # cluster_columns = color_branches(col_dend, k = 2)   # add color on the column tree branches
        ) -> he

he+side_anno


row_dend = hclust(dist(gene_list_expr_scale)) # row clustering
col_dend = hclust(dist(t(gene_list_expr_scale))) # column clustering
plot(col_dend,hang = -1,cex=.8)      #聚类树状图
re <- rect.hclust(col_dend, k = 10)    #用矩形画出分为3类的区域
re   # 样本分类信息