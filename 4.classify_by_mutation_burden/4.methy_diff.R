###############################################
# Methylation changes between groups

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

mthy <- readr::read_rds(file.path(tcga_path,"pancan33_meth.rds.gz"))
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list_type <- read.table(file.path(gene_list_path, "checkpoint.type"),header=T)


# methy difference calculate ------------------------------------------------
## t test to get the difference of methylation levels of genes between mutation groups
### data prepare
mutation_burden_class$cancer_types %>% unique() ->cancers_in_mutaion_burden_class
mthy %>%
  dplyr::filter(cancer_types %in% cancers_in_mutaion_burden_class) %>%    # some cancers in TCGA merged serveral cancers into one
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% gene_list$symbol) %>%
  dplyr::select(-gene) %>%
  tidyr::gather(-cancer_types,-symbol,key=barcode,value=methy) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  dplyr::filter(! is.na(methy)) %>%
  tidyr::nest(-symbol) -> genelist_methy_mutaion_class

### function to calculate difference

calculate_fc_pvalue <- function(.x) {
  # pvalue & fdr
  .x %>%
    dplyr::filter(! mutation_status=="Normal") -> .x
  broom::tidy(t.test(methy ~ mutation_status, data = .x)) %>%
    dplyr::mutate(`log2FC(h/l)`=log2(estimate1/estimate2)) %>%
    dplyr::mutate(high_mutation_burden_methy = estimate1,low_mutation_burden_methy = estimate2) %>%
    dplyr::select(high_mutation_burden_methy,low_mutation_burden_methy,`log2FC(h/l)`,p.value) -> df_pvalue
  return(df_pvalue)
}

### get result
genelist_methy_mutaion_class %>%
  dplyr::mutate(diff = purrr::map(data,calculate_fc_pvalue)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::filter(p.value <= 0.05) -> genelist_methy_diff_in_mutaion_groups.sig

genelist_methy_diff_in_mutaion_groups.sig %>%
  readr::write_tsv(file.path(result_path,"4.methy","genelist_methy_diff_in_mutaion_groups.sig"))

## specific gene violin plot 
genelist_methy_mutaion_class %>%
  dplyr::filter(symbol %in% c("BTNL9","ICOSLG","PDCD1LG2","BTN2A1")) %>%
  tidyr::unnest() %>%
  dplyr:::mutate(mutation_status=ifelse(mutation_status=="low_mutation_burden","Low","High")) %>%
  dplyr::rename("Mutation Load" = "mutation_status") %>%
  ggplot(aes(x= `Mutation Load`,y=methy)) +
  geom_violin(trim=F,aes(fill= `Mutation Load`)) +
  geom_boxplot(width=0.1) +
  theme(
    panel.background = element_blank(),
    panel.border =  element_rect(fill = NA)
  ) +
  facet_grid(~symbol)

## t test to get the difference of methylation levels of genes between mutation groups for each cancers
### function
cancers_with_enough_high_mutation <- c("LUAD","SKCM","BLCA","COAD","LUSC","STAD","HNSC","UCEC","BRCA","LIHC","READ","ACC")

fn_cancer_ttest <- function(.x){
  .x %>%
    dplyr::filter(cancer_types.x %in% cancers_with_enough_high_mutation) %>%
    dplyr::group_by(cancer_types.x) %>%
    dplyr::do(
      diff = broom::tidy(t.test(methy ~ mutation_status,data = .))
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest()
}

### calculation
genelist_methy_mutaion_class %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(diff = purrr::map(data,fn_cancer_ttest)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> cancers_methy_diff.result

### plot
cancers_methy_diff.result %>%
  dplyr::mutate(diff = estimate1 - estimate2) %>%
  dplyr::mutate(`-log10P` = -log10(p.value)) %>%
  dplyr::filter(p.value <= 0.05) -> plot_ready

plot_ready %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(sum=sum(diff)) %>%
  dplyr::select(symbol,sum) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::arrange(sum)-> gene_rank
plot_ready %>%
  dplyr::group_by(cancer_types.x) %>%
  dplyr::mutate(sum=sum(diff)) %>%
  dplyr::select(cancer_types.x,sum) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::arrange(sum)-> cancer_rank

plot_ready %>%
  ggplot(aes(x=cancer_types.x,y=symbol)) +
  geom_point(aes(color = diff, size = `-log10P`)) +
  scale_color_gradient2(
    low="blue",
    mid="white",
    high="red",
    midpoint = 0
  ) +
  scale_x_discrete(limit=cancer_rank$cancer_types.x) +
  scale_y_discrete(limit=gene_rank$symbol)



# cluster analysis --------------------------------------------------------

## data preparation
genelist_methy_mutaion_class %>%
  tidyr::unnest() %>%
  # dplyr::filter(symbol %in% gene_list_fc_sig$symbol) %>%
  # dplyr::filter(cancer_types.x %in% c("LUAD","SKCM","BLCA","COAD","LUSC","STAD","HNSC","UCEC")) %>%  # cancers with more samples with high mutation load
  dplyr::filter(cancer_types.x %in% c("LUAD")) %>%  # SKCM,LUAD,"SKCM","BLCA","COAD"
  dplyr::group_by(barcode,symbol) %>%
  dplyr::mutate(methy=mean(methy)) %>%          # use mean value to represent expression of the samples from same patient
  dplyr::ungroup() %>%
  unique() -> gene_list_methy_simplify

gene_list_methy_simplify %>%
  dplyr::select(-cancer_types.x,-cancer_types.y,-mutation_status,-sm_count) %>%
  tidyr::spread(key=barcode, value=methy)  %>%
  as.data.frame() -> gene_list_methy_scale


rownames(gene_list_methy_scale) <- gene_list_methy_scale$symbol
gene_list_methy_scale[,-1] -> gene_list_methy_scale
gene_list_methy_scale %>% t() %>% scale() %>% t() -> gene_list_methy_scale

## side bar preparation

library(ComplexHeatmap)

### Up annotation

gene_list_methy_simplify %>%
  dplyr::select(barcode,mutation_status,sm_count,cancer_types.x) %>%
  unique() %>%
  dplyr::arrange(barcode)-> sample_anno

# RColorBrewer::brewer.pal(8,"Set1") -> cancer_color
# gene_list_methy_simplify$cancer_types %>% unique() -> cancer_tyes
# data.frame(cancer_color=cancer_color,cancer_types=cancer_tyes) -> cancer_anno
# fn_get_cancer_color <- function(i){
#   cancer_anno[i,2] -> x
#   cancer_anno[i,1] -> y
#   print(paste(x, "=", y))
# }

up_anno <- HeatmapAnnotation(df = data.frame(mutation_status=sample_anno$mutation_status,
                                             # sm_count = log2(sample_anno$sm_count),
                                             cancer_types=sample_anno$cancer_types.x),
                             col = list(mutation_status = c("low_mutation_burden" = "#8C8C8C",
                                                            "high_muation_burden" = "#FFC1C1"),
                                        # sm_count = circlize::colorRamp2(c(0,
                                        #                                   max(log2(sample_anno$sm_count))),
                                        #                                 c("white","red")),
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
gene_list_methy_simplify %>%
  dplyr::inner_join(gene_list_type,by="symbol") %>%
  dplyr::select(symbol,type,functionWithImmune) %>%
  unique() %>%
  dplyr::arrange(symbol) -> gene_anno
side_anno <- rowAnnotation(df=data.frame(type=gene_anno$type,function_type=gene_anno$functionWithImmune),
                           col = list(type = c("Ligand" = "#66C2A5",
                                               "Receptor" = "#FC8D62",
                                               "Ligand&Receptor" = "#E78AC3"),
                                      function_type = c("TwoSide" = "#FFFF33",
                                                        "Activate" = "#B3B3B3",
                                                        "Inhibit" = "#A6D854")
                           ),
                           width = unit(0.5, "cm")
)
grid.newpage()
draw(side_anno,1:10)

### determining the optimal number of clusters
# set.seed(123)
# library(NbClust)
# gene_list_methy_scale %>%
#   NbClust(distance = "euclidean",
#           min.nc = 2, max.nc = 10, 
#           method = "complete", index ="all") -> res.nbclust

### get heatmap
library(circlize)
library(dendextend)
row_dend = hclust(dist(gene_list_methy_scale)) # row clustering
plot(row_dend)
col_dend = hclust(dist(t(gene_list_methy_scale))) # column clustering
plot(col_dend)

Heatmap(gene_list_methy_scale, name = "expression",  #km = 5, 
        col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
        top_annotation = up_anno, 
        show_row_names = T, show_column_names = FALSE,
        column_dend_height = unit(30, "mm"),
        clustering_distance_columns = "peason",
        cluster_rows = color_branches(row_dend, k = 5),     # add color on the row tree branches
        cluster_columns = color_branches(col_dend, k = 10)   # add color on the column tree branches
) -> he

he+side_anno

