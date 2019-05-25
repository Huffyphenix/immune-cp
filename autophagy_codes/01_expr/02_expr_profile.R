###########################################################################
# expression profile of ICP genes in cancers ------------------------------

library(magrittr)
library(ggplot2)

# data path config
basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
tcga_path <- file.path(immune_path,"/data/TCGA_data")
gene_list_path <- file.path(immune_path,"checkpoint/20171021_checkpoint")
expr_path <- file.path(basic_path,"immune_checkpoint/result_20171025/expr_rds")
out_path <- file.path(basic_path,"immune_checkpoint/result_20171025")

# load data
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol <- as.character(gene_list$symbol)
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))
ICP_expr_pattern <- readr::read_tsv(file.path(out_path,"ICP_exp_patthern","manual_edit_2_ICP_exp_pattern_in_immune_tumor_cell.tsv"))
fn_site_color <- function(.n,.x){
  print(.n)
  if(.x=="Mainly_Tumor"){
    "red"
  }else if(.x=="Mainly_Immune"){
    "Blue"
  }else if(.x=="Both"){
    c("#9A32CD")
  }else{
    "grey"
  }
}
gene_list %>%
  dplyr::inner_join(ICP_expr_pattern,by="symbol") %>%
  dplyr::rename("Exp_site"="Exp site") %>%
  dplyr::mutate(Exp_site=ifelse(is.na(Exp_site),"N",Exp_site)) %>%
  dplyr::mutate(site_col = purrr::map2(symbol,Exp_site,fn_site_color)) %>%
  tidyr::unnest() -> gene_list
# calculation of average expression of genes of samples-----
# 3 groups by function: activate, inhibit and two sides-----
# 样本中在各基因的平均表达量

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
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(sum = quantile(average_exp,0.5)) %>%
  dplyr::select(cancer_types,sum) %>%
  unique() %>%
  dplyr::arrange(sum) -> cancer_rank
ICP_mean_expr_in_cancers %>%
  tidyr::unnest() %>%
  dplyr::mutate(average_exp = log2(average_exp)) %>%
  ggplot(aes(x=cancer_types,y=average_exp)) +
  geom_violin(aes(fill=functionWithImmune),alpha=0.5) +
  facet_wrap(~functionWithImmune) +
  coord_flip() +
  scale_x_discrete(limits= cancer_rank$cancer_types) +
  scale_fill_manual(values = c( "#1C86EE", "#EE3B3B","#EE7600")) +
  theme_bw() +
  xlab("Cancer Types") +
  ylab("log2 (Expression)") +
  theme(
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(size = 12,color = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black"),
    panel.grid = element_line(linetype = "dashed")
  )
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_cancers.pdf"),device = "pdf", width = 6,height = 6)  
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_cancers.png"),device = "png",width = 6,height = 6)  

# 3 groups by function: exp site Mainly_Tumor, Mainly_Immune, both -----
# 样本中在各基因的平均表达量

fn_average.by_expsite <- function(.data){
  .data %>%
    tidyr::gather(-symbol,-entrez_id,key="barcode",value="expr") %>%
    dplyr::left_join(gene_list,by="symbol") %>%
    dplyr::filter(!is.na(expr)) %>%
    dplyr::group_by(barcode,Exp_site) %>%
    dplyr::mutate(average_exp = mean(expr)) %>%
    dplyr::select(barcode,Exp_site,average_exp) %>%
    unique() %>%
    dplyr::ungroup()
}

gene_list_expr %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(ICP_mean_expr = purrr::map(filter_expr,fn_average.by_expsite)) %>%
  dplyr::select(-filter_expr) %>%
  dplyr::ungroup() -> ICP_mean_expr_in_cancers.by_expsite

ICP_mean_expr_in_cancers.by_expsite %>%
  tidyr::unnest() %>%
  dplyr::mutate(average_exp = log2(average_exp)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(sum = quantile(average_exp,0.5)) %>%
  dplyr::select(cancer_types,sum) %>%
  unique() %>%
  dplyr::arrange(sum) -> cancer_rank.by_expsite
ICP_mean_expr_in_cancers.by_expsite %>%
  tidyr::unnest() %>%
  dplyr::filter(Exp_site!="N") %>%
  dplyr::mutate(average_exp = log2(average_exp)) %>%
  ggplot(aes(x=cancer_types,y=average_exp)) +
  geom_violin(aes(fill=Exp_site),alpha=0.5) +
  facet_wrap(~Exp_site) +
  coord_flip() +
  scale_x_discrete(limits= cancer_rank$cancer_types) +
  scale_fill_manual(values = c("#9A32CD", "Blue" ,"red")) +
  theme_bw() +
  xlab("Cancer Types") +
  ylab("log2 (Expression)") +
  theme(
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(size = 12,color = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black"),
    panel.grid = element_line(linetype = "dashed")
  )
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_cancers-by_expsite.pdf"),device = "pdf", width = 6,height = 6)  
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_cancers-by_expsite.png"),device = "png",width = 6,height = 6)  

# calculation of average expression of samples of genes ----
# 基因在各样本中的平均表达量
fn_average_a <- function(.data){
  .data %>%
    tidyr::gather(-symbol,-entrez_id,key="barcode",value="expr") %>%
    dplyr::filter(!is.na(expr)) %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(average_exp = mean(expr)) %>%
    dplyr::select(symbol,average_exp) %>%
    unique() %>%
    dplyr::ungroup()
}

gene_list_expr %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(ICP_mean_expr = purrr::map(filter_expr,fn_average_a)) %>%
  dplyr::select(-filter_expr) %>%
  dplyr::ungroup() -> mean_exp_of_ICP_in_cancers_samples

# heatmap
mean_exp_of_ICP_in_cancers_samples %>%
  tidyr::unnest() %>%
  tidyr::spread(key="symbol",value="average_exp") %>%
  as.data.frame() -> mean_exp_of_ICP_in_cancers_samples.df

mean_exp_of_ICP_in_cancers_samples %>%
  tidyr::unnest() %>%
  dplyr::select(symbol) %>%
  dplyr::inner_join(gene_list,by="symbol") %>%
  dplyr::arrange(symbol) %>%
  dplyr::select(symbol,type,functionWithImmune,Exp_site) %>%
  unique() %>%
  as.data.frame() -> symbol_anno
rownames(symbol_anno) <- symbol_anno[,1]
symbol_anno <- symbol_anno[,-1]
library(ComplexHeatmap)
# gene row annotation
gene_anno <- HeatmapAnnotation(df=symbol_anno,
                               col = list(functionWithImmune=c("Inhibit" = "#EE3B3B",
                                                               "Activate" = "#1C86EE",
                                                               "TwoSide" = "#EE7600"),
                                          Exp_site=c("Mainly_Tumor" = "red",
                                                     "Both" = "#9A32CD",
                                                     "Mainly_Immune" = "Blue",
                                                     "N" = "grey"),
                                          type = c("Receptor" = "black",
                                                   "Ligand" = "red",
                                                   "Ligand&Receptor" = "purple")),
                               width = unit(0.2, "cm"),
                               name = c("Immunity","Type"))

draw(gene_anno,1:20)

rownames(mean_exp_of_ICP_in_cancers_samples.df) <- mean_exp_of_ICP_in_cancers_samples.df[,1]
mean_exp_of_ICP_in_cancers_samples.df <- mean_exp_of_ICP_in_cancers_samples.df[,-1]

mean_exp_of_ICP_in_cancers_samples.df <- as.matrix(mean_exp_of_ICP_in_cancers_samples.df) # complexheatmap need matrix data form

mean_exp_of_ICP_in_cancers_samples.df.scaled <- apply(mean_exp_of_ICP_in_cancers_samples.df,1,scale) %>% t()
colnames(mean_exp_of_ICP_in_cancers_samples.df.scaled) <- colnames(mean_exp_of_ICP_in_cancers_samples.df)

library(circlize)
library(dendextend)

col_dend = hclust(dist(t(mean_exp_of_ICP_in_cancers_samples.df.scaled))) # column clustering
pdf(file.path(out_path,"e_6_exp_profile","ICP_cluster_tree.pdf"),width = 12,height = 4)
plot(col_dend)
dev.off()


he = Heatmap(mean_exp_of_ICP_in_cancers_samples.df.scaled,
             col = colorRamp2(c(-2, 0, 4), c("blue", "white", "red")),
             row_names_gp = gpar(fontsize = 8),
             show_row_names = T, 
             show_column_names = FALSE,
             cluster_rows = F,
             show_row_dend = T, # whether show row clusters.
             top_annotation = gene_anno,
             row_names_side = c("left"),
             cluster_columns = color_branches(col_dend, k = 6),   # add color on the column tree branches
             heatmap_legend_param = list(title = c("Scaled Exp.")))
pdf(file.path(out_path,"e_6_exp_profile","ICP_exp_in_cancers.pdf"),width = 6,height = 4)
he
dev.off()
