library(magrittr)
imm_path<-c("/project/huff/huff/immune_checkpoint/data")
gene_cellCount_pbmc_10kImmun.match<-readr::read_rds(file.path(imm_path,"10kImmun_gene_cellCount_pbmc.match.rds.gz"))

#out path
out_path<-"/project/huff/huff/immune_checkpoint/result_20171025/e_4_immunity"

#input gene list
gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)


###extract gene list data
gene_list %>%
  dplyr::select(symbol) %>%
  dplyr::inner_join(gene_cellCount_pbmc_10kImmun.match, by = "symbol") ->gene_list_pbmc.expr

gene_list_pbmc.expr %>%
  readr::write_rds(path = file.path(out_path,"gene_list_expr_cytof_pbmc.rds.gz"),compress = "gz")

#calculate correlation 
gene_list_pbmc.expr %>%
  dplyr::mutate(chrun=rnorm(length(.$count),sd=1e-6)) %>% #add a chrun to avoid ties in the data.
  #ties are equal number of two or more sampels, which will lead inaccuracy of pvalue in spearman correlation analysis.
  dplyr::mutate(count=count+chrun) %>%
  dplyr::mutate(expr=expr+chrun) %>%
  dplyr::group_by(symbol,Immune_cell) %>%
  dplyr::do(
    broom::tidy(
      tryCatch(
        cor.test(.$count, .$expr, methods="spearman"),
        error = function(e){1},
        warning = function(e){1})
    )
  )%>% 
  dplyr::mutate(spm_cor=estimate) %>% 
  dplyr::select(symbol, Immune_cell,spm_cor, p.value) %>%
  dplyr::group_by(symbol) %>% 
  dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
  dplyr::mutate(bfi = p.adjust(p.value, method = "bonferroni")) %>% 
  dplyr::ungroup() ->gene_list_pvalue

#filter with correlation and pvalue
gene_list_pvalue %>%
  dplyr::filter(abs(spm_cor)>=0.25 & p.value<=0.05) ->sig_gene_list_pvalue

sig_gene_list_pvalue %>%
  readr::write_csv(file.path(out_path,"sig_gene_immunity_pvalue.csv"))
#get rank
sig_gene_list_pvalue %>%
  dplyr::mutate(dir=ifelse(spm_cor>0,1,-1)) %>%
  dplyr::select(symbol,Immune_cell,dir) %>%
  tidyr::spread(Immune_cell,dir) %>%
  dplyr::rowwise() %>%
  dplyr::do(
    symbol=.$symbol,
    rank=unlist(.[-1],use.names = F) %>% sum(na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::inner_join(gene_list,by="symbol") %>%
  dplyr::arrange(rank) ->gene_rank

gene_rank$col %>% as.character() ->gene_rank$col
gene_rank$size %>% as.character() ->gene_rank$size

sig_gene_list_pvalue %>%
  dplyr::mutate(dir=ifelse(spm_cor>0,1,-1)) %>%
  dplyr::select(symbol,Immune_cell,dir) %>%
  tidyr::spread(symbol,dir) %>%
  dplyr::rowwise() %>%
  dplyr::do(
    cell_type=.$Immune_cell,
    rank=unlist(.[-1],use.names = F) %>% sum(na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::arrange(rank) ->cell_rank


#-----------------------------------------draw pic
library(ggplot2)


sig_gene_list_pvalue %>% 
  dplyr::mutate(dir=ifelse(spm_cor>0,1,2)) %>%
  ggplot(aes(x=symbol,y=Immune_cell,color=spm_cor)) +
  geom_point(aes(size=-log10(p.value)))+
  scale_x_discrete(limits=gene_rank$symbol)+
  scale_y_discrete(limits=cell_rank$cell_type)+
  scale_size_continuous(
    name ="P-value",
    breaks = c(-log10(0.05),3,5),
    limits = c(-log10(0.05), 10),
    labels = c("0.05", latex2exp::TeX("$10^{-3}$"), latex2exp::TeX("$10^{-5}$"))
  )+
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = c(0.6,0.4,0.2,0,-0.2,-0.4,0.5,-0.6,-0.8,-1),
    labels = c(0.6,0.4,0.2,0,-0.2,-0.4,0.5,-0.6,-0.8,-1),
    name = "Spearman"
  )+
  theme(
    panel.background = element_rect(colour = "black",fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(colour = gene_rank$col,face = gene_rank$size,angle = 45,hjust = 1),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key = element_rect(fill = "white", colour = "black")
  ) ->p;p
ggsave(
  filename = "fig_10_c_immunity_sig_genes.pdf",
  plot = p,
  device = "pdf",
  width = 6,
  height = 5,
  path = out_path
)

save.image(file = file.path(out_path, ".rda_03_c_immunity_gene_expr.rda"))
rm(list=ls())
load(file = file.path(out_path, ".rda_03_c_immunity_gene_expr.rda"))

