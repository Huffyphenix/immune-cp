
tcga_path<-"/project/huff/huff/immune_checkpoint/data"
exp_TP<-readr::read_rds(file.path(tcga_path,"pancancer23_exp_tumorPurity.ras.gz"))

gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)

out_path<-"/project/huff/huff/immune_checkpoint/result_20171025/e_1_tumor_purity"
#-----------------------------------------------------------get exp and TP of genes in gene list 
gene_list %>%
  dplyr::select(symbol) %>%
  dplyr::inner_join(exp_TP,by="symbol") ->genelist_exp_TP
#-----------------------------------------------------------get result
#calculate correlation between gene and purity
genelist_exp_TP %>%
  dplyr::mutate(chrun=rnorm(length(.$purity),sd=1e-6)) %>%
  dplyr::mutate(purity=purity+chrun) %>%
  dplyr::mutate(exp=exp+chrun) %>%
  dplyr::group_by(cancer_types,symbol) %>%
  dplyr::do(
    broom::tidy(
      tryCatch(
        cor.test(.$exp,.$purity,method = c("spearman")),
        error = function(e){1},
        warning = function(e){1})
    )
  ) %>%
  dplyr::mutate(spm_cor=estimate) %>%
  dplyr::select(cancer_types,symbol,spm_cor,p.value,method) %>%
  dplyr::ungroup()->genelist_exp_purity_spm_cor_result

#filter significant result by correlation and pvalue.
genelist_exp_purity_spm_cor_result %>%
  dplyr::filter(abs(spm_cor)>=0.2 & p.value<=0.05) ->sig_gene_cor_tumorPurity

sig_gene_cor_tumorPurity %>%
  readr::write_rds(file.path(out_path,"e_01_gene_sig_cor_tumorPurity.rds.gz"),compress = "gz")

#get rank
fn_gene_rank<-function(pattern){
  pattern %>%
    dplyr::rowwise() %>%
    dplyr::do(
      symbol=.$symbol,
      rank=unlist(.[-1],use.names = F) %>% sum(na.rm=T)
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest() %>%
    dplyr::arrange(rank)
}  
fn_cancer_rank<-function(pattern){
  pattern %>%
    dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(., na.rm = T))) %>%
    tidyr::gather(key = cancer_types, value = rank) %>%
    dplyr::arrange(rank)
}

sig_gene_cor_tumorPurity %>%
  dplyr::mutate(cor_1=ifelse(spm_cor>0,1,-1)) %>%
  dplyr::select(cancer_types,symbol,cor_1) %>%
  tidyr::spread(cancer_types,cor_1) %>%
  fn_gene_rank() %>%
  dplyr::inner_join(gene_list,by="symbol") %>%
  dplyr::select(symbol,rank,col,size) %>%
  dplyr::arrange(col,abs(rank)) ->gene_rank
gene_rank$col %>% as.character() ->gene_rank$col
gene_rank$size %>% as.character() ->gene_rank$size
sig_gene_cor_tumorPurity %>%
  dplyr::mutate(cor_1=ifelse(spm_cor>0,1,-1)) %>%
  dplyr::select(cancer_types,symbol,cor_1) %>%
  tidyr::spread(cancer_types,cor_1) %>%
  fn_cancer_rank() ->cancer_rank

#draw pic
library(ggplot2)
sig_gene_cor_tumorPurity %>%
  dplyr::mutate(dir=ifelse(spm_cor>0,"Pos","Neg")) %>%
  dplyr::mutate(p.value=ifelse(p.value<=1e-10,1e-10,p.value)) %>%
  ggplot(aes(x=cancer_types,y=symbol)) +
  geom_point(aes(size = -log10(p.value), col = spm_cor))+
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = c(0.6,0.4,0.2,0,-0.2,-0.4,0.5,-0.6,-0.8,-1),
    labels = c(0.6,0.4,0.2,0,-0.2,-0.4,0.5,-0.6,-0.8,-1),
    name = "Spearman"
  ) +
  scale_size_continuous(
    name = "-Log10(P)",
    breaks = c(2.5,5,7.5,10),
    labels = c(2.5,5,7.5,10)
  )+
  scale_y_discrete(limits=gene_rank$symbol) +
  scale_x_discrete(limits=cancer_rank$cancer_types)+
  
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
    axis.text.y = element_text(colour = gene_rank$col,face = gene_rank$size),
    axis.text.x = element_text(angle = 45,hjust = 1),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key = element_rect(fill = "white", colour = "black")
  )->p;p
ggsave(
  filename = "fig_01_e_tumorPurity_sig_genes.pdf",
  plot = p,
  device = "pdf",
  width = 8,
  height = 8,
  path = out_path
)
save(list=ls(),file=file.path(out_path,"e_01_tumorPurity.Rdata"))

