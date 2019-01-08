##############################
# get the correlation between gene expression and immune infiltration
##############################

# result path
result_path <- "/project/huff/huff/immune_checkpoint/result_20171025/e_5_immune_infiltration"

# load data ---------------------------------------------------------------

genelist_exp_inmmune_infiltration <-
  readr::read_rds(file.path("/project/huff/huff/immune_checkpoint/genelist_data","genelist_exp_immuneInfiltration.rds.gz"))

gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)

# function to do correlation ----------------------------------------------

fn_correlation <- function(.x){
  set.seed(12345)
  chrun=rnorm(length(.x$InfiltrationScore),sd=1e-6)
  .x %>%
    dplyr::mutate(InfiltrationScore = InfiltrationScore+chrun,expr=expr+chrun) -> .x
    
  tryCatch(cor.test(.x$expr,.x$InfiltrationScore,method = c("spearman")),
        error = function(e){1},
        warning = function(e){1}) %>%
    broom::tidy()
}

# the correlation of genes in each cancers.
genelist_exp_inmmune_infiltration %>%
  dplyr::mutate(spm = purrr::map(data,fn_correlation)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> genelist_exp_inmmune_infiltration_Cor

genelist_exp_inmmune_infiltration_Cor %>%
  readr::write_tsv(file.path(result_path,"Correlation_between_genelistExp_and_inmmuneInfiltration_each_cancers.tsv"))

# the correlation of gene in all samples.
genelist_exp_inmmune_infiltration %>%
  tidyr::unnest() %>%
  tidyr::nest(-symbol,-entrez_id) %>%
  dplyr::mutate(spm = purrr::map(data,fn_correlation)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> genelist_exp_inmmune_infiltration_Cor.allsamples

genelist_exp_inmmune_infiltration_Cor.allsamples %>%
  readr::write_tsv(file.path(result_path,"Correlation_between_genelistExp_and_inmmuneInfiltration_all_samples.tsv"))


# draw picture ------------------------------------------------------------
# fn get rank

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

# for each cancers --------------------------------------------------------

# get rank
genelist_exp_inmmune_infiltration_Cor %>%
  dplyr::filter(abs(estimate)>=0.2 & p.value<=0.05) %>%
  dplyr::mutate(cor_1=ifelse(estimate>0,1,-1)) %>%
  dplyr::select(cancer_types,symbol,cor_1) %>%
  tidyr::spread(cancer_types,cor_1) %>%
  fn_gene_rank() %>%
  dplyr::inner_join(gene_list,by="symbol") %>%
  dplyr::select(symbol,rank,col,size) %>%
  dplyr::arrange(col,abs(rank)) ->gene_rank
gene_rank$col %>% as.character() ->gene_rank$col
gene_rank$size %>% as.character() ->gene_rank$size

genelist_exp_inmmune_infiltration_Cor %>%
  dplyr::filter(abs(estimate)>=0.2 & p.value<=0.05) %>%
  dplyr::mutate(cor_1=ifelse(estimate>0,1,-1)) %>%
  dplyr::select(cancer_types,symbol,cor_1) %>%
  tidyr::spread(cancer_types,cor_1) %>%
  fn_cancer_rank() ->cancer_rank

#draw pic
library(ggplot2)
genelist_exp_inmmune_infiltration_Cor %>%
  dplyr::filter(abs(estimate)>=0.2 & p.value<=0.05) %>%
  dplyr::mutate(dir=ifelse(estimate>0,"Pos","Neg")) %>%
  dplyr::mutate(p.value=ifelse(p.value<=1e-10,1e-10,p.value)) %>%
  ggplot(aes(x=cancer_types,y=symbol)) +
  geom_point(aes(size = -log10(p.value), col = estimate))+
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = c(0.8,0.6,0.4,0.2,0,-0.2,-0.4),
    labels = c(0.8,0.6,0.4,0.2,0,-0.2,-0.4),
    name = "Spearman"
  ) +
  scale_size_continuous(
    name = "-Log10(P)",
    breaks = c(2.5,5,7.5,10),
    labels = c(2.5,5,7.5,10)
  )+
  scale_y_discrete(limits=gene_rank$symbol) +
  scale_x_discrete(limits=cancer_rank$cancer_types) +
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
  filename = "fig_05_e_immuneIfiltration_sig_genes.pdf",
  plot = p,
  device = "pdf",
  width = 8,
  height = 8,
  path = result_path
)


# For all samples ---------------------------------------------------------
genelist_exp_inmmune_infiltration_Cor.allsamples %>%
  dplyr::filter(abs(estimate)>=0.3 & p.value<=0.05) %>%
  .$symbol -> sig_genes_correlate_with_ImmuneIfiltration

genelist_exp_inmmune_infiltration %>%
  dplyr::filter(symbol %in% sig_genes_correlate_with_ImmuneIfiltration) %>%
  tidyr::unnest() %>%
  dplyr::mutate(expr = log2(expr)) %>%
  ggplot(aes(x=expr,y=InfiltrationScore)) +
  geom_jitter(aes(color = cancer_types),size=1) +
  facet_wrap(~symbol,scales = "free") +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    legend.position = "right",
    legend.background = element_blank()
  ) +
  ylab("Immune Infiltration Score") +
  xlab("log2 (Expression)")
ggsave(
  filename = "fig_05_e_immuneIfiltration_sig_genes_allSamples.pdf",
  device = "pdf",
  width = 12,
  height = 12,
  path = result_path
)
