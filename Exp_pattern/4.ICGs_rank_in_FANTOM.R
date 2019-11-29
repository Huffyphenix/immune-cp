##################################################################
# immune gene expresssion rank in cancer and immune cell
###################################################################
library(magrittr)
library(ggplot2)
# data path ---------------------------------------------------------------
basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
result_path <- file.path(immune_path,"result_20171025","ICP_exp_patthern-byMeanUQ")

# load data ---------------------------------------------------------------

ICP_fantom.gene_exp.cell_line.Immune_cell.combine <-
  readr::read_rds(file.path(immune_path,"genelist_data","FANTOM5","ICP_fantom.gene_exp.cell_line.Immune_cell.raw.exp.rds.gz")) %>%
  dplyr::filter(Group != "Stromal Cell")


TCGA_tissue <- readr::read_tsv(file.path(basic_path,"data/TCGA/TCGA_cancer_tissue_classification.txt"))
TCGA_tissue$Tissues %>% unique()

## get fold change of gene expression between tumor tissue cell line and immune cells ----------
fn_outlier_ratio <- function(cell_line_exp,.symbol){
  ICP_fantom.gene_exp.cell_line.Immune_cell.combine %>%
    dplyr::filter(Group == "Immune Cell") %>%
    dplyr::filter(symbol == .symbol) -> immune_exp
  
  broom::tidy(
    wilcox.test(immune_exp$gene_tpm, cell_line_exp$gene_tpm, alternative = "two.sided") #Comparing the means of two independent groups:Unpaired Two-Samples Wilcoxon Test (non-parametric) 
  ) %>%
    dplyr::mutate(mean_cell_line = mean(cell_line_exp$gene_tpm), 
                  mean_immune_exp = mean(immune_exp$gene_tpm),
                  UQ_cell_line = quantile(cell_line_exp$gene_tpm,0.75),
                  UQ_immune_exp = quantile(immune_exp$gene_tpm,0.75),
                  DQ_cell_line = quantile(cell_line_exp$gene_tpm,0.25),
                  DQ_immune_exp = quantile(immune_exp$gene_tpm,0.25))  -> FC
  
  FC
}

ICP_fantom.gene_exp.cell_line.Immune_cell.combine %>%
  dplyr::filter(! `Characteristics[Tissue]` %in% c("PrimaryCell","blood")) %>%
  tidyr::nest(-entrez_ID,-symbol,.key="cell_line_exp") %>%
  dplyr::mutate(DE_res = purrr::map2(cell_line_exp,symbol,fn_outlier_ratio)) %>%
  dplyr::select(-cell_line_exp) %>%
  tidyr::unnest() -> ICP_DE_FC_and_ratio_between_cellline_immune


# rank the gene -----------------------------------------------------------
ICP_classification <-
  readr::read_tsv(file.path("/home/huff/project/immune_checkpoint/checkpoint/20171021_checkpoint","ICPs_all_info_class-new.tsv")) %>% dplyr::select(symbol,Exp_site,family)

my_theme <-   theme(
  panel.background = element_rect(fill = "white",colour = "black"),
  panel.grid.major=element_line(colour=NA),
  axis.text.y = element_text(size = 8,colour = "black"),
  axis.text.x = element_text(size = 8,colour = "black"),
  # legend.position = "none",
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 12),
  legend.background = element_blank(),
  legend.key = element_rect(fill = "white", colour = "black"),
  plot.title = element_text(size = 20),
  axis.text = element_text(colour = "black"),
  strip.background = element_rect(fill = "white",colour = "black"),
  strip.text = element_text(size = 10),
  text = element_text(color = "black")
)

# UQ expression ---------
ICP_DE_FC_and_ratio_between_cellline_immune %>%
  dplyr::inner_join(ICP_classification, by="symbol") %>%
  dplyr::select(symbol,UQ_cell_line, UQ_immune_exp) %>%
  tidyr::gather(-symbol,key="data_type",value="exp") %>%
  dplyr::mutate(exp=log2(exp+1)) -> ready_draw
ready_draw %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(sum = mean(exp)) %>%
  dplyr::select(symbol,sum) %>%
  unique() %>%
  dplyr::arrange(sum) -> symbol_rank

ready_draw %>%
  ggplot(aes(x=symbol,y=data_type)) +
  geom_tile(aes(fill=exp),size=1,colour="grey",height=0.8) +
  scale_fill_gradient2(low="#00B2EE",high = "#CD2626",mid="white") +
  scale_x_discrete(limits = symbol_rank$symbol) +
  my_theme +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank()
  ) -> p1.1

ICP_classification %>%
  dplyr::mutate(type="Gene family") %>%
  ggplot(aes(y=type,x=symbol)) +
  geom_tile(aes(fill = family),color="grey",size=0.25) +
  scale_x_discrete(limit = symbol_rank$symbol) +
  scale_fill_manual(
    values = c( "#838B8B", "#000000","#0000FF", "#8B2323","#CDAA7D","#8EE5EE","#1C86EE", "#EE3B3B","#EE7600")
  ) +
  my_theme +
  theme(
    panel.grid = element_line(colour = "white"),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(angle = 45,hjust = 1),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    # legend.position = "none",
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.margin=unit(c(0,0,0,-0), "cm")
  ) -> p1.2;p1.2
library(ggpubr)
ggarrange(p1.2,p1.1,
          ncol = 1, nrow = 2,  align = "hv", 
          heights = c(1, 4),
          legend = "left",
          common.legend = TRUE) -> p;p
ggsave(file.path(result_path,"ICGs_rank_in_FANTOM_Sample.UQ.pdf"),plot = p,device = "pdf",height = 8,width = 8)
ggsave(file.path(result_path,"ICGs_rank_in_FANTOM_Sample.UQ.png"),plot = p,device = "png",height = 5,width = 8)

# mean expression ---------
ICP_DE_FC_and_ratio_between_cellline_immune %>%
  dplyr::inner_join(ICP_classification, by="symbol") %>%
  dplyr::select(symbol,mean_cell_line, mean_immune_exp) %>%
  tidyr::gather(-symbol,key="data_type",value="exp") %>%
  dplyr::mutate(exp=log2(exp+1)) -> ready_draw.mean
ready_draw.mean %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(sum = sum(exp)) %>%
  dplyr::select(symbol,sum) %>%
  unique() %>%
  dplyr::arrange(sum) -> symbol_rank.mean

ready_draw.mean %>%
  ggplot(aes(x=symbol,y=data_type)) +
  geom_tile(aes(fill=exp),size=1,colour="grey",height=0.8) +
  scale_fill_gradient2(low="#00B2EE",high = "#CD2626",mid="white") +
  scale_x_discrete(limits = symbol_rank.mean$symbol) +
  my_theme +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank()
  ) -> p2.1
ICP_classification %>%
  dplyr::mutate(type="Gene family") %>%
  ggplot(aes(y=type,x=symbol)) +
  geom_tile(aes(fill = family),color="grey",size=0.25) +
  scale_x_discrete(limit = symbol_rank.mean$symbol) +
  scale_fill_manual(
    values = c( "#838B8B", "#000000","#0000FF", "#8B2323","#CDAA7D","#8EE5EE","#1C86EE", "#EE3B3B","#EE7600")
  ) +
  my_theme +
  theme(
    panel.grid = element_line(colour = "white"),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(angle = 45,hjust = 1),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    # legend.position = "none",
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.margin=unit(c(0,0,0,-0), "cm")
  ) -> p2.2;p2.2
ggarrange(p2.2,p2.1,
          ncol = 1, nrow = 2,  align = "hv", 
          heights = c(1, 4),
          legend = "left",
          common.legend = TRUE) -> p;p
ggsave(file.path(result_path,"ICGs_rank_in_FASNTOM_Sample.mean.pdf"),device = "pdf",height = 5,width = 8)
ggsave(file.path(result_path,"ICGs_rank_in_FANTOM_Sample.mean.png"),device = "png",height = 5,width = 8)

