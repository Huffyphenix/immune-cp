###################################################################
# immune gene expresssion pattern in cancer and immune cell
###################################################################
library(magrittr)
library(ggplot2)
library(ggrepel)
library(latex2exp)
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

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}
check_outlier <- function(v, coef=1.5){
  quantiles <- quantile(v,probs=c(0.25,0.75))
  IQR <- quantiles[2]-quantiles[1]
  res <- v < (quantiles[1]-coef*IQR)|v > (quantiles[2]+coef*IQR)
  return(res)
}
max_3 <- function(.x){
  top4 <- sort(.x,decreasing = T)[2]
  res <- .x > top4 
  return(res)
}

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
                  mid_cell_line = quantile(cell_line_exp$gene_tpm,0.75),
                  mid_immune_exp = quantile(immune_exp$gene_tpm,0.75),
                  down_cell_line = quantile(cell_line_exp$gene_tpm,0.25),
                  down_immune_exp = quantile(immune_exp$gene_tpm,0.25)) %>%
    dplyr::mutate(`log2FC(I/T).mean` = log2((mean_immune_exp+1)/(mean_cell_line+1)))%>%
    dplyr::mutate(`log2FC(I/T).mid` = log2((mid_immune_exp+1)/(mid_cell_line+1)))%>%
    dplyr::mutate(`log2FC(I/T).down` = log2((down_immune_exp+1)/(down_cell_line+1))) -> FC
  
  
  FC
}
ICP_fantom.gene_exp.cell_line.Immune_cell.combine %>%
  dplyr::filter(! `Characteristics[Tissue]` %in% c("PrimaryCell","blood")) %>%
  tidyr::nest(-entrez_ID,-symbol,.key="cell_line_exp") %>%
  dplyr::mutate(DE_res = purrr::map2(cell_line_exp,symbol,fn_outlier_ratio)) %>%
  dplyr::select(-cell_line_exp) %>%
  tidyr::unnest() -> ICP_DE_FC_and_ratio_between_cellline_immune

my_theme <-   theme(
  panel.background = element_rect(fill = "white",colour = "black"),
  panel.grid.major=element_line(colour=NA),
  axis.text.y = element_text(size = 10,colour = "black"),
  axis.text.x = element_text(size = 10,colour = "black"),
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

# classification of IC expression pattern ---------------------------------

broom::tidy(
  cor.test(ICP_DE_FC_and_ratio_between_cellline_immune$`log2FC(I/T).mid`,
           ICP_DE_FC_and_ratio_between_cellline_immune$`log2FC(I/T).mean`,method = "spearman")) %>%
  dplyr::mutate(label = paste("r = ",signif(estimate,2),"\np = ",signif(p.value,2),sep="")) -> cor_label

ICP_DE_FC_and_ratio_between_cellline_immune %>%
  dplyr::mutate(Exp_site = ifelse(`log2FC(I/T).mean` >=1 & `log2FC(I/T).mid` >=2, "Immune cell dominate","Immune and tumor cell almost")) %>%
  dplyr::mutate(Exp_site = ifelse(`log2FC(I/T).mean` <=(-1) & `log2FC(I/T).mid` <=(-2), "Tumor cell dominate",Exp_site)) %>%
  ggplot(aes(x=`log2FC(I/T).mean`,y=`log2FC(I/T).mid`)) +
  geom_jitter(aes(color = Exp_site)) +
  geom_smooth(method = "lm") +
  # geom_text_repel(aes(x=`log2FC(I/T).mean`,y=`log2FC(I/T).mid`,label=symbol)) +
  geom_label(x=4,y=10,aes(label=label),data = cor_label) +
  geom_hline(yintercept = c(-2,2),linetype = 2) +
  geom_vline(xintercept = c(-1,1),linetype = 2) +
  labs(x=TeX("log_2 ($\\frac{mean(Immune)+1}{\\mean(Tumor)+1})"),
       y=TeX("log_2 ($\\frac{UQ(Immune)+1}{\\UQ(Tumor)+1})"),
       title = "Classification of ICPs' expression pattern") +
  scale_color_manual(values = c("#CD950C", "#66CD00", "#EE2C2C"),
                     name = "ICPs expression pattern") +
  my_theme +
  theme(
    plot.title = element_text(size=15)
  )
  
ggsave(file.path(result_path,"classify_ICP_exp_pattern.pdf"),device = "pdf",height = 4, width = 8)
ggsave(file.path(result_path,"classify_ICP_exp_pattern.png"),device = "png",height = 4, width = 8)




ICP_DE_FC_and_ratio_between_cellline_immune %>%
  dplyr::mutate(Exp_site = ifelse(`log2FC(I/T).mean` >=1 & `log2FC(I/T).mid` >=2, "Immune cell dominate","Immune and tumor cell almost")) %>%
  dplyr::mutate(Exp_site = ifelse(`log2FC(I/T).mean` <=(-1) & `log2FC(I/T).mid` <=(-2), "Tumor cell dominate",Exp_site)) -> ICP_Exp_site_by_DE_Fc_and_ratio_between_cellline_immune

ICP_Exp_site_by_DE_Fc_and_ratio_between_cellline_immune %>%
  readr::write_tsv(file.path(result_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv"))

# draw picture ------------------------------------------------------------
strip_color <- data.frame(Exp_site = unique(ICP_Exp_site_by_DE_Fc_and_ratio_between_cellline_immune$exp_site),
                          site_cplor = c("green", "yellow", "pink"),
                          rank = c(3,2,1))
ICP_fantom.gene_exp.cell_line.Immune_cell.combine  %>%
  dplyr::filter(! `Characteristics[Tissue]` %in% c("blood")) %>%
  # dplyr::mutate(Group = ifelse(Group=="Stromal Cell","Stromal",Group)) %>%
  dplyr::mutate(Group = ifelse(Group=="Immune Cell","Immune",Group)) %>%
  dplyr::mutate(Group = ifelse(Group=="Tumor Cell","Tumor",Group)) %>%
  dplyr::mutate(gene_tpm = log2(gene_tpm+1)) %>%
  dplyr::inner_join(ICP_Exp_site_by_DE_Fc_and_ratio_between_cellline_immune,by=c("symbol","entrez_ID")) %>%
  dplyr::inner_join(strip_color,by="Exp_site") %>%
  dplyr::filter(!is.na(`Characteristics[Tissue]`)) -> ready_for_draw

ready_for_draw %>%
  dplyr::arrange(rank,`log2FC(I/T).mean`) %>%
  .$symbol -> sort.symbol

ready_for_draw <- within(ready_for_draw, symbol <- factor(symbol, levels = unique(sort.symbol))) # change symbol's factor level to rank the facets on plot
with(ready_for_draw, levels(symbol))


ready_for_draw %>%
  dplyr::select(symbol,Exp_site) %>%
  # dplyr::filter(symbol %in% c("CD276","CTLA4")) %>%
  unique() -> color_bac
color_bac$Group <- color_bac$gene_tpm <- 1

library(ggbeeswarm)
ggplot(ready_for_draw,
       aes(x = Group,y = gene_tpm)) +
  # geom_violin() +
  geom_quasirandom(size=0.2) +
  geom_rect(data=color_bac,aes(fill = Exp_site),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.1) +
  facet_wrap(~symbol,scale = "free_y") +
  # ggpubr::stat_compare_means(method = "wilcox.test",label = "p.format") +
  # geom_text(aes(label = outlier), na.rm = TRUE, nudge_x = -0.25, nudge_y = 0.25,size=3) +
  scale_fill_manual(
    name = "ICPs expression pattern",
    values = c("yellow",  "green","pink"),
    # values = c("#008B00", "#00EE00", "#CD8500", "#FF4500"),
    breaks = c("Immune and tumor cell almost", "Immune cell dominate","Tumor cell dominate")
  ) +
  # theme_bw() +
  ylab(TeX("log_2 (TPM+1)")) +
  theme(
    panel.background = element_rect(fill = "white",colour = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_blank(),
    # legend.position = "none",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    legend.position = "top",
    plot.title = element_text(size = 20),
    axis.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "white",colour = "black"),
    strip.text = element_text(size = 10),
    text = element_text(color = "black")
  ) 

ggsave(file.path(result_path,"pattern_info","FANTOM5.ICP_exp_in_tumor(no-blood-tissue)_immune-TCGA.tissue.png"),device = "png",height = 10,width = 20)  
ggsave(file.path(result_path,"pattern_info","FANTOM5.ICP_exp_in_tumor(no-blood-tissue)_immune-TCGA.tissue.pdf"),device = "pdf",height = 10,width = 20)  


# detailed tissue exp pattern ------------------------------------------------------

fn_exp_detailed_pantern_classify <- function(.n,.x){
  print(.n)
  .x %>% 
    dplyr::filter(Group == "Stromal Cell") %>%
    .$gene_tpm -> Stromal_exp
  .x %>% 
    dplyr::filter(Group == "Immune Cell") %>%
    .$gene_tpm -> Immune_exp
  Immune_exp <- mean(Immune_exp)
  
  .x %>% 
    dplyr::filter(Group == "Tumor Cell") %>%
    dplyr::group_by(`Characteristics[Tissue]`) %>%
    dplyr::mutate(gene_tpm = mean(gene_tpm)) %>%
    dplyr::select(`Characteristics[Tissue]`,gene_tpm) %>%
    unique() %>%
    # dplyr::mutate(Immune_exp = ifelse(Immune_exp==0,0.01,Immune_exp)) %>%
    # dplyr::mutate(gene_tpm = ifelse(gene_tpm==0,0.01,gene_tpm)) %>%
    dplyr::mutate(log2FC = log2((Immune_exp+1)/(gene_tpm+1))) %>%
    dplyr::mutate(Exp_patern = purrr::map(log2FC,.f=function(log2fc){
      if(log2fc >= 1){
        "Immune Higher"
      } else if(log2fc <= -1){
        "Tumor Higher"
      } else {
        "Little diff"
      }
    })) %>%
    tidyr::unnest()
}

ICP_fantom.gene_exp.cell_line.Immune_cell.combine %>%
  dplyr::filter(! `Characteristics[Tissue]` %in% c("blood")) %>%
  dplyr::group_by(symbol) %>%
  dplyr::select(-entrez_ID) %>%
  tidyr::nest() %>%
  dplyr::mutate(Exp_pattern = purrr::map2(symbol,data,fn_exp_detailed_pantern_classify)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ICP_exp_pattern_in_immune_tumor_cell.detailed_tissue

ICP_exp_pattern_in_immune_tumor_cell.detailed_tissue %>%
  readr::write_tsv(file.path(result_path,"pattern_info","ICP_exp_pattern_in_immune_tumor(no-blood)_cell.detailed_tissues.tsv"))
# plot 
# classification of tcga cancers

# ICP_expr_pattern <- readr::read_tsv(file.path(result_path,"manual_edit_2_ICP_exp_pattern_in_immune_tumor_cell.tsv"))
fn_site_color <- function(.x){
  # print(.n)
 if(.x=="Tumor cell dominate"){
    "red"
  }else if(.x=="Immune cell dominate"){
    c("darkgreen")
  }else if(.x=="Immune and tumor cell almost"){
    c("darkorange")
  }else{
    "grey"
  }
}
ICP_Exp_site_by_DE_Fc_and_ratio_between_cellline_immune %>%
  dplyr::filter(!is.na(`Exp_site`)) %>%
  dplyr::mutate(site_col = purrr::map(`Exp_site`,fn_site_color)) %>%
  tidyr::unnest() %>%
  dplyr::inner_join(strip_color,by="Exp_site") %>%
  dplyr::arrange(rank,`log2FC(I/T).mean`) -> gene_rank

ICP_exp_pattern_in_immune_tumor_cell.detailed_tissue %>%
  dplyr::filter(`Characteristics[Tissue]` %in% c(TCGA_tissue$Tissues %>% unique())) %>%
  ggplot(aes(x=`Characteristics[Tissue]`,y=symbol)) +
  geom_tile(aes(fill=log2FC),color = "grey",size=0.5,width = 0.9) +
  scale_fill_gradient2(
    name = TeX("log_2 ($\\frac{mean(Immune+1)}{\\mean(Tumor+1)})"),
    low = "blue",
    high = "red",
    mid="white",
    midpoint = 0 
  ) +
  scale_y_discrete(limits = gene_rank$symbol) +
  # guides(fill = guide_colorbar(title.position = "left")) +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.text.y = element_text(size = 10,colour = gene_rank$site_col),
    axis.text.x = element_text(vjust = 0.5, hjust = 1, size = 12, angle = 90),
    axis.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(size = 20),
    axis.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "white",colour = "black"),
    strip.text = element_text(size = 12)
  ) -> p1;p1

ICP_exp_pattern_in_immune_tumor_cell.detailed_tissue %>%
  dplyr::rename("Tissues"="Characteristics[Tissue]") %>%
  dplyr::inner_join(TCGA_tissue,by="Tissues") %>%
  dplyr::arrange(Tissues) %>%
  dplyr::select(TCGA_Cancer) %>%
  unique()-> cancer_rank

ICP_exp_pattern_in_immune_tumor_cell.detailed_tissue %>%
  dplyr::rename("Tissues"="Characteristics[Tissue]") %>%
  dplyr::inner_join(TCGA_tissue,by="Tissues") %>%
  ggplot(aes(x=TCGA_Cancer),y=0) + 
  scale_x_discrete(limit = cancer_rank$TCGA_Cancer) +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "white"),
    axis.ticks = element_blank(),
    axis.text.x = element_text(vjust = 0.5, hjust = 1, size = 10, angle = 90),
    axis.title = element_blank(),
    axis.text = element_text(colour = "black"),
    plot.margin = unit(c(0,1,1,1),"cm")
  ) -> p2
library(ggpubr)
ggarrange(p1,
          p2, 
          ncol = 1, nrow = 2,  align = "hv", 
          heights = c(10, 1))

ggsave(file.path(result_path,"pattern_info","FANTOM5.ICP_exp_in_tumor_immune_stroma.detailed_tissue.png"),device = "png",height = 12,width = 6)
ggsave(file.path(result_path,"pattern_info","FANTOM5.ICP_exp_in_tumor_immune_stroma.detailed_tissue.pdf"),device = "pdf",height = 12,width = 6)

save.image(file.path(result_path,"FANTOM5.ICP_exp.pattern-byFC.rda"))
load(file.path(result_path,"FANTOM5.ICP_exp.pattern-byFC.rda"))