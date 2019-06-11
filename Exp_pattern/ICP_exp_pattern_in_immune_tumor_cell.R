###################################################################
# immune gene expresssion pattern in cancer and immune cell
###################################################################
library(magrittr)
library(ggplot2)

# data path ---------------------------------------------------------------
basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
result_path <- file.path(immune_path,"result_20171025","ICP_exp_patthern")

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
fn_DE_TI <- function(cell_line_exp,.symbol){
  ICP_fantom.gene_exp.cell_line.Immune_cell.combine %>%
    dplyr::filter(Group == "Immune Cell") %>%
    dplyr::filter(symbol == .symbol) -> immune_exp
  
  broom::tidy(
    wilcox.test(immune_exp$gene_tpm, cell_line_exp$gene_tpm, alternative = "two.sided") #Comparing the means of two independent groups:Unpaired Two-Samples Wilcoxon Test (non-parametric) 
  ) %>%
    dplyr::mutate(mean_cell_line = mean(cell_line_exp$gene_tpm), mean_immune_exp = mean(immune_exp$gene_tpm)) %>%
    dplyr::mutate(`log2FC(I/T)` = log2((mean_immune_exp+1)/(mean_cell_line+1)))
}

ICP_fantom.gene_exp.cell_line.Immune_cell.combine %>%
  dplyr::filter(! `Characteristics[Tissue]` %in% c("PrimaryCell","blood")) %>%
  tidyr::nest(-entrez_ID,-symbol,.key="cell_line_exp") %>%
  dplyr::mutate(DE_I_T = purrr::map2(cell_line_exp,symbol,fn_DE_TI)) %>%
  dplyr::select(-cell_line_exp) %>%
  tidyr::unnest() -> ICP_DE_FC_between_cellline_immune

## define genes exp site by fold change and pvalue ----
fn_define_exp_site <- function(symbol,fc,pvalue){
  print(symbol)
  if(is.na(pvalue)){
    tmp <- "Not_sure"
  }
  if(fc>=2){
    if(fc>=5){
      if(pvalue<=0.05){
        tmp <- "Only_exp_on_Immune"
      }
    }else if(fc<5){
      if(pvalue<=0.05){
        tmp <- "Mainly_exp_on_Immune"
      }
    }
  }else if(fc<=(-2)){
    if(fc<=(-5)){
      if(pvalue<=0.05){
        tmp <- "Only_exp_on_Tumor"
      }
    }else if(fc>(-5)){
      if(pvalue<=0.05){
        tmp <- "Mainly_exp_on_Tumor"
      }
    }
  } else {
    tmp <- "Both_exp_on_Tumor_Immune"
  }
  tmp
}

ICP_DE_FC_between_cellline_immune %>%
  dplyr::mutate(Exp_site = purrr::pmap(list(symbol,`log2FC(I/T)`,p.value),fn_define_exp_site)) %>%
  tidyr::unnest() -> ICP_Exp_site_by_DE_FC_between_cellline_immune

ICP_Exp_site_by_DE_FC_between_cellline_immune %>%
  readr::write_tsv(file.path(result_path,"ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv"))

# draw picture ------------------------------------------------------------
strip_color <- data.frame(Exp_site = c("Only_exp_on_Immune","Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor" ),
                          site_cplor = c("#458B74", "#76EEC6", "#EEC591", "#FF6A6A"))
ICP_fantom.gene_exp.cell_line.Immune_cell.combine  %>%
  dplyr::filter(! `Characteristics[Tissue]` %in% c("blood")) %>%
  # dplyr::mutate(Group = ifelse(Group=="Stromal Cell","Stromal",Group)) %>%
  dplyr::mutate(Group = ifelse(Group=="Immune Cell","Immune",Group)) %>%
  dplyr::mutate(Group = ifelse(Group=="Tumor Cell","Tumor",Group)) %>%
  dplyr::mutate(gene_tpm = log2(gene_tpm+1)) %>%
  dplyr::inner_join(ICP_Exp_site_by_DE_FC_between_cellline_immune,by=c("symbol","entrez_ID")) %>%
  dplyr::inner_join(strip_color,by="Exp_site") %>%
  dplyr::filter(!is.na(`Characteristics[Tissue]`)) %>%
  dplyr::mutate(symbol = paste(symbol, ", FC =",signif(`log2FC(I/T)`,2),sep = ""))-> ready_for_draw

ready_for_draw %>%
  dplyr::arrange(`log2FC(I/T)`) %>%
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
  ggpubr::stat_compare_means(method = "wilcox.test",label = "p.format") +
  # geom_text(aes(label = outlier), na.rm = TRUE, nudge_x = -0.25, nudge_y = 0.25,size=3) +
  scale_fill_manual(
    values = c("yellow",  "green","red", "blue"),
    # values = c("#008B00", "#00EE00", "#CD8500", "#FF4500"),
    breaks = c("Only_exp_on_Immune", "Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor")
  ) +
  # theme_bw() +
  ylab("log2(TPM)") +
  theme(
    panel.background = element_rect(fill = "white",colour = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(vjust = 1, hjust = 1, size = 10),
    axis.title.x = element_blank(),
    # legend.position = "none",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(size = 20),
    axis.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "white",colour = "black"),
    strip.text = element_text(size = 12),
    text = element_text(color = "black")
) 
ggsave(file.path(result_path,"FANTOM5.ICP_exp_in_tumor(no-blood-tissue)_immune-TCGA.tissue.png"),device = "png",height = 10,width = 15)  
ggsave(file.path(result_path,"FANTOM5.ICP_exp_in_tumor(no-blood-tissue)_immune-TCGA.tissue.pdf"),device = "pdf",height = 10,width = 15)  

ggsave(file.path(result_path,"FANTOM5.ICP_exp_in_tumor_immune_stroma-all.tissue.png"),device = "png",height = 10,width = 15)  
ggsave(file.path(result_path,"FANTOM5.ICP_exp_in_tumor_immune_stroma-all.tissue.pdf"),device = "pdf",height = 10,width = 15)

# global exp pattern ------------------------------------------------------

fn_exp_pantern_classify <- function(.n,.x){
  print(.n)
  .x %>% 
    dplyr::filter(Group == "Stromal Cell") %>%
    .$gene_mean_exp -> Stromal_exp
  .x %>% 
    dplyr::filter(Group == "Immune Cell") %>%
    .$gene_mean_exp -> Immune_exp
  .x %>% 
    dplyr::filter(Group == "Tumor Cell") %>%
    .$gene_mean_exp  %>%
    quantile(0.5) %>%
    as.numeric()-> Tumor_exp
  if(Tumor_exp==0){Tumor_exp <- 0.01}
  if(Immune_exp==0){Immune_exp <- 0.01}
  log2fc <- log2(Immune_exp/Tumor_exp)
  
  if(log2fc >= 1){
    patt <- "Immune Higher"
  } else if(log2fc <= -1){
    patt <- "Tumor Higher"
  } else if(log2fc==0){
    patt <- "Little diff"
  } else {
    patt <- "Little diff"
  }
  return(patt)
}

ICP_fantom.gene_exp.cell_line.Immune_cell.combine %>%
  dplyr::group_by(symbol) %>%
  dplyr::select(-hgnc_id) %>%
  tidyr::nest() %>%
  dplyr::mutate(Exp_pattern = purrr::map2(symbol,data,fn_exp_pantern_classify)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ICP_exp_pattern_in_immune_tumor_cell

ICP_exp_pattern_in_immune_tumor_cell %>%
  dplyr::arrange(symbol) %>%
  readr::write_tsv(file.path(result_path,"ICP_exp_pattern_in_immune_tumor_cell.tsv"))

# detailed tissue exp pattern ------------------------------------------------------

fn_exp_detailed_pantern_classify <- function(.n,.x){
  print(.n)
  .x %>% 
    dplyr::filter(Group == "Stromal Cell") %>%
    .$gene_mean_exp -> Stromal_exp
  .x %>% 
    dplyr::filter(Group == "Immune Cell") %>%
    .$gene_mean_exp -> Immune_exp
  if(Immune_exp==0){Immune_exp <- 0.01}
  
  .x %>% 
    dplyr::filter(Group == "Tumor Cell") %>%
    dplyr::group_by(`Characteristics[Tissue]`) %>%
    dplyr::mutate(gene_mean_exp = ifelse(gene_mean_exp==0,0.01,gene_mean_exp)) %>%
    dplyr::mutate(Immune_exp = Immune_exp) %>%
    dplyr::mutate(log2FC = log2(Immune_exp/gene_mean_exp)) %>%
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
  dplyr::group_by(symbol) %>%
  dplyr::select(-hgnc_id) %>%
  tidyr::nest() %>%
  dplyr::mutate(Exp_pattern = purrr::map2(symbol,data,fn_exp_detailed_pantern_classify)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ICP_exp_pattern_in_immune_tumor_cell.detailed_tissue

ICP_exp_pattern_in_immune_tumor_cell.detailed_tissue %>%
  readr::write_tsv(file.path(result_path,"ICP_exp_pattern_in_immune_tumor_cell.detailed_tissues.tsv"))

# plot 
# classification of tcga cancers

ICP_expr_pattern <- readr::read_tsv(file.path(result_path,"manual_edit_2_ICP_exp_pattern_in_immune_tumor_cell.tsv"))
fn_site_color <- function(.x){
  # print(.n)
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
ICP_expr_pattern %>%
  dplyr::filter(!is.na(`Exp site`)) %>%
  dplyr::mutate(site_col = purrr::map(`Exp site`,fn_site_color)) %>%
  tidyr::unnest() %>%
  dplyr::arrange(`Exp site`,symbol) -> gene_rank

ICP_exp_pattern_in_immune_tumor_cell.detailed_tissue %>%
  dplyr::filter(`Characteristics[Tissue]` %in% c(TCGA_tissue$Tissues %>% unique())) %>%
  ggplot(aes(x=`Characteristics[Tissue]`,y=symbol)) +
  geom_tile(aes(fill=log2FC),color = "grey",size=0.5,width = 0.9) +
  scale_fill_gradient2(
    name = "log2 (I/T)",
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

ggarrange(p1,
          p2, 
          ncol = 1, nrow = 2,  align = "hv", 
          heights = c(10, 1))

ggsave(file.path(result_path,"FANTOM5.ICP_exp_in_tumor_immune_stroma.detailed_tissue.png"),device = "png",height = 12,width = 6)
ggsave(file.path(result_path,"FANTOM5.ICP_exp_in_tumor_immune_stroma.detailed_tissue.pdf"),device = "pdf",height = 12,width = 6)


save.image(file.path(result_path,"FANTOM5.ICP_exp.pattern.rda"))
load(file.path(result_path,"FANTOM5.ICP_exp.pattern.rda"))
