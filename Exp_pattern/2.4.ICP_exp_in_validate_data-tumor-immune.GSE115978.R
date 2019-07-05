######################## verify the ICP expression site in melenoma single cell data set ############
######################## GSE115978_melanoma/
###### Paper:Single-cell RNA-seq of melanoma ecosystems reveals sources of T cells exclusion linked to immunotherapy clinical outcomes
###### Description:Single-cell RNA-seq of melanoma ecosystems reveals sources of T cells exclusion linked to immunotherapy clinical outcomes. We applied our approach to 7,186 high-quality scRNA-seqprofiles from 33 human melanoma tumors (from 31 patients)comprised of 2,987 cells from 17 newly collected patient tumorsand 4,199 cells from 16 patient tumors that we previouslyreported (Tirosh et al., 2016)

library(magrittr)

basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio")
data_path <- file.path(basic_path,"data/single_cell_RNAseq/GSE115978_melanoma")

# load image --------------------------------------------------------------
load(file.path(
  res_path,"pattern_validation","GSE115978_melanoma.TI.compare.Rdata")
)

#### gene list ------------------------------------------------------------------------
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list_exp_site <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(entrez_ID,symbol,Exp_site,`log2FC(I/T)`) %>%
  dplyr::inner_join(gene_list,by="symbol") 


# load expresion data -----------------------------------------------------

# sample info
sample_info.class <- readr::read_csv(file.path(data_path,"GSE115978_cell.annotations.csv")) %>%
  dplyr::rename("sample"="cells","patient"="samples") %>%
  dplyr::filter(treatment.group %in% "treatment.naive") %>%
  dplyr::mutate(cell_source = ifelse(cell.types %in% c("T.cel","T.CD8","T.CD4","NK","B.cell","Macrophage"),"Immune cells","Unknown")) %>%
  dplyr::mutate(cell_source = ifelse(cell.types %in% c("CAF","Endo."),"Stromal cells",cell_source)) %>%
  dplyr::mutate(cell_source = ifelse(cell.types %in% c("Mal"),"Cancer cells",cell_source)) %>%
  dplyr::select(sample,cell_source)

# exp data
ICP_exp_in_GSE115978 <- readr::read_csv(file.path(data_path,"GSE115978_tpm.csv"))
ICP_exp_in_GSE115978 %>%
  dplyr::rename("symbol" ="X1") %>%
  dplyr::filter(symbol %in% gene_list_exp_site$symbol) -> ICP_exp_in_GSE115978


# compare ICP between tumor and immune cells ------------------------------
### function to compare ICP exp between tumor and immune cells, FC and plot
fn_compare_TI_FC <- function(.data,cell_type){
  # data filter
  .data %>%
    dplyr::filter(cell_source  %in% cell_type) -> .data
  # mean exp
  .data %>%
    dplyr::filter(cell_source == cell_type[1]) %>%
    .$Exp %>%
    mean() -> mean_immune_exp
  .data %>%
    dplyr::filter(cell_source == cell_type[2]) %>%
    .$Exp %>%
    mean() -> mean_tumor_exp
  # test
  broom::tidy(
    wilcox.test(Exp ~ cell_source, data = .data, alternative = "two.sided") #Comparing the means of two independent groups:Unpaired Two-Samples Wilcoxon Test (non-parametric) 
  ) %>%
    dplyr::mutate(mean_immune_exp=mean_immune_exp,mean_tumor_exp=mean_tumor_exp) %>%
    dplyr::mutate(`log2FC(I/T)` = log2((mean_immune_exp+1)/(mean_tumor_exp+1)))
}

ICP_exp_in_GSE115978 %>%
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::mutate(Exp = as.numeric(Exp)) %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  tidyr::nest(-symbol) %>%
  dplyr::mutate(test = purrr::map(data,fn_compare_TI_FC,cell_type=c("Immune cells","Cancer cells"))) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ICP_exp_in_GSE115978.wilcox.test.FC.TI

# plot
strip_color <- data.frame(Exp_site = c("Only_exp_on_Immune","Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor","Only_exp_on_Tumor" ),
                          site_cplor = c("blue", "green", "orange", "pink","red"),
                          rank = c(5,4,3,2,1))

my_theme <-   theme(
  panel.background = element_rect(fill = "white",colour = "black"),
  panel.grid.major=element_blank(),
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
library(ggbeeswarm)
ICP_exp_in_GSE115978 %>%
  dplyr::filter(symbol != "BTNL3") %>%
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::mutate(Exp = as.numeric(Exp)) %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  dplyr::filter(cell_source %in% c("Cancer cells","Immune cells")) %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
  dplyr::inner_join(ICP_exp_in_GSE115978.wilcox.test.FC.TI,by="symbol") %>%
  dplyr::mutate(label = paste(symbol,"log2FC=",signif(`log2FC(I/T).y`,2))) -> ready_for_draw

ready_for_draw %>%
  dplyr::select(label,Exp_site,`log2FC(I/T).x`) %>%
  dplyr::inner_join(strip_color,by="Exp_site") %>%
  dplyr::arrange(rank,`log2FC(I/T).x`) %>%
  .$label -> symbol_rank

ready_for_draw <- within(ready_for_draw,label <- factor(label,levels = unique(symbol_rank)))  
with(ready_for_draw, levels(label))

ready_for_draw %>%
  dplyr::select(label,Exp_site) %>%
  unique() -> color_bac
color_bac$cell_source <- color_bac$Exp <- 1


ggplot(ready_for_draw,aes(x=cell_source, y=Exp)) +
  geom_quasirandom(size=0.1) +
  geom_rect(data=color_bac,aes(fill = Exp_site),xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.1) +
  # geom_violin() +
  facet_wrap(~label,scales = "free_y") +
  scale_fill_manual(
    # values = site_cplor,
    values = c("yellow",  "green","pink","blue", "red"),
    # values = c("#008B00", "#00EE00", "#CD8500", "#FF4500"),
    breaks = c("Only_exp_on_Immune", "Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor","Only_exp_on_Tumor")
  ) +
  my_theme +
  ylab("Expression") +
  theme(
    axis.title.x = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5)
  )
ggsave(file.path(res_path,"pattern_validation","7.2.GSE115978.ICP_exp-T-I_compare.pdf"),device = "pdf",height = 8,width = 14)
ggsave(file.path(res_path,"pattern_validation","7.2.GSE115978.ICP_exp-T-I_compare.png"),device = "png",height = 8,width = 14)
  
# correlation between FC got from fantom and melanoma ---------------------------------
fantom_res <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(symbol,mean_cell_line, mean_immune_exp,`log2FC(I/T)`) %>%
  tidyr::gather(-symbol,key="data_type",value="Fantom5") %>% 
  dplyr::mutate(data_type = ifelse(data_type == "mean_cell_line","mean_tumor_exp",data_type))

ICP_exp_in_GSE115978.wilcox.test.FC.TI %>%
  dplyr::select(symbol,mean_tumor_exp, mean_immune_exp,`log2FC(I/T)`) %>%
  tidyr::gather(-symbol,key="data_type",value="GSE115978") %>%
  dplyr::inner_join(fantom_res,by=c("symbol","data_type")) -> correlation.ready

# spearman correlation
correlation.ready %>%
  dplyr::filter(symbol != 'BTNL3') %>%
  tidyr::nest(-data_type) %>%
  dplyr::mutate(cpm_cor = purrr::map(data,.f=function(.x){
    broom::tidy(
      cor.test(.x$GSE115978,.x$Fantom5,method = "spearman")
    )
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> cor.res

# plot
library(ggplot2)
library(ggpubr)

correlation.ready %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>% 
  dplyr::group_by(data_type) %>%
  dplyr::mutate(y=(max(Fantom5)-min(Fantom5))*0.85+min(Fantom5),x=min(GSE115978)+(max(GSE115978)-min(GSE115978))*0.4) %>%
  dplyr::select(data_type,x,y) %>%
  unique() %>%
  dplyr::inner_join(cor.res,by="data_type") %>%
  dplyr::select(data_type,x,y,estimate,p.value) %>%
  dplyr::mutate(label = paste("r = ",signif(estimate,2),", p = ",signif(p.value,2),sep="")) %>%
  dplyr::ungroup() -> cor_text

correlation.ready %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
  dplyr::filter(Exp_site!="Not_sure") %>%
  dplyr::filter(data_type=="log2FC(I/T)") %>%
  ggplot(aes(x=GSE115978,y=Fantom5)) +
  geom_jitter(aes(color=Exp_site)) +
  geom_smooth(se = F, method = "lm") +
  geom_text(aes(x=x,y=y,label = label),
            data=cor_text %>%
              dplyr::filter(data_type=="log2FC(I/T)")) +
  # facet_wrap(~data_type,scales = "free") +
  scale_color_manual(values=c("#CD661D",  "#008B00", "#FF69B4", "#1874CD","#CD3333")) +
  my_theme +
  labs(x="Log2 expression fold change of immune regulators\nbetween immune cells and tumor cells [GSE115978]",
       y="Log2 expression fold change of immune regulators\nbetween immune cells and tumor cells [FANTOM5]") +
  theme(
    legend.position = "bottom"
  )
ggsave(file.path(res_path,"pattern_validation","7.1.GSE115978-Fantom5.correlation.pdf"),device = "pdf",height = 4,width = 5)
ggsave(file.path(res_path,"pattern_validation","7.1.GSE115978-Fantom5.correlation.png"),device = "png",height = 4,width = 5)

save.image(file.path(
  res_path,"pattern_validation","GSE115978_melanoma.TI.compare.Rdata")
)
