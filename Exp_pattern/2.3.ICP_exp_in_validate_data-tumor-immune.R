######################## verify the ICP expression site in melenoma single cell data set ############
######################## GSE72056
###### malignant tumor cell and non-malignant cell including: T cells, B cells, macrophages, CAFs(cancer-associated fibroblasts), and endothelial (Endo.) cells from preferentially expressed genes. NK, natural killer cells.

library(magrittr)

basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio")

#### gene list ------------------------------------------------------------------------
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list_exp_site <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(entrez_ID,symbol,Exp_site,`log2FC(I/T)`) %>%
  dplyr::inner_join(gene_list,by="symbol") 


# load expresion data -----------------------------------------------------

# sample info
sample_info <- readr::read_tsv(file.path(basic_path,"data/melanoma_single_cell_data","sample_anno.txt"),col_names = F) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::select(V1,V2) %>%
  dplyr::rename("sample"="V1","GSM"="V2") %>%
  .[-1,]

sample_info.class <- readr::read_tsv(file.path(basic_path,"data/melanoma_single_cell_data","GSE72056_melanoma_single_cell_revised_v2.txt")) %>%
  .[1:3,]  %>%
  tidyr::gather(-Cell,key="sample",value="Class") %>%
  tidyr::spread(key="Cell",value="Class") %>%
  # dplyr::rename("sample_type" = "malignant(1=no,2=yes,0=unresolved)","cell_type"="non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)","tumor_n" = "tumor") %>%
  dplyr::inner_join(tibble::tibble(sample_type=c("Malignant","Non-malignant","Unresolved"),
                                   `malignant(1=no,2=yes,0=unresolved)`=c(2,1,0)),by="malignant(1=no,2=yes,0=unresolved)") %>%
  dplyr::inner_join(tibble::tibble(cell_type = c("Tumor cell","T cell","B cell","Macrophage","Endothelial","CAF","NK"),
                                   `non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)`= c(0,1,2,3,4,5,6),
                                   cell_source = c("Tumor","Immune","Immune","Immune","Stromal","Stromal","Immune")),by="non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)") %>%
  dplyr::select(sample,sample_type,cell_type,cell_source,tumor) 

# exp data
ICP_exp_in_GSE72056 <- readr::read_tsv(file.path(basic_path,"data/melanoma_single_cell_data","GSE72056_melanoma_single_cell_revised_v2.txt"))
ICP_exp_in_GSE72056 %>%
  .[-c(1:3),] %>%
  dplyr::rename("symbol" ="Cell") %>%
  dplyr::filter(symbol %in% gene_list_exp_site$symbol) -> ICP_exp_in_GSE72056


# compare ICP between tumor and immune cells ------------------------------
### function to compare ICP exp between tumor and immune cells, FC and plot
fn_compare_TI_FC <- function(.data){
  # data filter
  .data %>%
    dplyr::filter(cell_source != "Stromal") -> .data
  # mean exp
  .data %>%
    dplyr::filter(cell_source == "Immune") %>%
    .$Exp %>%
    mean() -> mean_immune_exp
  .data %>%
    dplyr::filter(cell_source == "Tumor") %>%
    .$Exp %>%
    mean() -> mean_tumor_exp
  # test
  broom::tidy(
    wilcox.test(Exp ~ cell_source, data = .data, alternative = "two.sided") #Comparing the means of two independent groups:Unpaired Two-Samples Wilcoxon Test (non-parametric) 
  ) %>%
    dplyr::mutate(mean_immune_exp=mean_immune_exp,mean_tumor_exp=mean_tumor_exp) %>%
    dplyr::mutate(`log2FC(I/T)` = log2((mean_immune_exp+1)/(mean_tumor_exp+1)))
}

ICP_exp_in_GSE72056 %>%
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  tidyr::nest(-symbol) %>%
  dplyr::mutate(test = purrr::map(data,fn_compare_TI_FC)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ICP_exp_in_GSE72056.wilcox.test.FC.TI

# plot
library(ggbeeswarm)
ICP_exp_in_GSE72056 %>%
  dplyr::filter(symbol != "BTNL3") %>%
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  dplyr::filter(cell_source != "Stromal") %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") -> ready_for_draw

ready_for_draw %>%
  dplyr::select(symbol,Exp_site,`log2FC(I/T)`) %>%
  dplyr::inner_join(strip_color,by="Exp_site") %>%
  dplyr::arrange(rank,`log2FC(I/T)`) %>%
  .$symbol -> symbol_rank

ready_for_draw <- within(ready_for_draw,symbol <- factor(symbol,levels = unique(symbol_rank)))  
with(ready_for_draw, levels(symbol))

ready_for_draw %>%
  dplyr::select(symbol,Exp_site) %>%
  unique() -> color_bac
color_bac$cell_source <- color_bac$Exp <- 1


ggplot(ready_for_draw,aes(x=cell_source, y=Exp)) +
  geom_quasirandom(size=0.1) +
  geom_rect(data=color_bac,aes(fill = Exp_site),xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.1) +
  # geom_violin() +
  facet_wrap(~symbol,scales = "free_y") +
  scale_fill_manual(
    # values = site_cplor,
    values = c("yellow",  "green","pink","blue", "red"),
    # values = c("#008B00", "#00EE00", "#CD8500", "#FF4500"),
    breaks = c("Only_exp_on_Immune", "Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor","Only_exp_on_Tumor")
  ) +
  my_theme +
  ylab("Expression") +
  theme(
    axis.title.x = element_blank()
  )
ggsave(file.path(res_path,"pattern_validation","5.2.GSE72056.ICP_exp-T-I_compare.pdf"),device = "pdf",height = 8,width = 12)
ggsave(file.path(res_path,"pattern_validation","5.2.GSE72056.ICP_exp-T-I_compare.png"),device = "png",height = 8,width = 12)

# correlation between FC got from fantom and melanoma ---------------------------------
fantom_res <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(symbol,mean_cell_line, mean_immune_exp,`log2FC(I/T)`) %>%
  tidyr::gather(-symbol,key="data_type",value="Fantom5") %>%
  dplyr::mutate(data_type = ifelse(data_type == "mean_cell_line","mean_tumor_exp",data_type))

ICP_exp_in_GSE72056.wilcox.test.FC.TI %>%
  dplyr::select(symbol,mean_tumor_exp, mean_immune_exp,`log2FC(I/T)`) %>%
  tidyr::gather(-symbol,key="data_type",value="GSE72056") %>%
  dplyr::inner_join(fantom_res,by=c("symbol","data_type")) -> correlation.ready

# spearman correlation
correlation.ready %>%
  dplyr::filter(symbol != 'BTNL3') %>%
  tidyr::nest(-data_type) %>%
  dplyr::mutate(cpm_cor = purrr::map(data,.f=function(.x){
    broom::tidy(
      cor.test(.x$GSE72056,.x$Fantom5,method = "spearman")
    )
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> cor.res

# plot
library(ggplot2)
library(ggpubr)

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
strip_color <- data.frame(Exp_site = c("Only_exp_on_Immune","Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor","Only_exp_on_Tumor" ),
                          site_cplor = c("blue", "green", "yellow", "pink","red"),
                          rank = c(5,4,3,2,1))
correlation.ready %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>% 
  dplyr::group_by(data_type) %>%
  dplyr::mutate(y=(max(Fantom5)-min(Fantom5))*0.85+min(Fantom5),x=min(GSE72056)+(max(GSE72056)-min(GSE72056))*0.4) %>%
  dplyr::select(data_type,x,y) %>%
  unique() %>%
  dplyr::inner_join(cor.res,by="data_type") %>%
  dplyr::select(data_type,x,y,estimate,p.value) %>%
  dplyr::mutate(label = paste("r = ",signif(estimate,2),", p = ",signif(p.value,2),sep="")) -> cor_text

correlation.ready %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
  dplyr::filter(Exp_site!="Not_sure") %>%
  ggplot(aes(x=GSE72056,y=Fantom5)) +
  geom_jitter(aes(color=Exp_site)) +
  geom_smooth(se = F, method = "lm") +
  geom_text(aes(x=x,y=y,label = label),data=cor_text) +
  facet_wrap(~data_type,scales = "free") +
  scale_color_manual(values=c("#CD661D",  "#008B00", "#FF69B4", "#1874CD","#CD3333")) +
  my_theme 
ggsave(file.path(res_path,"pattern_validation","5.1.GSE72056-Fantom5.correlation.pdf"),device = "pdf",height = 4,width = 10)
ggsave(file.path(res_path,"pattern_validation","5.1.GSE72056-Fantom5.correlation.png"),device = "png",height = 4,width = 10)

save.image(file.path(
  res_path,"pattern_validation","GSE72056.melenoma.TI.compare.Rdata")
)
