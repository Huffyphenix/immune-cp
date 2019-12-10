######################## verify the ICP expression site in melenoma single cell data set ############
######################## GSE115978_melanoma/
###### Paper:Single-cell RNA-seq of melanoma ecosystems reveals sources of T cells exclusion linked to immunotherapy clinical outcomes
###### Description:Single-cell RNA-seq of melanoma ecosystems reveals sources of T cells exclusion linked to immunotherapy clinical outcomes. We applied our approach to 7,186 high-quality scRNA-seqprofiles from 33 human melanoma tumors (from 31 patients)comprised of 2,987 cells from 17 newly collected patient tumorsand 4,199 cells from 16 patient tumors that we previouslyreported (Tirosh et al., 2016)

library(magrittr)

basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio")
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byMeanUQ")
data_path <- file.path(basic_path,"data/single_cell_RNAseq/GSE115978_melanoma")

# load image --------------------------------------------------------------
load(file.path(
  res_path,"pattern_validation","GSE115978_melanoma.TI.compare.Rdata")
)

#### gene list ------------------------------------------------------------------------
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list_exp_site <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(entrez_ID,symbol,Exp_site,`log2FC(I/T).mean`,`log2FC(I/T).mid`) %>%
  dplyr::inner_join(gene_list,by="symbol") 


# load expresion data -----------------------------------------------------

# sample info
sample_info.class <- readr::read_csv(file.path(data_path,"GSE115978_cell.annotations.csv")) %>%
  dplyr::rename("sample"="cells","patient"="samples") %>%
  dplyr::filter(treatment.group %in% "treatment.naive") %>%
  dplyr::mutate(cell_source = ifelse(cell.types %in% c("T.cel","T.CD8","T.CD4","NK","B.cell","Macrophage"),"Immune","Unknown")) %>%
  dplyr::mutate(cell_source = ifelse(cell.types %in% c("CAF","Endo."),"Stromal",cell_source)) %>%
  dplyr::mutate(cell_source = ifelse(cell.types %in% c("Mal"),"Tumor",cell_source)) %>%
  dplyr::select(sample,cell_source)

# exp data
ICP_exp_in_GSE115978 <- readr::read_csv(file.path(data_path,"GSE115978_tpm.csv"))
ICP_exp_in_GSE115978 %>%
  dplyr::rename("symbol" ="X1") %>%
  dplyr::filter(symbol %in% gene_list_exp_site$symbol) -> ICP_exp_in_GSE115978


# compare ICP between tumor and immune cells ------------------------------
### function to compare ICP exp between tumor and immune cells, FC and plot
fn_compare_TI_FC <- function(.data){
  # data filter
  .data %>%
    dplyr::filter(cell_source  %in% c("Tumor","Immune")) -> .data
  # mean exp
  .data %>%
    dplyr::filter(cell_source == "Immune") %>%
    .$Exp %>%
    mean() -> mean_immune_exp
  .data %>%
    dplyr::filter(cell_source == "Tumor") %>%
    .$Exp %>%
    mean() -> mean_tumor_exp
  # UQ exp
  .data %>%
    dplyr::filter(cell_source == "Immune") %>%
    .$Exp %>%
    quantile(0.75) -> UQ_immune_exp
  .data %>%
    dplyr::filter(cell_source == "Tumor") %>%
    .$Exp %>%
    quantile(0.75) -> UQ_tumor_exp
  # test
  broom::tidy(
    wilcox.test(Exp ~ cell_source, data = .data, alternative = "two.sided") #Comparing the means of two independent groups:Unpaired Two-Samples Wilcoxon Test (non-parametric) 
  ) %>%
    dplyr::mutate(mean_immune_exp=mean_immune_exp,mean_tumor_exp=mean_tumor_exp,
                  UQ_immune_exp=UQ_immune_exp, UQ_tumor_exp=UQ_tumor_exp) %>%
    dplyr::mutate(`log2FC(I/T).mean` = log2((mean_immune_exp+0.01)/(mean_tumor_exp+0.01)),
                  `log2FC(I/T).UQ` = log2((UQ_immune_exp+0.01)/(UQ_tumor_exp+0.01)))
}
ICP_exp_in_GSE115978 %>%
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::mutate(Exp = as.numeric(Exp)) %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  tidyr::nest(-symbol) %>%
  dplyr::mutate(test = purrr::map(data,fn_compare_TI_FC)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ICP_exp_in_GSE115978.wilcox.test.FC.TI

ICP_exp_in_GSE115978.wilcox.test.FC.TI %>%
  dplyr::mutate(Exp_site = ifelse(`log2FC(I/T).UQ` >=0.585, "Immune cell dominate","Immune and tumor cell almost")) %>%
  dplyr::mutate(Exp_site = ifelse(`log2FC(I/T).UQ` <=(-0.585), "Tumor cell dominate",Exp_site)) %>%
  dplyr::inner_join(gene_list_exp_site, by="symbol") %>%
  dplyr::select(symbol, Exp_site.x, Exp_site.y,`log2FC(I/T).UQ`) %>% 
  dplyr::mutate(fit_fantom=ifelse(Exp_site.x==Exp_site.y,"yes","no")) %>%
  readr::write_tsv(file.path(res_path,"predict_res_validate_by_GSE115978.tsv"))

## define genes exp site by fold change and pvalue ----
# fn_define_exp_site <- function(symbol,fc,pvalue,tumor_ratio,immune_ratio,mean_cell_line, mean_immune_exp){
#   print(symbol)
#   if(is.na(pvalue)){
#     tmp <- "Not_sure"
#   } else {
#     if(fc>=1 && pvalue<=0.05){
#       if(tumor_ratio<0.25){
#         if(immune_ratio>=0.5){
#           tmp <- "Mainly_exp_on_Immune"
#         } else{
#           tmp <- "Both_exp_on_Tumor_Immune"
#         }
#       } else if(tumor_ratio>=0.25){
#         tmp <- "Both_exp_on_Tumor_Immune"
#       }
#     }else if(fc<=(-1) && pvalue<=0.05){
#       if(immune_ratio<0.25){
#         if(tumor_ratio>=0.5){
#           tmp <- "Mainly_exp_on_Tumor"
#         } else{
#           tmp <- "Both_exp_on_Tumor_Immune"
#         }
#       } else{
#         tmp <- "Both_exp_on_Tumor_Immune"
#       }
#     }else if(fc>(-1) && fc<1){
#       tmp <- "Both_exp_on_Tumor_Immune"
#     } else {
#       tmp <- "Both_exp_on_Tumor_Immune"
#     }
#   }
#   tmp
# }
# ICP_exp_in_GSE115978.wilcox.test.FC.TI %>%
#   dplyr::mutate(Exp_site = purrr::pmap(list(symbol,`log2FC(I/T)`,p.value,`tumor_ratio_diff(U-D)`,`immune_ratio_diff(U-D)`,mean_tumor_exp,mean_immune_exp),fn_define_exp_site)) %>%
#   tidyr::unnest() -> ICP_Exp_site_by_DE_Fc_and_ratio_in_GSE115978

# get p value of ICP pattern in validation data -----

fantom_res.expsite <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(symbol,Exp_site) %>%
  dplyr::rename("FANTOM_res"="Exp_site")

# ICP_Exp_site_by_DE_Fc_and_ratio_in_GSE115978 %>%
#   dplyr::select(symbol,Exp_site) %>%
#   dplyr::rename("validation_res"="Exp_site") %>%
#   dplyr::inner_join(fantom_res.expsite, by ="symbol") %>%
#   dplyr::mutate(true_pos = ifelse(validation_res==FANTOM_res,"Ture","False")) %>%
#   .$true_pos %>%
#   table() %>%
#   as.data.frame() %>%
#   readr::write_tsv(file.path(res_path,"pattern_validation","7.3.validation_accuracy.tsv"))

# plot
strip_color <- data.frame(Exp_site = unique(gene_list_exp_site$Exp_site),
                          site_cplor = c("green", "orange", "pink"),
                          rank = c(3,2,1))
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
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  dplyr::filter(cell_source %in% c("Tumor","Immune")) %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
  dplyr::inner_join(ICP_exp_in_GSE115978.wilcox.test.FC.TI,by="symbol")  -> ready_for_draw


ready_for_draw %>%
  dplyr::select(symbol,Exp_site,`log2FC(I/T).mean.x`) %>%
  dplyr::inner_join(strip_color,by="Exp_site") %>%
  dplyr::arrange(rank,`log2FC(I/T).mean.x`)-> symbol_rank

ready_for_draw <- within(ready_for_draw,symbol <- factor(symbol,levels = unique(symbol_rank$symbol)))  
with(ready_for_draw, levels(symbol))

ready_for_draw %>%
  dplyr::select(symbol,Exp_site) %>%
  unique() -> color_bac
color_bac$cell_source <- color_bac$Exp <- 1

ready_for_draw %>%
  # dplyr::mutate(cell_source = ifelse(cell_source=="Cancer cells", "Tumor","Immune")) %>%
  ggplot(aes(x=cell_source, y=Exp)) +
  geom_violin(size = 0.25) +
  geom_rect(data=color_bac,aes(fill = Exp_site),xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.1) +
  # geom_violin() +
  facet_wrap(~symbol,scales = "free_y", ncol = 7) +
  scale_fill_manual(
    # values = site_cplor,
    values = c("yellow",  "green","pink","blue", "red"),
    # values = c("#008B00", "#00EE00", "#CD8500", "#FF4500"),
    breaks = c("Only_exp_on_Immune", "Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor","Only_exp_on_Tumor")
  ) +
  my_theme +
  labs(y="Expression",title="GSE115978, melanoma") +
  theme(
    axis.title.x = element_blank(),
    legend.position = "bottom"
  )
ggsave(file.path(res_path,"pattern_validation","7.2.GSE115978.ICP_exp-T-I_compare.pdf"),device = "pdf",height =20,width = 16, units = c("cm"))
ggsave(file.path(res_path,"pattern_validation","7.2.GSE115978.ICP_exp-T-I_compare.png"),device = "png",height = 20,width = 16, units = c("cm"))
  
# correlation between FC got from fantom and melanoma ---------------------------------
fantom_res <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(symbol,mean_cell_line, mean_immune_exp,mid_cell_line,mid_immune_exp,`log2FC(I/T).mean`,`log2FC(I/T).mid`) %>%
  dplyr::rename("log2FC(I/T).UQ"="log2FC(I/T).mid","UQ_tumor_exp"="mid_cell_line","UQ_immune_exp"="mid_immune_exp","mean_tumor_exp"="mean_cell_line") %>%
  tidyr::gather(-symbol,key="data_type",value="Fantom5") 

ICP_exp_in_GSE115978.wilcox.test.FC.TI %>%
  dplyr::select(symbol,mean_tumor_exp, mean_immune_exp,UQ_immune_exp, UQ_tumor_exp,`log2FC(I/T).mean`,`log2FC(I/T).UQ`) %>%
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
  dplyr::filter(data_type %in% c("log2FC(I/T).UQ","log2FC(I/T).mean")) %>%
  ggplot(aes(x=GSE115978,y=Fantom5)) +
  geom_jitter(aes(color=Exp_site)) +
  geom_smooth(se = F, method = "lm") +
  geom_text(aes(x=x,y=y,label = label),
            data=cor_text %>% dplyr::filter(data_type %in% c("log2FC(I/T).UQ","log2FC(I/T).mean"))) +
  facet_wrap(~data_type,scales = "free") +
  scale_color_manual(values=c("#CD661D",  "#008B00", "#FF69B4", "#1874CD","#CD3333")) +
  my_theme +
  labs(x="GSE72056",
       y="FANTOM5",
       title = "FANTOM5 vs. GSE115978") +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.2,"inches"),
    legend.key.height=unit(0.2,"inches"),
    legend.text = element_text(size=8),
    legend.title = element_blank()
  )
ggsave(file.path(res_path,"pattern_validation","7.1.GSE115978-Fantom5.correlation.pdf"),device = "pdf",height = 4,width = 6)
ggsave(file.path(res_path,"pattern_validation","7.1.GSE115978-Fantom5.correlation.png"),device = "png",height = 4,width = 6)


ICP_exp_in_GSE115978.wilcox.test.FC.TI %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
  dplyr::mutate(log2Immune.UQ=log2(UQ_immune_exp+0.01),log2Tumor.UQ=log2(UQ_tumor_exp+0.01)) %>%
  ggplot(aes(x=`log2Immune.UQ`,y=`log2Tumor.UQ`)) +
  geom_jitter(aes(color = Exp_site),width = 0.1,height = 0.1) +
  geom_abline(intercept = 0.585, slope = 1,linetype=2) +
  geom_abline(intercept = -0.585, slope = 1,linetype=2) +
  geom_text(aes(x=x,y=y,label=label),
            data=tibble::tibble(x=c(2,2),
                                y=c(7,-4),
                                label=c("log2(I/T)<-0.585","log2(I/T)>0.585"))) +
  # geom_smooth(method = "lm") +
  # geom_text_repel(aes(x=`log2FC(I/T).mean`,y=`log2FC(I/T).mid`,label=symbol)) +
  # geom_label(x=4,y=10,aes(label=label),data = cor_label) +
  # geom_hline(yintercept = c(-2,2),linetype = 2) +
  # geom_vline(xintercept = c(-1,1),linetype = 2) +
  labs(x=TeX("log_2 (UQ(Immune)+0.01)"),
       y=TeX("log_2 (UQ(Tumor)+0.01)"),
       title = "Classification of ICPs' expression pattern, GSE115978") +
  scale_color_manual(values = c("#CD950C", "#66CD00", "#EE2C2C"),
                     name = "ICPs expression pattern") +
  my_theme +
  theme(
    plot.title = element_text(size=15)
  )
ggsave(file.path(res_path,"classify_ICP_exp_pattern_onlybyUQ.GSE115978.pdf"),device = "pdf",height = 4, width = 8)
ggsave(file.path(res_path,"classify_ICP_exp_pattern_onlybyUQ.GSE115978.png"),device = "png",height = 4, width = 8)

save.image(file.path(
  res_path,"pattern_validation","GSE115978_melanoma.TI.compare.Rdata")
)
