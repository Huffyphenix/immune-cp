######################## verify the ICP expression site in head and neck single cell data set ############
######################## GSE103322_HNSCC
###### To understand the diversity of expression states within head and neck cancers, we profiled 5902 single cells from 18 patients with oral cavity tumors by single cell RNA-seq.

library(magrittr)

basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio")
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byMeanUQ")
data_path <- file.path(basic_path,"data/single_cell_RNAseq/GSE103322_HNSCC")

# load image --------------------------------------------------------------
load(file.path(
  res_path,"pattern_validation","GSE103322_HNSCC.TI.compare.Rdata")
)

#### gene list ------------------------------------------------------------------------
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list_exp_site <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(entrez_ID,symbol,Exp_site,`log2FC(I/T).mean`,`log2FC(I/T).mid`) %>%
  dplyr::inner_join(gene_list,by="symbol") 


# load expresion data -----------------------------------------------------

# sample info
sample_info <- readr::read_tsv(file.path(data_path,"sample_anno.txt"),col_names = F) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::select(V1,V2) %>%
  dplyr::rename("sample"="V1","GSM"="V2") %>%
  .[-1,]

sample_info.class <- readr::read_tsv(file.path(data_path,"cell_types.class.txt")) %>%
  tidyr::gather(-X1,key="sample",value="Class") %>%
  tidyr::spread(key="X1",value="Class") %>%
  dplyr::mutate(cell_source = ifelse(`classified  as cancer cell`==1,"Tumor","Non cancer")) %>%
  dplyr::mutate(cell_source = ifelse(`non-cancer cell type` %in% c("T cell","B cell","Macrophage","Mast", "Dendritic","myocyte"),"Immune",cell_source)) %>%
  dplyr::mutate(cell_source = ifelse(`non-cancer cell type` %in% c( "-Fibroblast","Endothelial"),"Stromal",cell_source)) %>%
  dplyr::select(sample,cell_source)

# exp data
ICP_exp_in_GSE103322 <- readr::read_tsv(file.path(data_path,"GSE103322_HNSCC_all_data.txt"))
ICP_exp_in_GSE103322 %>%
  .[-c(1:5),] %>%
  dplyr::rename("symbol" ="X1") %>%
  dplyr::mutate(symbol = purrr::map(symbol,.f=function(.x){gsub("\\'","",.x)})) %>%
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% gene_list_exp_site$symbol) -> ICP_exp_in_GSE103322


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

ICP_exp_in_GSE103322 %>%
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::mutate(Exp = as.numeric(Exp)) %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  tidyr::nest(-symbol) %>%
  dplyr::mutate(test = purrr::map(data,fn_compare_TI_FC)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ICP_exp_in_GSE103322.wilcox.test.FC.TI

ICP_exp_in_GSE103322.wilcox.test.FC.TI  %>%
  dplyr::mutate(Exp_site = ifelse(`log2FC(I/T).UQ` >=1, "Immune cell dominate","Immune and tumor cell almost")) %>%
  dplyr::mutate(Exp_site = ifelse(`log2FC(I/T).UQ` <=(-1), "Tumor cell dominate",Exp_site)) %>%
  dplyr::inner_join(gene_list_exp_site, by="symbol") %>%
  dplyr::select(symbol, Exp_site.x, Exp_site.y,`log2FC(I/T).UQ`) %>% 
  dplyr::mutate(fit_fantom=ifelse(Exp_site.x==Exp_site.y,"yes","no")) %>%
  readr::write_tsv(file.path(res_path,"predict_res_validate_by_GSE103322.tsv"))

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
# ICP_exp_in_GSE103322.wilcox.test.FC.TI %>%
#   dplyr::mutate(Exp_site = purrr::pmap(list(symbol,`log2FC(I/T)`,p.value,`tumor_ratio_diff(U-D)`,`immune_ratio_diff(U-D)`,mean_tumor_exp,mean_immune_exp),fn_define_exp_site)) %>%
#   tidyr::unnest() -> ICP_Exp_site_by_DE_Fc_and_ratio_in_GSE103322

# get p value of ICP pattern in validation data -----

fantom_res.expsite <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(symbol,Exp_site) %>%
  dplyr::rename("FANTOM_res"="Exp_site")

# ICP_Exp_site_by_DE_Fc_and_ratio_in_GSE103322 %>%
#   dplyr::select(symbol,Exp_site) %>%
#   dplyr::rename("validation_res"="Exp_site") %>%
#   dplyr::inner_join(fantom_res.expsite, by ="symbol") %>%
#   dplyr::mutate(true_pos = ifelse(validation_res==FANTOM_res,"Ture","False")) %>%
#   .$true_pos %>%
#   table() %>%
#   as.data.frame() %>%
#   readr::write_tsv(file.path(res_path,"pattern_validation","6.3.validation_accuracy.tsv"))

# boxplot ------------------------
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
ICP_exp_in_GSE103322 %>%
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  dplyr::filter(cell_source %in% c("Tumor","Immune")) %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
  dplyr::inner_join(ICP_exp_in_GSE103322.wilcox.test.FC.TI,by="symbol")  -> ready_for_draw

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
  # geom_quasirandom(size=0.1) +
  geom_violin(size = 0.25) +
  geom_rect(data=color_bac,aes(fill = Exp_site),xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.1) +
  # geom_violin() +
  facet_wrap(~symbol,scales = "free_y", ncol = 7) +
  scale_fill_manual(
    # values = site_cplor,
    values = c("yellow","green","pink","blue", "red"),
    # values = c("#008B00", "#00EE00", "#CD8500", "#FF4500"),
    breaks =c("Immune and tumor cell almost", "Immune cell dominate","Tumor cell dominate")
  ) +
  my_theme +
  labs(y="Expression",title="GSE103322, HNSCC") +
  theme(
    axis.title.x = element_blank(),
    legend.position = "bottom"
  )
ggsave(file.path(res_path,"pattern_validation","6.2.GSE103322.ICP_exp-T-I_compare.pdf"),device = "pdf",height = 20,width = 16, units = c("cm"))
ggsave(file.path(res_path,"pattern_validation","6.2.GSE103322.ICP_exp-T-I_compare.png"),device = "png",height = 20,width = 16, units = c("cm"))

# correlation between FC got from fantom and melanoma ---------------------------------
fantom_res <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(symbol,mean_cell_line, mean_immune_exp,mid_cell_line,mid_immune_exp,`log2FC(I/T).mean`,`log2FC(I/T).mid`) %>%
  dplyr::rename("log2FC(I/T).UQ"="log2FC(I/T).mid","UQ_tumor_exp"="mid_cell_line","UQ_immune_exp"="mid_immune_exp","mean_tumor_exp"="mean_cell_line") %>%
  tidyr::gather(-symbol,key="data_type",value="Fantom5") 

ICP_exp_in_GSE103322.wilcox.test.FC.TI %>%
  dplyr::select(symbol,mean_tumor_exp, mean_immune_exp,UQ_immune_exp, UQ_tumor_exp,`log2FC(I/T).mean`,`log2FC(I/T).UQ`) %>%
  tidyr::gather(-symbol,key="data_type",value="GSE103322") %>%
  dplyr::inner_join(fantom_res,by=c("symbol","data_type")) -> correlation.ready

# spearman correlation
correlation.ready %>%
  # dplyr::filter(symbol != 'BTNL3') %>%
  tidyr::nest(-data_type) %>%
  dplyr::mutate(cpm_cor = purrr::map(data,.f=function(.x){
    broom::tidy(
      cor.test(.x$GSE103322,.x$Fantom5,method = "spearman")
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
  dplyr::mutate(y=(max(Fantom5)-min(Fantom5))*0.85+min(Fantom5),x=min(GSE103322)+(max(GSE103322)-min(GSE103322))*0.4) %>%
  dplyr::select(data_type,x,y) %>%
  unique() %>%
  dplyr::inner_join(cor.res,by="data_type") %>%
  dplyr::select(data_type,x,y,estimate,p.value) %>%
  dplyr::mutate(label = paste("r = ",signif(estimate,2),", p = ",signif(p.value,2),sep="")) -> cor_text

correlation.ready %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
  dplyr::filter(Exp_site!="Not_sure", data_type %in% c("log2FC(I/T).mean","log2FC(I/T).UQ")) %>%
  ggplot(aes(x=GSE103322,y=Fantom5)) +
  geom_jitter(aes(color=Exp_site)) +
  geom_smooth(se = F, method = "lm") +
  geom_text(aes(x=x,y=y,label = label),data=cor_text %>% dplyr::filter(data_type == "log2FC(I/T)")) +
  facet_wrap(~data_type,scales = "free") +
  scale_color_manual(values=c("#CD661D",  "#008B00", "#FF69B4", "#1874CD","#CD3333")) +
  my_theme +
  labs(x="GSE72056",
       y="FANTOM5",
       title = "FANTOM5 vs. GSE103322") +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.2,"inches"),
    legend.key.height=unit(0.2,"inches"),
    legend.text = element_text(size=8),
    legend.title = element_blank()
  )
ggsave(file.path(res_path,"pattern_validation","6.1.GSE103322-Fantom5.correlation.pdf"),device = "pdf",height = 4,width = 6)
ggsave(file.path(res_path,"pattern_validation","6.1.GSE103322-Fantom5.correlation.png"),device = "png",height = 4,width = 6)

ICP_exp_in_GSE103322.wilcox.test.FC.TI %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
  dplyr::mutate(log2Immune.UQ=log2(UQ_immune_exp+0.01),log2Tumor.UQ=log2(UQ_tumor_exp+0.01)) %>%
  ggplot(aes(x=`log2Immune.UQ`,y=`log2Tumor.UQ`)) +
  geom_jitter(aes(color = Exp_site),width = 0.1,height = 0.1) +
  geom_abline(intercept = 1, slope = 1,linetype = 2) +
  geom_abline(intercept = -1, slope = 1,linetype = 2) +
  geom_text(aes(x=x,y=y,label=label),
            data=tibble::tibble(x=c(2,2),
                                y=c(7,-4),
                                label=c("log2(I/T)<-1","log2(I/T)>1"))) +
  # geom_smooth(method = "lm") +
  # geom_text_repel(aes(x=`log2FC(I/T).mean`,y=`log2FC(I/T).mid`,label=symbol)) +
  # geom_label(x=4,y=10,aes(label=label),data = cor_label) +
  # geom_hline(yintercept = c(-2,2),linetype = 2) +
  # geom_vline(xintercept = c(-1,1),linetype = 2) +
  labs(x=TeX("log_2 (UQ(Immune)+0.01)"),
       y=TeX("log_2 (UQ(Tumor)+0.01)"),
       title = "Classification of ICPs' expression pattern, GSE103322") +
  scale_color_manual(values = c("#CD950C", "#66CD00", "#EE2C2C"),
                     name = "ICPs expression pattern") +
  my_theme +
  theme(
    plot.title = element_text(size=15)
  )
ggsave(file.path(res_path,"classify_ICP_exp_pattern_onlybyUQ.GSE103322.pdf"),device = "pdf",height = 4, width = 8)
ggsave(file.path(res_path,"classify_ICP_exp_pattern_onlybyUQ.GSE103322.png"),device = "png",height = 4, width = 8)

# save image --------------------------------------------------------------
save.image(file.path(
  res_path,"pattern_validation","GSE103322_HNSCC.TI.compare.Rdata")
)
#>>>>>>>>>>>>>>>>>>>>>>> HAVE NOT RUN 
# tSNE: use ICP exp to distingrush tumor and immune cells -----------------
library("Rtsne")
# filter repeat samples
ICP_exp_in_GSE103322 %>%
  tidyr::gather(-symbol,key="sample",value="exp") %>% 
  dplyr::group_by(sample) %>% 
  dplyr::mutate(mean = mean(exp)) %>% 
  dplyr::select(sample,mean) %>% 
  dplyr::filter(mean == 0) %>% 
  unique() -> duplicated_samples # with same ICP exp
# data prepare
ICP_exp_in_GSE103322 %>%
  tidyr::gather(-symbol,key="sample",value="exp") %>%
  dplyr::filter(!sample %in% duplicated_samples$sample) %>%
  tidyr::spread(key="symbol",value="exp") %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  dplyr::filter(cell_source %in% c("Cancer cells","Immune cells")) -> data_for_tSNE

# data_for_tSNE %>%
#   dplyr::select(-sample,-sample_type,-cell_type,-cell_source,-tumor) %>%
#   as.data.frame() %>%
#   as.matrix() -> data_for_tSNE.matx
# normalization
data_for_tSNE.matx.normalize <- normalize_input(data_for_tSNE[,2:66] %>% as.matrix())
# colMeans(data_for_tSNE.matx.normalize)
# range(data_for_tSNE.matx.normalize)
# do tSNE
tsne_res <- Rtsne(data_for_tSNE.matx.normalize,dims=2,pca=FALSE,theta=0.0)
# Show the objects in the 2D tsne representation
plot(tsne_res$Y,col=factor(data_for_tSNE$cell_source), asp=1)
tsne_res$Y %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  # dplyr::rename("tSNE 1"="V1","tSNE 2"="V2") %>%
  dplyr::mutate(sample = data_for_tSNE$sample,
                cell_source = data_for_tSNE$cell_source,
                cell_type = data_for_tSNE$cell_type,
                tumor = data_for_tSNE$tumor) %>%
  ggplot(aes(x=V1,y= V2)) +
  geom_jitter(aes(color=cell_type),size=0.5) +
  xlab("tSNE 1") +
  ylab("tSNE 2") +
  ggpubr::color_palette("jco") +
  my_theme
ggsave(filename = file.path(res_path,"pattern_validation","5.3.GSE72056.ICP_exp-T-I_tSNE.png"),device = "png",width = 6,height = 4)
ggsave(filename = file.path(res_path,"pattern_validation","5.3.GSE72056.ICP_exp-T-I_tSNE.pdf"),device = "pdf",width = 6,height = 4)

# why some tumor cells overlap with immune cell?
tsne_res$Y %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  # dplyr::rename("tSNE 1"="V1","tSNE 2"="V2") %>%
  dplyr::mutate(sample = data_for_tSNE$sample,
                cell_source = data_for_tSNE$cell_source,
                cell_type = data_for_tSNE$cell_type,
                tumor = paste("Mel",data_for_tSNE$tumor)) %>%
  dplyr::mutate(cell_type = ifelse(cell_type!="Tumor cell","Immune cell",cell_type)) %>%
  dplyr::mutate(tumor = ifelse(cell_type== "Immune cell","Immune cell",tumor)) %>%
  ggplot(aes(x=V1,y= V2)) +
  geom_jitter(aes(color=tumor),size=0.2) +
  xlab("tSNE 1") +
  ylab("tSNE 2") +
  ggpubr::color_palette(palette = c("grey77", "#000000", "#0000FF", "#8B2323", "#CDAA7D", "#8EE5EE", "#76EE00","#D2691E", "#8B008B", "#6E8B3D","#FAEBD7", "#006400", "#FFD700", "#EE00EE", "#FFB6C1", "#FFBBFF", "#00F5FF", "#76EEC6","#DB7093", "#FF3030"),name="Sample type") +
  my_theme
ggsave(filename = file.path(res_path,"pattern_validation","5.3.GSE72056.ICP_exp-T-I_tSNE-patients-colored.png"),device = "png",width = 10,height = 6)
ggsave(filename = file.path(res_path,"pattern_validation","5.3.GSE72056.ICP_exp-T-I_tSNE-patients-colored.pdf"),device = "pdf",width = 10,height = 6)
# PCA analysis ------------------------------------------------------------

library("FactoMineR")
library("factoextra")
ICP_exp_in_GSE72056 %>%
  tidyr::gather(-symbol,key="sample",value="exp") %>%
  dplyr::arrange(symbol) %>%
  dplyr::filter(!sample %in% duplicated_samples$sample) %>%
  tidyr::spread(key="symbol",value="exp") %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  dplyr::filter(cell_source %in% c("Tumor","Immune")) -> data_for_PCA

res.pca <- PCA(data_for_PCA[,2:66] %>% as.matrix(), scale.unit = T,graph = FALSE)
# plot -- variables
fviz_pca_var(res.pca, 
             # geom.var = "text",
             labelsize = 2, 
             col.var = ICP_exp_in_GSE72056 %>%
               dplyr::select(symbol) %>%
               dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
               dplyr::inner_join(data.frame(Exp_site = c("Only_exp_on_Immune","Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor","Only_exp_on_Tumor","Not_sure"),
                                            rank = c(5,4,3,2,1,0)), by = "Exp_site") %>%
               # dplyr::filter(Exp_site!="Not_sure") %>%
               dplyr::arrange(symbol) %>%
               .$Exp_site,
             repel = TRUE,
             legend.title = list(color = "Exp. site")
) +
  my_theme +
  ggpubr::color_palette("npg")     # Variable colors
ggsave(filename = file.path(res_path,"pattern_validation","5.4.1.GSE72056.ICP_exp-T-I_PCA-variables.png"),device = "png",width = 6,height = 4)
ggsave(filename = file.path(res_path,"pattern_validation","5.4.1.GSE72056.ICP_exp-T-I_PCA-variables.pdf"),device = "pdf",width = 6,height = 4)

# plot -- Individuas
fviz_pca_ind(res.pca, 
             geom.ind = "point",
             pointsize = 2.5,
             pointshape = 21,
             fill.ind = data_for_PCA$cell_type,
             legend.title = list(fill = "Cell type")
)+
  my_theme +
  ggpubr::fill_palette("jco")
ggsave(filename = file.path(res_path,"pattern_validation","5.4.2.GSE72056.ICP_exp-T-I_PCA-Individuals.png"),device = "png",width = 6,height = 4)
ggsave(filename = file.path(res_path,"pattern_validation","5.4.2.GSE72056.ICP_exp-T-I_PCA-Individuals.pdf"),device = "pdf",width = 6,height = 4)


# vertical comparison of ICPs in tumor and immune -------------------------
fn_plot_ICP_exp_in_dataset <- function(.data,ylab,facet,title,filename){
  .data %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(mid = quantile(Exp,0.5)) %>%
    dplyr::arrange(mid) %>%
    dplyr::select(symbol,mid) %>%
    unique() %>%
    dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
    dplyr::inner_join(strip_color,by="Exp_site") -> .symbol_rank
  .data %>%
    dplyr::mutate(Exp = log2(Exp+1)) %>%
    ggplot(aes(x=symbol,y=Exp)) +
    geom_boxplot(outlier.colour = "grey",outlier.size = 0.5) +
    facet_wrap(as.formula(facet)) +
    rotate() +
    ggtitle(title) +
    ylab(ylab) +
    xlab("Symbol") +
    scale_x_discrete(limits = .symbol_rank$symbol) +
    my_theme +
    theme(
      axis.text.y = element_text(size = 10,colour = .symbol_rank$site_cplor),
      plot.title = element_text(size=12)
    )
  ggsave(file.path(res_path,"pattern_validation",paste(filename,"pdf",sep=".")),device = "pdf",width = 4, height = 10)
  ggsave(file.path(res_path,"pattern_validation",paste(filename,"png",sep=".")),device = "png", width = 4, height = 10)
}

ICP_exp_in_GSE72056 %>%
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::arrange(symbol) %>%
  dplyr::filter(!sample %in% duplicated_samples$sample) %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  dplyr::filter(cell_source %in% c("Tumor","Immune")) %>%
  fn_plot_ICP_exp_in_dataset(ylab="Expression",facet="~cell_source",title="GSE72056, ICP expression",filename="5.5.GSE72056.ICP_exp.T-I")


# correlation between ICP exp in tumor and immune -------------------------
ICP_exp_in_GSE72056 %>%
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::arrange(symbol) %>%
  dplyr::filter(!sample %in% duplicated_samples$sample) %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  dplyr::filter(cell_source %in% c("Tumor","Immune")) %>%
  dplyr::group_by(symbol,cell_source) %>%
  dplyr::mutate(Mean_exp = mean(Exp)) %>%
  dplyr::select(symbol,cell_source,Mean_exp) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="cell_source",value="Mean_exp") -> ready_for_cor

ready_for_cor %>%
  readr::write_tsv(file.path(res_path,"pattern_validation","ICP_mean_exp_in_GSE72056.tsv"))

broom::tidy(
  cor.test(ready_for_cor$Immune,ready_for_cor$Tumor,method = "spearman")
) %>%
  dplyr::mutate(y=(max(ready_for_cor$Immune)-min(ready_for_cor$Immune))*0.85+min(ready_for_cor$Immune),x=min(ready_for_cor$Tumor)+(max(ready_for_cor$Tumor)-min(ready_for_cor$Tumor))*0.4) %>%
  dplyr::mutate(label = paste("r = ",signif(estimate,2),", p = ",signif(p.value,2))) -> cor_anno

ready_for_cor %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
  tidyr::nest(-Exp_site) %>%
  dplyr::mutate(spm = purrr::map(data,.f=function(.x){
    if(nrow(.x)>5){
      broom::tidy(
        cor.test(.x$Immune,.x$Tumor,method = "spearman")) %>%
        dplyr::mutate(y=(max(.x$Immune)-min(.x$Immune))*0.85+min(.x$Immune),
                      x=2.72) %>%
        dplyr::mutate(label = paste("r = ",signif(estimate,2),", p = ",signif(p.value,2)))
    }else{
      tibble::tibble()
    }
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> cor_anno.specific
broom::tidy(
  cor.test(ready_for_cor$Immune,ready_for_cor$Tumor,method = "spearman")
) %>%
  dplyr::mutate(y=(max(ready_for_cor$Immune)-min(ready_for_cor$Immune))*0.75+min(ready_for_cor$Immune),x=min(ready_for_cor$Tumor)+(max(ready_for_cor$Tumor)-min(ready_for_cor$Tumor))*0.4) %>%
  dplyr::mutate(label = paste("r = ",signif(estimate,2),", p = ",signif(p.value,2))) -> cor_anno

ready_for_cor %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
  dplyr::filter(Exp_site!="Not_sure") %>%
  ggplot(aes(x=Tumor,y=Immune)) +
  geom_jitter(aes(color=Exp_site),size=0.5) +
  geom_smooth(se= F, method = "lm", aes(color=Exp_site,group=Exp_site)) +
  geom_smooth(se = F, method = "lm",color= "black") +
  geom_text(aes(x=x,y=y,label = label),data=cor_anno,color="black") +
  geom_text(aes(x=x,y=y,label = label,color=Exp_site),data = cor_anno.specific) +
  scale_color_manual(values=c("#CD661D",  "#008B00", "#FF69B4", "#1874CD","#CD3333")) +
  my_theme +
  xlab("Mean exppression in tumor cells") +
  ylab("Mean exppression in immune cells") 
ggsave(file.path(res_path,"pattern_validation","5.6.GSE72056-T-I-meanExp.correlation.pdf"),device = "pdf",height = 6,width = 8)
ggsave(file.path(res_path,"pattern_validation","5.6.GSE72056-T-I-meanExp.correlation.png"),device = "png",height = 6,width = 8)



