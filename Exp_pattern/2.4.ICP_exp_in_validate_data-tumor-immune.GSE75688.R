######################## verify the ICP expression site in melenoma single cell data set ############
######################## GSE75688
###### malignant tumor cell and non-malignant cell including: T cells, B cells, macrophages, CAFs(cancer-associated fibroblasts), and endothelial (Endo.) cells from preferentially expressed genes. NK, natural killer cells.

library(magrittr)

basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byMeanUQ")
# res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio.new")

# load image --------------------------------------------------------------
load(file.path(
  res_path,"pattern_validation","GSE75688.melenoma.TI.compare.Rdata")
)

#### gene list ------------------------------------------------------------------------
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list_exp_site <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(entrez_ID,symbol,Exp_site,`log2FC(I/T).mean`,`log2FC(I/T).mid`) %>%
  dplyr::inner_join(gene_list,by="symbol") 


# load expresion data -----------------------------------------------------

# sample info
sample_info.class <- readr::read_tsv(file.path(basic_path,"data/single_cell_RNAseq/GSE75688_breast_cancer","GSE75688_final_sample_information.txt")) %>%
  dplyr::filter(index2 != "Stromal") %>%
  dplyr::rename("cell_source"="index2")


# exp data
ICP_exp_in_GSE75688 <- readr::read_tsv(file.path(basic_path,"data/single_cell_RNAseq/GSE75688_breast_cancer","GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt"))

ICP_exp_in_GSE75688 %>%
  dplyr::select(-gene_type,-gene_id) %>%
  dplyr::filter(gene_name %in% gene_list_exp_site$symbol) -> ICP_exp_in_GSE75688


# Single cell data: compare ICP between tumor and immune cells ------------------------------
ICP_exp_in_GSE75688 %>%
  dplyr::rename("symbol"="gene_name") %>%
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  dplyr::filter(type == "SC") %>%
  dplyr::select(-index,-index3) -> ICP_exp_in_GSE75688.SC

ICP_exp_in_GSE75688 %>%
  dplyr::rename("symbol"="gene_name") %>%
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  dplyr::filter(type == "Bulk") %>%
  dplyr::select(-index,-index3) -> ICP_exp_in_GSE75688.Bulk

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

ICP_exp_in_GSE75688.SC %>%
  tidyr::nest(-symbol) %>%
  dplyr::mutate(test = purrr::map(data,fn_compare_TI_FC)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ICP_exp_in_GSE75688.wilcox.test.FC.TI.SC

ICP_exp_in_GSE75688.Bulk %>%
  tidyr::nest(-symbol) %>%
  dplyr::mutate(test = purrr::map(data,fn_compare_TI_FC)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ICP_exp_in_GSE75688.wilcox.test.FC.TI.Bulk

ICP_exp_in_GSE75688.wilcox.test.FC.TI.Bulk %>%
  dplyr::mutate(Exp_site = ifelse(`log2FC(I/T).mean` >=1 & `log2FC(I/T).UQ` >=2, "Immune cell dominate","Immune and tumor cell almost")) %>%
  dplyr::mutate(Exp_site = ifelse(`log2FC(I/T).mean` <=(-1) & `log2FC(I/T).UQ` <=(-2), "Tumor cell dominate",Exp_site)) %>%
  dplyr::inner_join(gene_list_exp_site, by="symbol") %>%
  dplyr::select(symbol, Exp_site.x, Exp_site.y) %>% 
  dplyr::filter(Exp_site.x==Exp_site.y)

ICP_exp_in_GSE75688.wilcox.test.FC.TI.SC %>%
  dplyr::mutate(Exp_site = ifelse(`log2FC(I/T).UQ` >=1, "Immune cell dominate","Immune and tumor cell almost")) %>%
  dplyr::mutate(Exp_site = ifelse(`log2FC(I/T).UQ` <=(-1), "Tumor cell dominate",Exp_site)) %>%
  dplyr::inner_join(gene_list_exp_site, by="symbol") %>%
  dplyr::select(symbol, Exp_site.x, Exp_site.y, `log2FC(I/T).UQ`) %>% 
  dplyr::mutate(fit_fantom=ifelse(Exp_site.x==Exp_site.y,"yes","no")) %>%
  readr::write_tsv(file.path(res_path,"predict_res_validate_by_GSE75688.tsv"))

strip_color <- data.frame(Exp_site = unique(gene_list_exp_site$Exp_site),
                          site_cplor = c("green", "orange", "pink"),
                          rank = c(3,2,1))
# plot
library(ggbeeswarm)
ICP_exp_in_GSE75688 %>%
  dplyr::rename("symbol"="gene_name") %>%
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  dplyr::filter(cell_source %in% c("Tumor","Immune")) %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol")  %>%
  dplyr::inner_join(ICP_exp_in_GSE75688.wilcox.test.FC.TI.SC,by="symbol")  -> ready_for_draw.SC

ready_for_draw.SC %>%
  dplyr::select(symbol,Exp_site,`log2FC(I/T).mean.x`) %>%
  dplyr::inner_join(strip_color,by="Exp_site") %>%
  dplyr::arrange(rank,`log2FC(I/T).mean.x`)-> symbol_rank.SC

ready_for_draw.SC <- within(ready_for_draw.SC,symbol <- factor(symbol,levels = unique(symbol_rank.SC$symbol)))  
with(ready_for_draw.SC, levels(symbol))

ready_for_draw.SC %>%
  dplyr::select(symbol,Exp_site) %>%
  unique() -> color_bac
color_bac$cell_source <- color_bac$Exp <- 1


ggplot(ready_for_draw.SC,aes(x=cell_source, y=log2(Exp+0.01))) +
  # geom_quasirandom(size=0.1) +
  geom_violin(size = 0.25) +
  geom_rect(data=color_bac,aes(fill = Exp_site),xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.1) +
  # geom_violin() +
  facet_wrap(~symbol,scales = "free_y", ncol = 7) +
  scale_fill_manual(
    name = "Expression pattern",
    values = c("yellow",  "green","pink"),
    # values = c("#008B00", "#00EE00", "#CD8500", "#FF4500"),
    breaks = c("Immune and tumor cell almost", "Immune cell dominate","Tumor cell dominate")
  ) +
  my_theme +
  labs(y="Expression",title="GSE75688, breast cancer") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 6),
    axis.text.y = element_text(size = 6),
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    legend.key.width = unit(0.1,"inches"),
    legend.key.height=unit(0.1,"inches"),
    strip.text = element_text(size = 6),
    plot.title = element_text(size = 12, face = "bold"),
    panel.grid=element_blank()
  )
ggsave(file.path(res_path,"pattern_validation","5.2.GSE75688.ICP_exp-T-I_compare-violin.pdf"),device = "pdf",height = 20,width = 16, units = c("cm"))
ggsave(file.path(res_path,"pattern_validation","5.2.GSE75688.ICP_exp-T-I_compare-violin.png"),device = "png",height = 20,width = 16, units = c("cm"))

# correlation between FC got from fantom and melanoma ---------------------------------

fantom_res <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(symbol,mean_cell_line, mean_immune_exp,mid_cell_line,mid_immune_exp,`log2FC(I/T).mean`,`log2FC(I/T).mid`) %>%
  dplyr::rename("log2FC(I/T).UQ"="log2FC(I/T).mid","UQ_tumor_exp"="mid_cell_line","UQ_immune_exp"="mid_immune_exp","mean_tumor_exp"="mean_cell_line") %>%
  tidyr::gather(-symbol,key="data_type",value="Fantom5") 

# single cells ------
ICP_exp_in_GSE75688.wilcox.test.FC.TI.SC %>%
  dplyr::select(symbol,mean_tumor_exp, mean_immune_exp,UQ_immune_exp, UQ_tumor_exp,`log2FC(I/T).mean`,`log2FC(I/T).UQ`) %>%
  tidyr::gather(-symbol,key="data_type",value="GSE75688") %>%
  dplyr::inner_join(fantom_res,by=c("symbol","data_type")) -> correlation.ready.SC

# spearman correlation
correlation.ready.SC %>%
  # dplyr::filter(symbol != "BTNL3") %>%
  tidyr::nest(-data_type) %>%
  dplyr::mutate(cpm_cor = purrr::map(data,.f=function(.x){
    broom::tidy(
      cor.test(.x$GSE75688,.x$Fantom5,method = "spearman")
    )
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> cor.res.SC

# plot
library(ggplot2)
library(ggpubr)

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

correlation.ready.SC %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>% 
  dplyr::group_by(data_type) %>%
  dplyr::mutate(y=(max(Fantom5)-min(Fantom5))*0.85+min(Fantom5),x=min(GSE75688)+(max(GSE75688)-min(GSE75688))*0.4) %>%
  dplyr::select(data_type,x,y) %>%
  unique() %>%
  dplyr::inner_join(cor.res.SC,by="data_type") %>%
  dplyr::select(data_type,x,y,estimate,p.value) %>%
  dplyr::mutate(label = paste("r = ",signif(estimate,2),", p = ",signif(p.value,2),sep="")) -> cor_text.SC

library(latex2exp)
correlation.ready.SC %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
  dplyr::filter(Exp_site!="Not_sure", data_type %in% c("log2FC(I/T).mean","log2FC(I/T).UQ")) %>%
  ggplot(aes(x=GSE75688,y=Fantom5)) +
  geom_jitter(aes(color=Exp_site)) +
  geom_smooth(se = F, method = "lm") +
  geom_text(aes(x=x,y=y,label = label),data=cor_text.SC %>% dplyr::filter(data_type%in% c("log2FC(I/T).mean","log2FC(I/T).UQ"))) +
  facet_wrap(~data_type,scales = "free") +
  scale_color_manual(values=c("#CD661D",  "#008B00", "#FF69B4", "#1874CD","#CD3333")) +
  my_theme +
  labs(x="GSE75688",
       y="FANTOM5",
       title = "FANTOM5 vs. GSE75688") +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.2,"inches"),
    legend.key.height=unit(0.2,"inches"),
    legend.text = element_text(size=8),
    legend.title = element_blank()
  ) 
ggsave(file.path(res_path,"pattern_validation","8.1.GSE75688-Fantom5.correlation.pdf"),device = "pdf",height = 4,width = 6)
ggsave(file.path(res_path,"pattern_validation","8.1.GSE75688-Fantom5.correlation.png"),device = "png",height = 4,width = 6)

# bulk  ------
ICP_exp_in_GSE75688.wilcox.test.FC.TI.Bulk %>%
  dplyr::select(symbol,mean_tumor_exp, mean_immune_exp,UQ_immune_exp, UQ_tumor_exp,`log2FC(I/T).mean`,`log2FC(I/T).UQ`) %>%
  tidyr::gather(-symbol,key="data_type",value="GSE75688") %>%
  dplyr::inner_join(fantom_res,by=c("symbol","data_type")) -> correlation.ready.Bulk

# spearman correlation
correlation.ready.Bulk %>%
  # dplyr::filter(symbol != "BTNL3") %>%
  tidyr::nest(-data_type) %>%
  dplyr::mutate(cpm_cor = purrr::map(data,.f=function(.x){
    broom::tidy(
      cor.test(.x$GSE75688,.x$Fantom5,method = "spearman")
    )
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> cor.res.Bulk

# plot

correlation.ready.Bulk %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>% 
  dplyr::group_by(data_type) %>%
  dplyr::mutate(y=(max(Fantom5)-min(Fantom5))*0.85+min(Fantom5),x=min(GSE75688)+(max(GSE75688)-min(GSE75688))*0.4) %>%
  dplyr::select(data_type,x,y) %>%
  unique() %>%
  dplyr::inner_join(cor.res.Bulk,by="data_type") %>%
  dplyr::select(data_type,x,y,estimate,p.value) %>%
  dplyr::mutate(label = paste("r = ",signif(estimate,2),", p = ",signif(p.value,2),sep="")) -> cor_text.Bulk

correlation.ready.Bulk %>%
  dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
  dplyr::filter(Exp_site!="Not_sure", data_type %in% c("log2FC(I/T).mean","log2FC(I/T).UQ")) %>%
  ggplot(aes(x=GSE75688,y=Fantom5)) +
  geom_jitter(aes(color=Exp_site)) +
  geom_smooth(se = F, method = "lm") +
  geom_text(aes(x=x,y=y,label = label),data=cor_text.Bulk %>% dplyr::filter(data_type%in% c("log2FC(I/T).mean","log2FC(I/T).UQ"))) +
  facet_wrap(~data_type,scales = "free") +
  scale_color_manual(values=c("#CD661D",  "#008B00", "#FF69B4", "#1874CD","#CD3333")) +
  my_theme +
  labs(x="GSE75688",
       y="FANTOM5",
       title = "FANTOM5 vs. GSE75688") +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.2,"inches"),
    legend.key.height=unit(0.2,"inches"),
    legend.text = element_text(size=8),
    legend.title = element_blank()
  ) 
ggsave(file.path(res_path,"pattern_validation","5.1.GSE75688-Fantom5.correlation.pdf"),device = "pdf",height = 4,width = 6)
ggsave(file.path(res_path,"pattern_validation","5.1.GSE75688-Fantom5.correlation.png"),device = "png",height = 4,width = 6)




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

ICP_exp_in_GSE75688 %>%
  tidyr::gather(-symbol,key="sample",value="Exp") %>%
  dplyr::arrange(symbol) %>%
  dplyr::filter(!sample %in% duplicated_samples$sample) %>%
  dplyr::inner_join(sample_info.class,by="sample") %>%
  dplyr::filter(cell_source %in% c("Tumor","Immune")) %>%
  fn_plot_ICP_exp_in_dataset(ylab="Expression",facet="~cell_source",title="GSE75688, ICP expression",filename="5.5.GSE75688.ICP_exp.T-I")


# correlation between ICP exp in tumor and immune -------------------------
ICP_exp_in_GSE75688 %>%
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
  readr::write_tsv(file.path(res_path,"pattern_validation","ICP_mean_exp_in_GSE75688.tsv"))

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
ggsave(file.path(res_path,"pattern_validation","5.6.GSE75688-T-I-meanExp.correlation.pdf"),device = "pdf",height = 6,width = 8)
ggsave(file.path(res_path,"pattern_validation","5.6.GSE75688-T-I-meanExp.correlation.png"),device = "png",height = 6,width = 8)

ICP_exp_in_GSE75688.wilcox.test.FC.TI.SC %>%
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
       title = "Classification of ICPs' expression pattern, GSE75688") +
  scale_color_manual(values = c("#CD950C", "#66CD00", "#EE2C2C"),
                     name = "ICPs expression pattern") +
  my_theme +
  theme(
    plot.title = element_text(size=15)
  )
ggsave(file.path(res_path,"classify_ICP_exp_pattern_onlybyUQ.GSE75688.pdf"),device = "pdf",height = 4, width = 8)
ggsave(file.path(res_path,"classify_ICP_exp_pattern_onlybyUQ.GSE75688.png"),device = "png",height = 4, width = 8)
# save image --------------------------------------------------------------
save.image(file.path(
  res_path,"pattern_validation","GSE75688.melenoma.TI.compare.Rdata")
)
