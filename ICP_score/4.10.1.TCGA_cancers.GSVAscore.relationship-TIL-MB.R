########################## GSVA score of features are: expression site of ICPs, gene family
########################## get the correlation between GSVA score and TIL, mutation burden
library(magrittr)
library(tidyverse)
library(survival)

# data path ---------------------------------------------------------------
# server 1
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path(basic_path,"immune_checkpoint/data/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"TCGA_GSVAScore")
tcga_path <- file.path(basic_path,"/data/TCGA")

# load image --------------------------------------------------------------
load(file.path(res_path,"TCGA.GSVA_score.ICP_features.rda"))

# laod ICP feature info ---------------------------------------------------

gene_list <- readr::read_tsv(file.path(gene_list_path, "ICPs_all_info_class.tsv")) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Only_exp_on_Immune","Mainly_exp_on_Immune"),"Mainly_exp_on_Immune",Exp_site)) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Only_exp_on_Tumor","Mainly_exp_on_Tumor" ),"Mainly_exp_on_Tumor",Exp_site.1)) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("N"),"Not_sure",Exp_site.1))

# load GSVA score ---------------------------------------------------------

GSVA.score.onlytumor <- 
  readr::read_rds(file.path(res_path,"TCGA_cancer_specific.onlyTumor_GSVA.score_ICPs_features.rds.gz"))

# plot theme --------------------------------------------------------------

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
# 1. GSVA score correlation with TIL -----------------------------------------
# 1.1. load TIL ----
immunity_path_2 <- "/project/huff/huff/immune_checkpoint/data/immunity"
TIMER_immunity <- readr::read_tsv(file.path(immunity_path_2,"immuneEstimation.txt")) %>%
  dplyr::mutate(TIL = B_cell+CD4_Tcell+CD8_Tcell+Neutrophil+Macrophage+Dendritic)

# 1.2.function --------
fn_correlation_and_DE <- function(GSVA,TIL){
  GSVA %>%
    dplyr::mutate(barcode = substr(barcode,1,15)) %>%
    tidyr::gather(-barcode,key="Features",value="GSVA_score")%>%
    dplyr::inner_join(TIL,by="barcode") %>%
    dplyr::group_by(Features) %>%
    dplyr::mutate(group = ifelse(GSVA_score >= quantile(GSVA_score,0.5),"GSVA_high","GSVA_low")) %>%
    tidyr::nest(-Features) -> .tmp
  
  # correlation
  .tmp %>%
    dplyr::mutate(cor = purrr::map(data,.f=function(.x){
      cor.test(.x$TIL,.x$GSVA_score,method = "spearman") %>%
        broom::tidy()
    })) %>%
    dplyr::select(-data) -> cor.res
  
  # DE
  .tmp %>%
    dplyr::mutate(DE = purrr::map2(data,Features,.f=function(.x,.y){
      print(.y)
      .x %>%
        dplyr::filter(group == "GSVA_high") %>%
        .$GSVA_score %>%
        mean() -> mean_high
      .x %>%
        dplyr::filter(group == "GSVA_low") %>%
        .$GSVA_score %>%
        mean() -> mean_low
      broom::tidy(wilcox.test(TIL  ~ group, data = .x, alternative = "two.sided")) %>%
        dplyr::mutate(mean_high = mean_high, mean_low = mean_low) %>%
        dplyr::mutate(`log2FC(High/Low)` = log2((mean_high+1)/(mean_low+1))) 
    })) %>%
    dplyr::select(-data)  -> DE.res
  
  cor.res %>%
    dplyr::inner_join(DE.res,by="Features")
}

# 1.3.calculation ------
GSVA.score.onlytumor %>% 
  # head(1) %>%
  dplyr::mutate(res = purrr::map(GSVA,.f=fn_correlation_and_DE,
                                 TIL = TIMER_immunity %>% 
                                   dplyr::select(barcode,TIL))) %>%
  dplyr::select(-GSVA) %>%
  tidyr::unnest() -> GSVA.TIL.res

# 1.4.ploting ---------
# 1.4.1.TIL correlation by cancer types --------
set.seed(123)
GSVA.TIL.res %>%
  dplyr::select(-DE) %>%
  tidyr::unnest() %>%
  dplyr::mutate(`-log10P` = -log10(p.value)) %>%
  dplyr::mutate(`-log10P` = purrr::map(`-log10P`,.f=function(.x){
    ifelse(.x>10,10+runif(1,0,1),.x)
  })) %>%
  tidyr::unnest() %>%
  dplyr::rename("Correlation"="estimate") -> GSVA.TIL.cor.plot

readr::read_tsv(file.path(tcga_path,"02_pcc.tsv")) %>%
  dplyr::inner_join(GSVA.TIL.cor.plot,by="cancer_types") %>%
  dplyr::select(cancer_types,color) %>%
  unique() -> cancer.color

tibble::tibble(x=c(0.3,1,1,0.3,-0.3,-1,-1,-0.3),
               y=c(1.3,1.3,11,11,1.3,1.3,11,11),
               xend=c(1,1,0.3,0.3,-1,-1,-0.3,-0.3),
               yend=c(1.3,11,11,1.3,1.3,11,11,1.3)) %>%
  dplyr::mutate(sig = ifelse(x<0,"neg_sig","pos_sig")) -> frame.coord

GSVA.TIL.cor.plot %>%
  dplyr::mutate(sig = ifelse(p.value<=0.05 & Correlation>=0.3,"pos_sig","not")) %>%
  dplyr::mutate(sig = ifelse(p.value<=0.05 & Correlation<=(-0.3),"neg_sig",sig)) %>%
  dplyr::group_by(sig) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n_all=n()) %>%
  dplyr::mutate(percentage = paste(signif((n/n_all)*100,2),"%",sep="")) %>%
  dplyr::select(sig,n, n_all,percentage) %>%
  unique() %>%
  dplyr::mutate(x=c(0.8,0.2,-0.8),y=c(3,3,3)) %>%
  dplyr::filter(sig!="not")-> percentage.label

GSVA.TIL.cor.plot %>%
  ggplot(aes(x=Correlation,y=`-log10P`)) +
  geom_jitter(aes(color=cancer_types)) +
  scale_color_manual(
    values = cancer.color$color,
    name = "Cancer types"
  ) +
  labs(x="Spearman correlation between TIL and GSVA score",
       y=latex2exp::TeX("-log_{10} (P value)")) +
  my_theme +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=frame.coord) +
  geom_text(aes(x=x,y=y,label=percentage),data = percentage.label) 
ggsave(file.path(res_path,"1.TIL_with_GSVA_score","1.spm.cor.TIL-GSVAscore.allcancers.pdf"),device = "pdf",height = 6, width=8)
ggsave(file.path(res_path,"1.TIL_with_GSVA_score","1.spm.cor.TIL-GSVAscore.allcancers.png"),device = "png",height = 6, width=8)

# 1.4.2.TIL DE cancer types  --------
GSVA.TIL.res %>%
  dplyr::select(-cor) %>%
  tidyr::unnest() %>%
  dplyr::mutate(`-log10P` = -log10(p.value)) %>%
  dplyr::mutate(`-log10P` = purrr::map(`-log10P`,.f=function(.x){
    ifelse(.x>10,10+runif(1,0,1),.x)
  })) %>%
  tidyr::unnest()  -> GSVA.TIL.DE.plot

GSVA.TIL.DE.plot %>%
  ggplot(aes(x=`log2FC(High/Low)`,y=`-log10P`)) +
  geom_jitter(aes(color=cancer_types)) +
  scale_color_manual(
    values = cancer.color$color,
    name = "Cancer types"
  ) +
  labs(x="log2(TIL fold change between high and low \n GSVA score groups)",
       y="-log10(P value)") +
  my_theme
ggsave(file.path(res_path,"1.TIL_with_GSVA_score","2.DE.TIL-GSVAscore.allcancers.pdf"),device = "pdf",height = 6, width=8)
ggsave(file.path(res_path,"1.TIL_with_GSVA_score","2.DE.TIL-GSVAscore.allcancers.png"),device = "png",height = 6, width=8)

# 1.4.3.mean TIL DE in cancer types  --------
GSVA.TIL.res %>%
  tidyr::unnest() %>%
  dplyr::group_by(Features) %>%
  dplyr::mutate(mean_cor = mean(estimate),mean_FC = mean(`log2FC(High/Low)`)) %>%
  dplyr::select(Features,mean_FC,mean_cor) %>%
  dplyr::ungroup() %>%
  unique() -> GSVA.TIL.cor.DE.plot

gene_list %>%
  dplyr::mutate(All_gene = "All_gene") %>%
  dplyr::select(symbol,functionWithImmune,family,Recepter_pairs,Exp_site.1,All_gene) %>%
  tidyr::gather(-symbol,key="Groups",value="Features") %>%
  tidyr::nest(-Features) %>%
  dplyr::filter(! Features %in% "Other") %>%
  dplyr::filter(!is.na(Features)) %>%
  dplyr::mutate(Feature_symbol = purrr::map(data,.f=function(.x){
    tibble::tibble(
      Feature_symbol=paste(.x$symbol,collapse="/"),
      n=length(.x$symbol))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(label = ifelse(n<5,Feature_symbol,Features)) %>%
  dplyr::inner_join(GSVA.TIL.cor.DE.plot,by="Features") -> label.data

library(ggrepel)

GSVA.TIL.cor.DE.plot %>%
  ggplot(aes(x=mean_cor,y=mean_FC)) +
  geom_jitter(aes(color=Features)) +
  geom_text_repel(aes(x=mean_cor,y=mean_FC,label=label,color=Features),data=label.data,size=3) +
  scale_color_manual(
    values = c(cancer.color$color[1:nrow(GSVA.TIL.cor.DE.plot)],cancer.color$color[1:nrow(GSVA.TIL.cor.DE.plot)])
  ) +
  labs(x="Average spearman correlation between \nTIL and GSVA score in tumors",
       y="Average log2(TIL fold change between high\nand low  GSVA score groups) in tumors") +
  my_theme +
  theme(legend.position = "none")


ggsave(file.path(res_path,"1.TIL_with_GSVA_score","3.DE.Cor.TIL-GSVAscore.allcancers-Mean.pdf"),device = "pdf",height = 4, width=5)
ggsave(file.path(res_path,"1.TIL_with_GSVA_score","3.DE.Cor.TIL-GSVAscore.allcancers-Mean.png"),device = "png",height = 4, width=5)

# 1.4.4.all TIL DE in cancer types  --------
GSVA.TIL.res %>%
  tidyr::unnest() %>%
  dplyr::rename("Correlation"="estimate") -> GSVA.TIL.cor.DE.plot.all

tibble::tibble(x=c(0.3,1,1,0.3,-0.3,-1,-1,-0.3),
               y=c(1,1,3,3,1,1,3,3),
               xend=c(1,1,0.3,0.3,-1,-1,-0.3,-0.3),
               yend=c(1,3,3,1,1,3,3,1)) %>%
  dplyr::mutate(sig = ifelse(x<0,"neg_sig","pos_sig")) -> frame.coord.all

GSVA.TIL.cor.DE.plot.all %>%
  dplyr::mutate(sig = ifelse(`log2FC(High/Low)`>=1 & Correlation>=0.3,"pos_sig","not")) %>%
  dplyr::mutate(sig = ifelse(`log2FC(High/Low)`>=1 & Correlation<=(-0.3),"neg_sig",sig)) %>%
  dplyr::group_by(sig) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n_all=n()) %>%
  dplyr::mutate(percentage = paste(signif((n/n_all)*100,2),"%",sep="")) %>%
  dplyr::select(sig,n, n_all,percentage) %>%
  unique() %>%
  dplyr::mutate(x=c(0.8,0.2,-0.8),y=c(2.8,2.7,2.8)) %>%
  dplyr::filter(sig!="not")-> percentage.label.all

GSVA.TIL.cor.DE.plot.all %>%
  ggplot(aes(x=Correlation,y=`log2FC(High/Low)`)) +
  geom_jitter(aes(color=Features)) +
  scale_color_manual(
    values = c(cancer.color$color[1:nrow(GSVA.TIL.cor.DE.plot)],cancer.color$color[1:nrow(GSVA.TIL.cor.DE.plot)])
  ) +
  labs(x="Spearman correlation between \nTIL and GSVA score in tumors",
       y=latex2exp::TeX("log_{2} (TIL fold change between high\nand low  GSVA score groups) in tumors")) +
  my_theme +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,1),"cm"))+
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=frame.coord.all) +
  geom_text(aes(x=x,y=y,label=percentage),data = percentage.label.all) 


ggsave(file.path(res_path,"1.TIL_with_GSVA_score","4.DE.Cor.TIL-GSVAscore.allcancers.pdf"),device = "pdf",height = 4, width=5)
ggsave(file.path(res_path,"1.TIL_with_GSVA_score","4.DE.Cor.TIL-GSVAscore.allcancers.png"),device = "png",height = 4, width=5)

# 2. GSVA score with mutation burden ----------------------------
# 2.1. load mutation data ------
TCGA_mutation_burden <- readr::read_rds(file.path("/project/huff/huff/data/TCGA/snv","pancan33_maf_syn7824274.sample.mut_burden.rds.gz"))

# 2.2. calculation --------
GSVA.score.onlytumor %>% 
  # head(1) %>%
  dplyr::mutate(res = purrr::map(GSVA,.f=fn_correlation_and_DE,
                                 TIL = TCGA_mutation_burden %>% 
                                   dplyr::select(sample,mutation_durden) %>%
                                   dplyr::rename("barcode"="sample","TIL"="mutation_durden"))) %>%
  dplyr::select(-GSVA) %>%
  tidyr::unnest() -> GSVA.MB.res

# 2.3.ploting ---------
# 2.3.1.Mutation Burden correlation by cancer types --------
GSVA.MB.res %>%
  dplyr::select(-DE) %>%
  tidyr::unnest() %>%
  dplyr::mutate(`-log10P` = -log10(p.value)) %>%
  dplyr::mutate(`-log10P` = purrr::map(`-log10P`,.f=function(.x){
    ifelse(.x>10,10+runif(1,0,1),.x)
  })) %>%
  tidyr::unnest() %>%
  dplyr::rename("Correlation"="estimate") -> GSVA.MB.cor.plot

readr::read_tsv(file.path(tcga_path,"02_pcc.tsv")) %>%
  dplyr::inner_join(GSVA.MB.res,by="cancer_types") %>%
  dplyr::select(cancer_types,color) %>%
  unique() -> cancer.color.MB

gene_list %>%
  dplyr::mutate(All_gene = "All_gene") %>%
  dplyr::select(symbol,functionWithImmune,family,Recepter_pairs,Exp_site.1,All_gene) %>%
  tidyr::gather(-symbol,key="Groups",value="Features") %>%
  tidyr::nest(-Features) %>%
  dplyr::filter(! Features %in% "Other") %>%
  dplyr::filter(!is.na(Features)) %>%
  dplyr::mutate(Feature_symbol = purrr::map(data,.f=function(.x){
    tibble::tibble(
      Feature_symbol=paste(.x$symbol,collapse="/"),
      n=length(.x$symbol))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(label = ifelse(n<5,Feature_symbol,Features)) %>%
  dplyr::inner_join(GSVA.MB.cor.plot,by="Features") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(Max_min = purrr::map(data,.f=function(.x){
    .max <- max(.x$Correlation)
    .min <- min(.x$Correlation)
    .x %>%
      dplyr::mutate(Max_min=ifelse(Correlation==.max,"yes","no")) %>%
      dplyr::mutate(Max_min=ifelse(Correlation==.min,"yes",Max_min)) 
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(label = ifelse(Max_min=="yes",label,NA)) %>%
  dplyr::ungroup() -> label.data.MB.cor

tibble::tibble(x=c(0.3,1,1,0.3,-0.3,-1,-1,-0.3),
               y=c(1.3,1.3,11,11,1.3,1.3,11,11),
               xend=c(1,1,0.3,0.3,-1,-1,-0.3,-0.3),
               yend=c(1.3,11,11,1.3,1.3,11,11,1.3)) %>%
  dplyr::mutate(sig = ifelse(x<0,"neg_sig","pos_sig")) -> frame.coord.cor.MB

GSVA.MB.cor.plot %>%
  dplyr::mutate(sig = ifelse(p.value<=0.05 & Correlation>=0.3,"pos_sig","not")) %>%
  dplyr::mutate(sig = ifelse(p.value<=0.05 & Correlation<=(-0.3),"neg_sig",sig)) %>%
  dplyr::group_by(sig) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n_all=n()) %>%
  dplyr::mutate(percentage = paste(signif((n/n_all)*100,2),"%",sep="")) %>%
  dplyr::select(sig,n, n_all,percentage) %>%
  unique() %>%
  dplyr::mutate(x=c(0.8,0.8,-0.8),y=c(3,3,3)) %>%
  dplyr::filter(sig!="not") -> percentage.label.cor.MB

GSVA.MB.cor.plot %>%
  ggplot(aes(x=Correlation,y=`-log10P`)) +
  geom_jitter(aes(color=cancer_types)) +
  scale_color_manual(
    values = cancer.color.MB$color,
    name = "Cancer types"
  ) +
  # geom_text_repel(aes(x=Correlation,y=`-log10P`,label=label,color=cancer_types),
  #                 data=label.data.MB.cor, size = 2) +
  labs(x="Spearman correlation between \nmutation burden and GSVA score",
       y="-log10(P value)") +
  my_theme +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=frame.coord.cor.MB) +
  geom_text(aes(x=x,y=y,label=percentage),data = percentage.label.cor.MB) 
  
ggsave(file.path(res_path,"2.MB_with_GSVA_score","1.spm.cor.MB-GSVAscore.allcancers.pdf"),device = "pdf",height = 6, width=8)
ggsave(file.path(res_path,"2.MB_with_GSVA_score","1.spm.cor.MB-GSVAscore.allcancers.png"),device = "png",height = 6, width=8)

# 2.3.2.MB DE cancer types  --------
GSVA.MB.res %>%
  dplyr::select(-cor) %>%
  tidyr::unnest() %>%
  dplyr::mutate(`-log10P` = -log10(p.value)) %>%
  dplyr::mutate(`-log10P` = purrr::map(`-log10P`,.f=function(.x){
    ifelse(.x>10,10+runif(1,0,1),.x)
  })) %>%
  tidyr::unnest()  -> GSVA.MB.DE.plot

GSVA.MB.DE.plot %>%
  ggplot(aes(x=`log2FC(High/Low)`,y=`-log10P`)) +
  geom_jitter(aes(color=cancer_types)) +
  scale_color_manual(
    values = cancer.color$color,
    name = "Cancer types"
  ) +
  labs(x="log2(Mutation burden fold change between\nhigh and low GSVA score groups)",
       y="-log10(P value)") +
  my_theme
ggsave(file.path(res_path,"2.MB_with_GSVA_score","2.DE.MB-GSVAscore.allcancers.pdf"),device = "pdf",height = 6, width=8)
ggsave(file.path(res_path,"2.MB_with_GSVA_score","2.DE.MB-GSVAscore.allcancers.png"),device = "png",height = 6, width=8)

# 2.3.3.mean MB DE and Cor. in cancer types  --------
GSVA.MB.res %>%
  tidyr::unnest() %>%
  dplyr::group_by(Features) %>%
  dplyr::mutate(mean_cor = mean(estimate),mean_FC = mean(`log2FC(High/Low)`)) %>%
  dplyr::select(Features,mean_FC,mean_cor) %>%
  dplyr::ungroup() %>%
  unique() -> GSVA.MB.cor.DE.plot

gene_list %>%
  dplyr::mutate(All_gene = "All_gene") %>%
  dplyr::select(symbol,functionWithImmune,family,Recepter_pairs,Exp_site.1,All_gene) %>%
  tidyr::gather(-symbol,key="Groups",value="Features") %>%
  tidyr::nest(-Features) %>%
  dplyr::filter(! Features %in% "Other") %>%
  dplyr::filter(!is.na(Features)) %>%
  dplyr::mutate(Feature_symbol = purrr::map(data,.f=function(.x){
    tibble::tibble(
      Feature_symbol=paste(.x$symbol,collapse="/"),
      n=length(.x$symbol))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(label = ifelse(n<5,Feature_symbol,Features)) %>%
  dplyr::inner_join(GSVA.MB.cor.DE.plot,by="Features") -> label.data.meanMB

library(ggrepel)

GSVA.MB.cor.DE.plot %>%
  ggplot(aes(x=mean_cor,y=mean_FC)) +
  geom_jitter(aes(color=Features)) +
  geom_text_repel(aes(x=mean_cor,y=mean_FC,label=label,color=Features),data=label.data.meanMB,size=3) +
  scale_color_manual(
    values = c(cancer.color$color[1:nrow(GSVA.MB.cor.DE.plot)],cancer.color$color[1:nrow(GSVA.MB.cor.DE.plot)])
  ) +
  labs(x="Average spearman correlation between\n
       mutation burden and GSVA score in tumors",
       y="Average log2(Mutation burden fold change between high\n
       and low  GSVA score groups) in tumors") +
  my_theme +
  theme(legend.position = "none")

ggsave(file.path(res_path,"2.MB_with_GSVA_score","3.DE.Cor.MB-GSVAscore.allcancers-Mean.pdf"),device = "pdf",height = 4, width=5)
ggsave(file.path(res_path,"2.MB_with_GSVA_score","3.DE.Cor.MB-GSVAscore.allcancers-Mean.png"),device = "png",height = 4, width=5)


# 2.3.4.all TIL DE in cancer types  --------
GSVA.MB.res %>%
  tidyr::unnest() %>%
  dplyr::rename("Correlation"="estimate") -> GSVA.MB.DE.plot.all

tibble::tibble(x=c(0.3,1,1,0.3,-0.3,-1,-1,-0.3),
               y=c(1,1,3,3,1,1,3,3),
               xend=c(1,1,0.3,0.3,-1,-1,-0.3,-0.3),
               yend=c(1,3,3,1,1,3,3,1)) %>%
  dplyr::mutate(sig = ifelse(x<0,"neg_sig","pos_sig")) -> frame.coord.MB.all

GSVA.MB.DE.plot.all %>%
  dplyr::mutate(sig = ifelse(`log2FC(High/Low)`>=1 & Correlation>=0.3,"pos_sig","not")) %>%
  dplyr::mutate(sig = ifelse(`log2FC(High/Low)`>=1 & Correlation<=(-0.3),"neg_sig",sig)) %>%
  dplyr::group_by(sig) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n_all=n()) %>%
  dplyr::mutate(percentage = paste(signif((n/n_all)*100,2),"%",sep="")) %>%
  dplyr::select(sig,n, n_all,percentage) %>%
  unique() %>%
  dplyr::mutate(x=c(0.8,0.8,-0.8),y=c(2.8,2.7,2.8)) %>%
  dplyr::filter(sig!="not") -> percentage.label.MB.all

GSVA.MB.DE.plot.all %>%
  ggplot(aes(x=Correlation,y=`log2FC(High/Low)`)) +
  geom_jitter(aes(color=Features)) +
  scale_color_manual(
    values = c(cancer.color$color[1:nrow(GSVA.MB.DE.plot)],cancer.color$color[1:nrow(GSVA.MB.DE.plot)])
  ) +
  labs(x="Spearman correlation between \nmutation burden and GSVA score in tumors",
       y=latex2exp::TeX("log_{2} (mutation burden fold change between \n high and low GSVA score groups) in tumors")) +
  my_theme +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,1),"cm"))+
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=frame.coord.MB.all) +
  geom_text(aes(x=x,y=y,label=percentage),data = percentage.label.MB.all) 


ggsave(file.path(res_path,"1.TIL_with_GSVA_score","4.DE.Cor.TIL-GSVAscore.allcancers.pdf"),device = "pdf",height = 4, width=5)
ggsave(file.path(res_path,"1.TIL_with_GSVA_score","4.DE.Cor.TIL-GSVAscore.allcancers.png"),device = "png",height = 4, width=5)

# 3.survival analysis -------------------------------------------------------
# 3.1.load survival data ----
clinical_tcga <- readr::read_rds(file.path(basic_path,"TCGA_survival/data","Pancan.Merge.clinical.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::select(-cancer_types) %>%
  unique() %>%
  dplyr::mutate(OS=as.numeric(OS),Status=as.numeric(Status),Age=as.numeric(Age)) %>%
  dplyr::group_by(barcode) %>%
  dplyr::mutate(OS= max(OS)) %>%
  dplyr::mutate(Status =  max(Status)) %>%
  dplyr::ungroup()

clinical <- readr::read_rds(file.path(basic_path,"data/TCGA-survival-time/cell.2018.survival","TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::select(type,bcr_patient_barcode,PFS,PFS.time) %>%
  tidyr::drop_na() %>%
  dplyr::rename("barcode" = "bcr_patient_barcode") %>%
  unique()

clinical %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::rename("cancer_types"= "type") -> clinical_all_data

# 3.2.function to do cox survival analysis  ----
# 3.2.1.univariable cox analysis ---------
fn_survival_test <- function(data,feature){
  print(feature)
  .cox <- survival::coxph(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)
  summary(.cox) -> .z
  
  # KM pvalue
  kmp <- 1 - pchisq(survival::survdiff(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)$chisq,df = length(levels(as.factor(data$group))) - 1)
  
  # mearge results
  tibble::tibble(
    group = rownames(.z$conf.int),
    n = .z$n,
    hr = .z$conf.int[1],
    hr_l = .z$conf.int[3],
    hr_h = .z$conf.int[4],
    coxp = .z$waldtest[3],
    kmp = kmp) %>%
    dplyr::mutate(status = ifelse(hr > 1, "High_risk", "Low_risk"))
  # multi-variable cox analysis
}

# 3.2.2.multi-variable cox analysis ---------
fn_survival_test.multiCox <- function(data,uni_sig_feature){
  covariates <- uni_sig_feature$Features %>% unique()
  if(length(covariates)>=1){
    multi_formulas <- as.formula(paste('Surv(time, status)~', paste(covariates,collapse = " + ")))
    model <- coxph(multi_formulas, data = data)
    model.s <- summary(model)
    
    # Extract data
    tibble::tibble(
      Features = rownames(model.s$coefficients),
      n = model.s$n,
      hr = model.s$conf.int[,1],
      hr_l = model.s$conf.int[,3],
      hr_h = model.s$conf.int[,4],
      coxp = model.s$coefficients[,5]) %>%
      dplyr::mutate(status = ifelse(hr > 1, "High_risk", "Low_risk"))
  }else{
    tibble::tibble()
  }
}

# 3.2.3.drawing cox plot ------
fn_cox_plot <- function(data,filename,title,facet, dir,w=4,h=4){
  data %>% 
    # dplyr::mutate(functionWithImmune=functionWithImmune) %>%
    dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
    ggplot(aes(y = hr, x = cancer_types, ymin=hr_l,ymax=hr_h)) +
    geom_pointrange(aes(color=cox_sig),size=0.5) +
    scale_color_manual(values=c("red","black")) +
    geom_hline(aes(yintercept = 1), linetype =2) +
    scale_size(name = "p-value") +
    scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),
                       labels = c("1/16","1/8","1/4","1/2",1,2,3)) +
    facet_grid(as.formula(facet),scales = "free", space = "free") +
    # facet_wrap(as.formula(facet)) +
    coord_flip() +
    ggthemes::theme_gdocs() +
    my_theme +
    theme(
      legend.position = "none",
      axis.line.y = element_line(color="black"),
      axis.text = element_text(color = "black",size=8),
      axis.title = element_text(color = "black",size=10),
      text = element_text(color = "black")
    ) +
    labs(y = "Hazard Ratio (High vs. low expression)", x = "Cancers",title = title) -> p
  ggsave(file.path(out_path,"e_6_exp_profile",dir,paste(filename,"png",sep=".")),device = "png",width = w,height = h)
  ggsave(file.path(out_path,"e_6_exp_profile",dir,paste(filename,"pdf",sep=".")),device = "pdf",width = w,height = h)
}
# 3.3. univariate survival analysis ------
# 3.3.1.PFS ----
stage_class <- tibble::tibble(Stage = c("i/ii nos","is","Not_applicable",
                                        "stage i","stage ia","stage ib",
                                        "stage ii","stage iia","stage iib","stage iic",
                                        "stage iii","stage iiia","stage iiib","stage iiic",
                                        "stage iv","stage iva","stage ivb","stage ivc",
                                        "stage x",NA),
                              group = c("Not_applicable","Not_applicable","Not_applicable",
                                        "stage i","stage i","stage i",
                                        "stage ii","stage ii","stage ii","stage ii",
                                        "stage iii","stage iii","stage iii","stage iii",
                                        "stage iv","stage iv","stage iv","stage iv",
                                        "Not_applicable","Not_applicable"))
GSVA.score.onlytumor %>%
  dplyr::mutate(clinical_gsva = purrr::map2(GSVA,cancer_types,.f=function(.x,.y){
    print(.y)
    .x %>%
      dplyr::mutate(barcode = substr(barcode,1,12)) %>%
      dplyr::inner_join(clinical_all_data,by="barcode")  %>%
      dplyr::rename("status"="PFS","time"="PFS.time") %>%
      dplyr::select(barcode,Stage,time,status) %>%
      dplyr::inner_join(stage_class,by="Stage") %>%
      dplyr::filter(group != "Not_applicable") %>%
      tidyr::gather(-barcode,-status,-time,-group,key="Features",value="value") %>%
      dplyr::select(barcode,status,time,Features,value,group) -> stage.surv.data
    .x %>%
      dplyr::mutate(barcode = substr(barcode,1,12)) %>%
      dplyr::inner_join(clinical_all_data,by="barcode") %>%
      dplyr::rename("status"="PFS","time"="PFS.time") %>%
      dplyr::select(-OS,-Status,-cancer_types,-Stage) %>%
      tidyr::gather(-barcode,-status,-time,key="Features",value="value")  %>%
      dplyr::filter(!is.na(value)) %>%
      dplyr::group_by(Features) %>%
      dplyr::mutate(group = ifelse(value >= quantile(value,0.5),"2high","1low")) %>%
      dplyr::mutate(n=n()) %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(nn=n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n>10,nn>=2) %>%
      dplyr::select(-n,-nn) %>%
      rbind(stage.surv.data) %>%
      tidyr::nest(-Features) 
  })) %>%
  dplyr::select(-GSVA) %>%
  tidyr::unnest() %>%
  dplyr::mutate(surv_res = purrr::map2(data,Features,fn_survival_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS
GSVA.score.univarite.surv.PFS %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score","GSVA.score.univarite.surv.PFS.tsv"))

# 3.4. multi-variate survival analysis ------
# 3.4.1.PFS ----
GSVA.score.univarite.surv.PFS %>%
  dplyr::filter(coxp<=0.05 | kmp <= 0.05) %>%
  dplyr::mutate(Features = gsub(" ","_",Features)) %>%
  dplyr::select(cancer_types,Features) %>%
  tidyr::nest(-cancer_types,.key="sig_features") -> GSVA.score.univarite.surv.PFS.sig

GSVA.score.onlytumor %>%
  dplyr::mutate(clinical_gsva = purrr::map2(GSVA,cancer_types,.f=function(.x,.y){
    print(.y)
    colnames(.x) <- gsub(" ","_",colnames(.x))
    .x %>%
      dplyr::mutate(barcode = substr(barcode,1,12)) %>%
      dplyr::inner_join(clinical_all_data,by="barcode")  %>%
      dplyr::rename("status"="PFS","time"="PFS.time") %>%
      dplyr::inner_join(stage_class,by="Stage") %>%
      dplyr::mutate(Stage = ifelse(group!="Not_applicable",group,NA)) 
  })) %>%
  dplyr::select(-GSVA) %>%
  dplyr::inner_join(GSVA.score.univarite.surv.PFS.sig,by="cancer_types") %>%
  dplyr::mutate(surv_res.multi = purrr::map2(clinical_gsva,sig_features,fn_survival_test.multiCox)) %>%
  dplyr::select(-clinical_gsva,-sig_features) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.multi
GSVA.score.univarite.surv.PFS.multi %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score","GSVA.score.multi-varite.surv.PFS.tsv"))

# 3.4.2.PFS cox plot -------
fn_cox_plot(data = GSVA.score.univarite.surv.PFS.multi %>%
              dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
              dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1),
            hr="hr",hr_l="hr_l",hr_h = "hr_h",title="Multi-variable cox model of GSVA score of features (PFS)",facet = ".~cancer_types", filename="1.PFS.multi-variable.cox",dir=file.path(res_path,"3.survival_with_GSVA_score"),w=20,h=15)

# save.image --------------------------------------------------------------

save.image(file.path(res_path,"TCGA.GSVA_score.ICP_features.rda"))
