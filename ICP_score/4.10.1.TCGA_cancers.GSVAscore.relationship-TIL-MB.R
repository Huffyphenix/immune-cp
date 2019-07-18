########################## GSVA score of features are: expression site of ICPs, gene family
########################## get the correlation between GSVA score and TIL, mutation burden
########################## cancer type specific
library(magrittr)
library(tidyverse)
library(survival)

# data path ---------------------------------------------------------------
# server 1
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path(basic_path,"immune_checkpoint/data/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"TCGA_GSVAScore/cancer_specific")
tcga_path <- file.path(basic_path,"/data/TCGA")

# load image --------------------------------------------------------------
load(file.path(res_path,"TCGA.GSVA_score.ICP_features-cancer.specific.rda"))

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
GSVA.TIL.res %>%
  tidyr::unnest() %>%
  readr::write_tsv(file.path(res_path,"1.TIL_with_GSVA_score","1.spm-cor.DE.TIL-GSVAscore.allcancers.tsv"))
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
  theme(
        plot.margin = unit(c(0,0,0,1),"cm"))+
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=frame.coord.all) +
  geom_text(aes(x=x,y=y,label=percentage),data = percentage.label.all) 


ggsave(file.path(res_path,"1.TIL_with_GSVA_score","4.DE.Cor.TIL-GSVAscore.allcancers.pdf"),device = "pdf",height = 4, width=8)
ggsave(file.path(res_path,"1.TIL_with_GSVA_score","4.DE.Cor.TIL-GSVAscore.allcancers.png"),device = "png",height = 4, width=8)

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

GSVA.MB.res %>%
  tidyr::unnest() %>%
  readr::write_tsv(file.path(res_path,"2.MB_with_GSVA_score","1.spm-cor.DE.MB-GSVAscore.allcancers.tsv"))

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


# 2.3.4.all MB DE in cancer types  --------
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
  theme(
        plot.margin = unit(c(0,0,0,1),"cm"))+
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=frame.coord.MB.all) +
  geom_text(aes(x=x,y=y,label=percentage),data = percentage.label.MB.all) 


ggsave(file.path(res_path,"2.MB_with_GSVA_score","4.DE.Cor.MB-GSVAscore.allcancers.pdf"),device = "pdf",height = 4, width=8)
ggsave(file.path(res_path,"2.MB_with_GSVA_score","4.DE.Cor.MB-GSVAscore.allcancers.png"),device = "png",height = 4, width=8)

# 4. GSVA score with CTL --------------------------------------------------
# 4.1.load CTL data -----
CTL_path <- "/home/huff/project/data/TCGA/CTL_level_estimated"
CTL_data <- readr::read_rds(file.path(CTL_path,"CTL_estimated_from_CD8A_CD8B_GZMA_GZMB_PRF1.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::rename("barcode" = "sample","CTL"="CTL(average_exp)") %>%
  dplyr::select(cancer_types,barcode,CTL)

# 4.2.calculate correlation ----
GSVA.score.onlytumor %>% 
  # head(1) %>%
  dplyr::mutate(res = purrr::map(GSVA,.f=fn_correlation_and_DE,
                                 TIL = CTL_data %>% 
                                   dplyr::select(barcode,`CTL(average_exp)`) %>%
                                   dplyr::mutate(barcode = substr(barcode,1,15)) %>%
                                   dplyr::rename("TIL"="CTL"))) %>%
  dplyr::select(-GSVA) %>%
  tidyr::unnest() -> GSVA.CTL.res

GSVA.CTL.res %>%
  tidyr::unnest() %>%
  readr::write_tsv(file.path(res_path,"4.CTL_with_GSVA_score","1.spm-cor.DE.CTL-GSVAscore.allcancers.tsv"))

# 4.3.all MB DE in cancer types  --------
GSVA.CTL.res %>%
  tidyr::unnest() %>%
  dplyr::rename("Correlation"="estimate") -> GSVA.CTL.DE.plot.all

tibble::tibble(x=c(0.3,1,1,0.3,-0.3,-1,-1,-0.3),
               y=c(1,1,3,3,1,1,3,3),
               xend=c(1,1,0.3,0.3,-1,-1,-0.3,-0.3),
               yend=c(1,3,3,1,1,3,3,1)) %>%
  dplyr::mutate(sig = ifelse(x<0,"neg_sig","pos_sig")) -> frame.coord.CTL.all

GSVA.CTL.DE.plot.all %>%
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
  dplyr::filter(sig!="not") -> percentage.label.CTL.all

GSVA.CTL.DE.plot.all %>%
  ggplot(aes(x=Correlation,y=`log2FC(High/Low)`)) +
  geom_jitter(aes(color=Features)) +
  scale_color_manual(
    values = c(cancer.color$color[1:length(unique(GSVA.CTL.DE.plot.all$Features))])
  ) +
  labs(x="Spearman correlation between \nCTL and GSVA score in tumors",
       y=latex2exp::TeX("log_{2} (CTL fold change between \n high and low GSVA score groups) in tumors")) +
  my_theme +
  theme(
        plot.margin = unit(c(0,0,0,1),"cm"))+
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=frame.coord.CTL.all) +
  geom_text(aes(x=x,y=y,label=percentage),data = percentage.label.CTL.all) 


ggsave(file.path(res_path,"4.CTL_with_GSVA_score","4.DE.Cor.TIL-GSVAscore.allcancers.pdf"),device = "pdf",height = 4, width=8)
ggsave(file.path(res_path,"4.CTL_with_GSVA_score","4.DE.Cor.TIL-GSVAscore.allcancers.png"),device = "png",height = 4, width=8)

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
    coef = .z$coefficients[,1],
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
      coef = model.s$coefficients[,1],
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
fn_cox_plot_1 <- function(data,filename,title,facet, dir,w=4,h=4){
  data %>% 
    # dplyr::mutate(functionWithImmune=functionWithImmune) %>%
    # dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
    ggplot(aes(y = hr, x = Features, ymin=hr_l,ymax=hr_h)) +
    geom_pointrange(aes(color=cox_sig),size=0.5) +
    scale_color_manual(values=c("red","black")) +
    geom_hline(aes(yintercept = 1), linetype =2) +
    scale_size(name = "p-value") +
    scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6),
                       labels = c("1/64","1/32","1/16","1/8","1/4","1/2",1,2,3,4,5,6)) +
    # facet_grid(as.formula(facet),scales = "free", space = "free") +
    facet_wrap(as.formula(facet),scales = "free") +
    coord_flip() +
    # ggthemes::theme_gdocs() +
    my_theme +
    theme(
      legend.position = "none",
      axis.line.y = element_line(color="black"),
      axis.text = element_text(color = "black",size=8),
      axis.title = element_text(color = "black",size=10),
      text = element_text(color = "black")
    ) +
    labs(y = "Hazard Ratio (High vs. low value)", x = "Features",title = title) -> p;p
  ggsave(file.path(dir,paste(filename,"png",sep=".")),device = "png",width = w,height = h)
  ggsave(file.path(dir,paste(filename,"pdf",sep=".")),device = "pdf",width = w,height = h)
}

fn_cox_plot_2 <- function(data,filename,title,facet, dir,w=4,h=4){
  data %>% 
    # dplyr::mutate(functionWithImmune=functionWithImmune) %>%
    # dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
    ggplot(aes(y = hr, x = Features, ymin=hr_l,ymax=hr_h)) +
    geom_pointrange(aes(color=cox_sig),size=0.5) +
    scale_color_manual(values=c("red","black")) +
    geom_hline(aes(yintercept = 1), linetype =2) +
    scale_size(name = "p-value") +
    scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6),
                       labels = c("1/64","1/32","1/16","1/8","1/4","1/2",1,2,3,4,5,6)) +
    facet_grid(as.formula(facet),scales = "free", space = "free") +
    # facet_wrap(as.formula(facet),scales = "free") +
    coord_flip() +
    # ggthemes::theme_gdocs() +
    my_theme +
    theme(
      legend.position = "none",
      axis.line.y = element_line(color="black"),
      axis.text = element_text(color = "black",size=8),
      axis.title = element_text(color = "black",size=10),
      text = element_text(color = "black")
    ) +
    labs(y = "Hazard Ratio (High vs. low value)", x = "Features",title = title) -> p;p
  ggsave(file.path(dir,paste(filename,"png",sep=".")),device = "png",width = w,height = h)
  ggsave(file.path(dir,paste(filename,"pdf",sep=".")),device = "pdf",width = w,height = h)
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
      dplyr::left_join(CTL_data,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,15)) %>%
      dplyr::left_join(TIMER_immunity,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,12)) %>%
      dplyr::left_join(TCGA_mutation_burden %>%
                         dplyr::mutate(barcode=substr(Tumor_Sample_Barcode,1,12)) %>%
                         dplyr::select(barcode,mutation_durden),by="barcode") %>%
      dplyr::inner_join(clinical_all_data,by=c("barcode","cancer_types")) -> .tmp
    if(nrow(.tmp)!=0){
      .tmp %>%
        dplyr::rename("status"="PFS","time"="PFS.time") %>%
        dplyr::select(-OS,-Status,-cancer_types,-Stage) %>%
        tidyr::gather(-barcode,-status,-time,key="Features",value="value")  %>%
        dplyr::filter(!is.na(value)) %>%
        dplyr::group_by(Features) %>%
        dplyr::mutate(group = ifelse(value >= quantile(value,0.5),"2high","1low")) %>%
        dplyr::mutate(n=n()) %>%
        dplyr::ungroup() %>%
        tidyr::nest(-Features,.key = group) %>%
        dplyr::mutate(group=purrr::map(group,.f=function(.x){
          .x %>%
            dplyr::mutate(nn=.x$group %>% unique() %>% length())
        })) %>%
        tidyr::unnest() %>%
        dplyr::filter(n>10,nn>=2) %>%
        dplyr::select(-n,-nn) %>%
        # rbind(stage.surv.data) %>%
        tidyr::nest(-Features) 
    } else{
      tibble::tibble()
    }
  })) %>%
  dplyr::select(-GSVA) %>%
  tidyr::unnest() %>%
  dplyr::mutate(surv_res = purrr::map2(data,Features,fn_survival_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS
GSVA.score.univarite.surv.PFS %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","GSVA.score.univarite.surv.PFS.tsv"))

# 3.3.2.OS ----
GSVA.score.onlytumor %>%
  dplyr::mutate(clinical_gsva = purrr::map2(GSVA,cancer_types,.f=function(.x,.y){
    print(.y)
    .x %>%
      dplyr::left_join(CTL_data,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,15)) %>%
      dplyr::left_join(TIMER_immunity,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,12)) %>%
      dplyr::left_join(TCGA_mutation_burden %>%
                         dplyr::mutate(barcode=substr(Tumor_Sample_Barcode,1,12)) %>%
                         dplyr::select(barcode,mutation_durden),by="barcode") %>%
      dplyr::inner_join(clinical_all_data,by=c("barcode","cancer_types")) -> .tmp
    if(nrow(.tmp)!=0){
      .tmp %>%
        dplyr::rename("status"="Status","time"="OS") %>%
        dplyr::select(-PFS,-PFS.time,-cancer_types,-Stage) %>%
        tidyr::gather(-barcode,-status,-time,key="Features",value="value")  %>%
        dplyr::filter(!is.na(value)) %>%
        dplyr::group_by(Features) %>%
        dplyr::mutate(group = ifelse(value >= quantile(value,0.5),"2high","1low")) %>%
        dplyr::mutate(n=n()) %>%
        dplyr::ungroup() %>%
        tidyr::nest(-Features,.key = group) %>%
        dplyr::mutate(group=purrr::map(group,.f=function(.x){
          .x %>%
            dplyr::mutate(nn=.x$group %>% unique() %>% length())
        })) %>%
        tidyr::unnest() %>%
        dplyr::filter(n>10,nn>=2) %>%
        dplyr::select(-n,-nn) %>%
        # rbind(stage.surv.data) %>%
        tidyr::nest(-Features) 
    } else{
      tibble::tibble()
    }
  })) %>%
  dplyr::select(-GSVA) %>%
  tidyr::unnest() %>%
  dplyr::mutate(surv_res = purrr::map2(data,Features,fn_survival_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS
GSVA.score.univarite.surv.OS %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","2.GSVA.score.univarite.surv.OS.tsv"))

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
      dplyr::left_join(CTL_data,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,15)) %>%
      dplyr::left_join(TIMER_immunity,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,12)) %>%
      dplyr::left_join(TCGA_mutation_burden %>%
                         dplyr::mutate(barcode=substr(Tumor_Sample_Barcode,1,12)) %>%
                         dplyr::select(barcode,mutation_durden),by="barcode") %>%
      dplyr::inner_join(clinical_all_data,by="barcode")  %>%
      dplyr::rename("status"="PFS","time"="PFS.time") 
  })) %>%
  dplyr::select(-GSVA) %>%
  dplyr::inner_join(GSVA.score.univarite.surv.PFS.sig,by="cancer_types") %>%
  dplyr::mutate(surv_res.multi = purrr::map2(clinical_gsva,sig_features,fn_survival_test.multiCox)) %>%
  dplyr::select(-clinical_gsva,-sig_features) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.multi
GSVA.score.univarite.surv.PFS.multi %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","GSVA.score.multi-varite.surv.PFS.tsv"))

# 3.4.2.OS ----
GSVA.score.univarite.surv.OS %>%
  dplyr::filter(coxp<=0.05 | kmp <= 0.05) %>%
  dplyr::mutate(Features = gsub(" ","_",Features)) %>%
  dplyr::select(cancer_types,Features) %>%
  tidyr::nest(-cancer_types,.key="sig_features") -> GSVA.score.univarite.surv.OS.sig

GSVA.score.onlytumor %>%
  dplyr::mutate(clinical_gsva = purrr::map2(GSVA,cancer_types,.f=function(.x,.y){
    print(.y)
    colnames(.x) <- gsub(" ","_",colnames(.x))
    .x %>%
      dplyr::left_join(CTL_data,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,15)) %>%
      dplyr::left_join(TIMER_immunity,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,12)) %>%
      dplyr::left_join(TCGA_mutation_burden %>%
                         dplyr::mutate(barcode=substr(Tumor_Sample_Barcode,1,12)) %>%
                         dplyr::select(barcode,mutation_durden),by="barcode") %>%
      dplyr::inner_join(clinical_all_data,by="barcode")  %>%
      dplyr::rename("status"="Status","time"="OS") 
  })) %>%
  dplyr::select(-GSVA) %>%
  dplyr::inner_join(GSVA.score.univarite.surv.OS.sig,by="cancer_types") %>%
  dplyr::mutate(surv_res.multi = purrr::map2(clinical_gsva,sig_features,fn_survival_test.multiCox)) %>%
  dplyr::select(-clinical_gsva,-sig_features) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS.multi
GSVA.score.univarite.surv.OS.multi %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","2.GSVA.score.multi-varite.surv.OS.tsv"))

# 3.5.1.OS cox plot -------
# 3.5.1.1. cancer groups 
GSVA.score.univarite.surv.OS.multi %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  dplyr::filter(abs(hr)<10) %>%
  dplyr::mutate(hr_l=ifelse(hr_l<(-10),-10,hr_l),hr_h=ifelse(hr_h>10,10,hr_h)) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(Group = c(rep("Group1",7),rep("Group2",7),rep("Group3",7),rep("Group4",8))) %>%
  tidyr::unnest() %>%
  tidyr::nest(-Group) %>%
  dplyr::mutate(filename = paste("2.OS.multi-variable.cox",Group,sep="_")) %>%
  dplyr::mutate(res = purrr::map2(data,filename,fn_cox_plot_2,facet = "cancer_types~.",title="Multi-variable cox model of GSVA score of features (OS)",dir=file.path(res_path,"3.survival_with_GSVA_score.new"),w=10,h=12))

# 3.5.1.2. all cancer togather
GSVA.score.univarite.surv.OS.multi %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  dplyr::filter(abs(hr)<6) %>%
  dplyr::mutate(hr_l=ifelse(hr_l<(-6),-6,hr_l),hr_h=ifelse(hr_h>6,6,hr_h)) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(Group = "All") %>%
  tidyr::unnest() %>%
  tidyr::nest(-Group) %>%
  dplyr::mutate(filename = paste("2.OS.multi-variable.cox",Group,sep="_")) %>%
  dplyr::mutate(res = purrr::map2(data,filename,fn_cox_plot_1,facet = "cancer_types~.",title="Multi-variable cox model of GSVA score of features (OS)",dir=file.path(res_path,"3.survival_with_GSVA_score.new"),w=25,h=20))

# 3.5.2.PFS cox plot -------
# 3.5.2.1. cancer groups 
GSVA.score.univarite.surv.PFS.multi %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  dplyr::filter(abs(hr)<10) %>%
  dplyr::mutate(hr_l=ifelse(hr_l<(-10),-10,hr_l),hr_h=ifelse(hr_h>10,10,hr_h)) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(Group = c(rep("Group1",7),rep("Group2",7),rep("Group3",7),rep("Group4",8))) %>%
  tidyr::unnest() %>%
  tidyr::nest(-Group) %>%
  dplyr::mutate(filename = paste("1.PFS.multi-variable.cox",Group,sep="_")) %>%
  dplyr::mutate(res = purrr::map2(data,filename,fn_cox_plot_2,facet = "cancer_types~.",title="Multi-variable cox model of GSVA score of features (PFS)",dir=file.path(res_path,"3.survival_with_GSVA_score.new"),w=10,h=12))

# 3.5.2.2. all cancer togather
GSVA.score.univarite.surv.PFS.multi %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  dplyr::filter(abs(hr)<6) %>%
  dplyr::mutate(hr_l=ifelse(hr_l<(-6),-6,hr_l),hr_h=ifelse(hr_h>6,6,hr_h)) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(Group = "All") %>%
  tidyr::unnest() %>%
  tidyr::nest(-Group) %>%
  dplyr::mutate(filename = paste("1.PFS.multi-variable.cox",Group,sep="_")) %>%
  dplyr::mutate(res = purrr::map2(data,filename,fn_cox_plot_1,facet = "cancer_types~.",title="Multi-variable cox model of GSVA score of features (PFS)",dir=file.path(res_path,"3.survival_with_GSVA_score.new"),w=25,h=20))

# 5. score from cox model -------------------------------------------------
# 5.1.function--------
fn_cox_score <- function(GSVA,survival){
  GSVA %>%
    tidyr::gather(-barcode,key="Features",value="GSVAscore") %>%
    dplyr::filter(Features %in% survival$Features) %>%
    dplyr::inner_join(survival,by="Features") %>%
    dplyr::mutate(GSVA_muliple_coef = GSVAscore*coef) %>%
    dplyr::group_by(barcode) %>%
    dplyr::mutate(cox_score = sum(GSVA_muliple_coef)) %>%
    dplyr::ungroup() %>%
    dplyr::select(barcode,cox_score) %>%
    unique()
}
# 5.2. multi-cox model coef ------
# 5.2.1.PFS --------
GSVA.score.univarite.surv.PFS.multi %>%
  tidyr::nest(-cancer_types,.key="PFS") %>%
  dplyr::inner_join(GSVA.score.onlytumor,by="cancer_types") %>%
  dplyr::mutate(cox_score = purrr::map2(GSVA,PFS,fn_cox_score)) %>%
  dplyr::select(cancer_types,cox_score) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.multi.coxscore

GSVA.score.univarite.surv.PFS.multi.coxscore %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","PFS_cox_score(GSVA*coef_of_multi_cox).tsv"))

GSVA.score.univarite.surv.PFS.multi.coxscore %>%
  dplyr::mutate(barcode = substr(barcode,1,15)) %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$TIL,method = "spearman"))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.multi.coxscore.correlation.TIL

GSVA.score.univarite.surv.PFS.multi.coxscore %>%
  # dplyr::mutate(barcode = substr(barcode,1,15)) %>%
  dplyr::inner_join(CTL_data,by=c("barcode","cancer_types")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$CTL,method = "spearman"))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.multi.coxscore.correlation.CTL

GSVA.score.univarite.surv.PFS.multi.coxscore %>%
  dplyr::mutate(sample = substr(barcode,1,15)) %>%
  dplyr::inner_join(TCGA_mutation_burden,
                    by=c("sample")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$mutation_durden,method = "spearman"))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.multi.coxscore.correlation.mutation_durden

# 5.2.2.OS --------

GSVA.score.univarite.surv.OS.multi %>%
  tidyr::nest(-cancer_types,.key="OS") %>%
  dplyr::inner_join(GSVA.score.onlytumor,by="cancer_types") %>%
  dplyr::mutate(cox_score = purrr::map2(GSVA,OS,fn_cox_score)) %>%
  dplyr::select(cancer_types,cox_score) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS.multi.coxscore
GSVA.score.univarite.surv.OS.multi.coxscore %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","OS_cox_score(GSVA*coef_of_multi_cox).tsv"))

GSVA.score.univarite.surv.OS.multi.coxscore %>%
  dplyr::mutate(barcode = substr(barcode,1,15)) %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$TIL,method = "spearman"))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS.multi.coxscore.correlation.TIL

GSVA.score.univarite.surv.OS.multi.coxscore %>%
  # dplyr::mutate(barcode = substr(barcode,1,15)) %>%
  dplyr::inner_join(CTL_data,by=c("barcode","cancer_types")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$CTL,method = "spearman"))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS.multi.coxscore.correlation.CTL

GSVA.score.univarite.surv.OS.multi.coxscore %>%
  dplyr::mutate(sample = substr(barcode,1,15)) %>%
  dplyr::inner_join(TCGA_mutation_burden,
                    by=c("sample")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$mutation_durden,method = "spearman"))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS.multi.coxscore.correlation.mutation_durden

# 5.3 uni-cox coef to get cox_score --------
# 5.2.1.PFS --------
GSVA.score.univarite.surv.PFS %>%
  tidyr::nest(-cancer_types,.key="PFS") %>%
  dplyr::inner_join(GSVA.score.onlytumor,by="cancer_types") %>%
  dplyr::mutate(cox_score = purrr::map2(GSVA,PFS,fn_cox_score)) %>%
  dplyr::select(cancer_types,cox_score) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.coxscore

GSVA.score.univarite.surv.PFS.coxscore %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","3.PFS_cox_score(GSVA*coef_of_multi_cox)-all_features.tsv"))

GSVA.score.univarite.surv.PFS.coxscore %>%
  dplyr::mutate(barcode = substr(barcode,1,15)) %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$TIL,method = "spearman"))
  })) %>%
  dplyr::select(-data)  -> GSVA.score.univarite.surv.PFS.coxscore.correlation.TIL

GSVA.score.univarite.surv.PFS.coxscore %>%
  dplyr::inner_join(CTL_data,by=c("barcode","cancer_types")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(CTL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$CTL,method = "spearman"))
  })) %>%
  dplyr::select(-data)-> GSVA.score.univarite.surv.PFS.coxscore.correlation.CTL

GSVA.score.univarite.surv.PFS.coxscore %>%
  dplyr::mutate(sample = substr(barcode,1,15)) %>%
  dplyr::inner_join(TCGA_mutation_burden,
                    by=c("sample")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(MB_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$mutation_durden,method = "spearman"))
  })) %>%
  dplyr::select(-data) -> GSVA.score.univarite.surv.PFS.coxscore.correlation.MB

GSVA.score.univarite.surv.PFS.coxscore.correlation.TIL %>%
  dplyr::inner_join(GSVA.score.univarite.surv.PFS.coxscore.correlation.CTL,by="cancer_types") %>%
  dplyr::inner_join(GSVA.score.univarite.surv.PFS.coxscore.correlation.MB,by="cancer_types") -> GSVA.score.univarite.surv.PFS.coxscore.correlation.combine

GSVA.score.univarite.surv.PFS.coxscore.correlation.combine %>%
  readr::write_rds(file.path(res_path,"3.survival_with_GSVA_score.new","3.PFS_cox_score(GSVA*coef_of_multi_cox)-all_features-correlation.combine(TIL.CTL.MB).rds.gz"),compress = "gz")

# 5.2.2.OS --------
GSVA.score.univarite.surv.OS %>%
  tidyr::nest(-cancer_types,.key="PFS") %>%
  dplyr::inner_join(GSVA.score.onlytumor,by="cancer_types") %>%
  dplyr::mutate(cox_score = purrr::map2(GSVA,PFS,fn_cox_score)) %>%
  dplyr::select(cancer_types,cox_score) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS.coxscore

GSVA.score.univarite.surv.OS.coxscore %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","3.OS_cox_score(GSVA*coef_of_multi_cox)-all_features.tsv"))

GSVA.score.univarite.surv.OS.coxscore %>%
  dplyr::mutate(barcode = substr(barcode,1,15)) %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$TIL,method = "spearman"))
  })) %>%
  dplyr::select(-data)  -> GSVA.score.univarite.surv.OS.coxscore.correlation.TIL

GSVA.score.univarite.surv.OS.coxscore %>%
  dplyr::inner_join(CTL_data,by=c("barcode","cancer_types")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(CTL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$CTL,method = "spearman"))
  })) %>%
  dplyr::select(-data)-> GSVA.score.univarite.surv.OS.coxscore.correlation.CTL

GSVA.score.univarite.surv.OS.coxscore %>%
  dplyr::mutate(sample = substr(barcode,1,15)) %>%
  dplyr::inner_join(TCGA_mutation_burden,
                    by=c("sample")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(MB_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$mutation_durden,method = "spearman"))
  })) %>%
  dplyr::select(-data) -> GSVA.score.univarite.surv.OS.coxscore.correlation.MB

GSVA.score.univarite.surv.OS.coxscore.correlation.TIL %>%
  dplyr::inner_join(GSVA.score.univarite.surv.OS.coxscore.correlation.CTL,by="cancer_types") %>%
  dplyr::inner_join(GSVA.score.univarite.surv.OS.coxscore.correlation.MB,by="cancer_types") -> GSVA.score.univarite.surv.OS.coxscore.correlation.combine

GSVA.score.univarite.surv.OS.coxscore.correlation.combine %>%
  readr::write_rds(file.path(res_path,"3.survival_with_GSVA_score.new","3.OS_cox_score(GSVA*coef_of_multi_cox)-all_features-correlation.combine(TIL.CTL.MB).rds.gz"),compress = "gz")

# 6. filter Features for each cancers -------------------------------------
GSVA.TIL.res %>%
  tidyr::unnest() %>%
  dplyr::filter(abs(estimate)>=0.3 & p.value<=0.05) %>%
  dplyr::filter(abs(`log2FC(High/Low)`)>=1 & p.value1<=0.05) %>%
  dplyr::mutate(Analysis = "TIL_cor_DE") %>%
  dplyr::select(cancer_types,Features,Analysis)-> GSVA.TIL.res.sig
GSVA.MB.res %>%
  tidyr::unnest() %>%
  dplyr::filter(abs(estimate)>=0.3 & p.value<=0.05) %>%
  dplyr::filter(abs(`log2FC(High/Low)`)>=1 & p.value1<=0.05) %>%
  dplyr::mutate(Analysis = "MB_cor_DE") %>%
  dplyr::select(cancer_types,Features,Analysis)-> GSVA.MB.res.sig
GSVA.CTL.res %>%
  tidyr::unnest() %>%
  dplyr::filter(abs(estimate)>=0.3 & p.value<=0.05) %>%
  dplyr::filter(abs(`log2FC(High/Low)`)>=1 & p.value1<=0.05)%>%
  dplyr::mutate(Analysis = "CTL_cor_DE") %>%
  dplyr::select(cancer_types,Features,Analysis) -> GSVA.CTL.res.sig
GSVA.score.univarite.surv.OS.multi %>%
  dplyr::filter(coxp<=0.1) %>%
  dplyr::mutate(Analysis = "Multi-cox,OS") %>%
  dplyr::select(cancer_types,Features,Analysis) -> GSVA.score.univarite.surv.OS.multi.sig
GSVA.score.univarite.surv.PFS.multi %>%
  dplyr::filter(coxp<=0.1) %>%
  dplyr::mutate(Analysis = "Multi-cox,PFS") %>%
  dplyr::select(cancer_types,Features,Analysis) -> GSVA.score.univarite.surv.PFS.multi.sig

GSVA.TIL.res.sig %>%
  rbind(GSVA.MB.res.sig) %>%
  rbind(GSVA.CTL.res.sig) %>%
  rbind(GSVA.score.univarite.surv.OS.multi.sig) %>%
  rbind(GSVA.score.univarite.surv.PFS.multi.sig) %>%
  dplyr::group_by(Features) %>%
  dplyr::mutate(Sig_n = n()) %>%
  dplyr::select(Features,Sig_n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Sig_n) %>%
  dplyr::rename("Counts_of_significant_analysis_in_all_cancers"="Sig_n") %>%
  dplyr::filter(!Features %in% c(colnames(TIMER_immunity),"CTL","Age")) %>%
  dplyr::mutate(color = ifelse(Counts_of_significant_analysis_in_all_cancers>60,"red","black")) -> GSVA_score_res.counts

GSVA_score_res.counts %>%
  readr::write_tsv(file.path(res_path,"overall_feature_filter.tsv"))

GSVA_score_res.counts %>%
  # dplyr::filter(cancer_types %in% c("SKCM","STAD","LUSC","LUAD","KIRC","KICH","KIRP")) %>%
  ggplot(aes(x=Features,y=Counts_of_significant_analysis_in_all_cancers)) +
  geom_bar(aes(color=color),stat="identity", fill="orange",alpha = 0.5,
           position=position_dodge()) +
  scale_x_discrete(limits=GSVA_score_res.counts$Features) +
  scale_color_manual(values=c("black","red")) +
  coord_flip() +
  my_theme +
  theme(
    axis.text.y = element_text(color = GSVA_score_res.counts$color),
    legend.position = "none"
  ) +
  labs(y="Counts of significant analysis among cancers") 
ggsave(file.path(res_path,"overall_feature_filter.pdf"),device = "pdf",height = 8,width =6)
ggsave(file.path(res_path,"overall_feature_filter.png"),device = "png",height = 8,width = 6)

# save.image --------------------------------------------------------------

save.image(file.path(res_path,"TCGA.GSVA_score.ICP_features-cancer.specific.rda"))
