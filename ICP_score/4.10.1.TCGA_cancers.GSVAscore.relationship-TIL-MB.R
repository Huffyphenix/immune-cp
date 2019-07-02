########################## GSVA score of features are: expression site of ICPs, gene family
########################## get the correlation between GSVA score and TIL, mutation burden
library(magrittr)
library(tidyverse)

# data path ---------------------------------------------------------------
# server 1
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path(basic_path,"immune_checkpoint/data/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"TCGA_GSVAScore")
tcga_path <- file.path(basic_path,"/data/TCGA")


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

GSVA.TIL.cor.plot %>%
  ggplot(aes(x=Correlation,y=`-log10P`)) +
  geom_jitter(aes(color=cancer_types)) +
  scale_color_manual(
    values = cancer.color$color,
    name = "Cancer types"
  ) +
  labs(x="Spearman correlation between TIL and GSVA score",
       y="-log10P") +
  my_theme
ggsave(file.path(res_path,"1.TIL_with_GSVA_score","1.spm.cor.TIL-GSVAscore.allcancers.pdf"),device = "pdf",height = 8, width=6)
ggsave(file.path(res_path,"1.TIL_with_GSVA_score","1.spm.cor.TIL-GSVAscore.allcancers.png"),device = "png",height = 8, width=6)

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
ggsave(file.path(res_path,"1.TIL_with_GSVA_score","2.DE.TIL-GSVAscore.allcancers.pdf"),device = "pdf",height = 8, width=6)
ggsave(file.path(res_path,"1.TIL_with_GSVA_score","2.DE.TIL-GSVAscore.allcancers.png"),device = "png",height = 8, width=6)

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
  dplyr::ungroup()-> label.data.MB.cor

GSVA.MB.cor.plot %>%
  ggplot(aes(x=Correlation,y=`-log10P`)) +
  geom_jitter(aes(color=cancer_types)) +
  scale_color_manual(
    values = cancer.color.MB$color,
    name = "Cancer types"
  ) +
  geom_text_repel(aes(x=Correlation,y=`-log10P`,label=label,color=cancer_types),
                  data=label.data.MB.cor, size = 2) +
  labs(x="Spearman correlation between \nmutation burden and GSVA score",
       y="-log10(P value)") +
  my_theme
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


# save.image --------------------------------------------------------------

save.image(file.path(res_path,"TCGA.GSVA_score.ICP_features.rda"))
