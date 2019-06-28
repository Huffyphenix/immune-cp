###########################################################################
# survival situation of each cancers [with diff of TIL(T-N)]------------------------------

library(magrittr)
library(ggplot2)

# data path config ---------------------------------------------------
basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
tcga_path <- file.path(basic_path,"/data/TCGA")
# gene_list_path <- file.path(immune_path,"checkpoint/20171021_checkpoint")
# expr_path <- file.path(basic_path,"immune_checkpoint/result_20171025/expr_rds")
result_path <- file.path(basic_path,"immune_checkpoint/result_20171025/e_5_immune_infiltration")

load(
  file.path(out_path,"pan-can_survival-TIL(T-N).rds")
)
# load data ---------------------------------------------------------------
## survival data
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
  dplyr::rename("cancer_types"= "type")-> clinical_all_data

## TIL data
# 2. data from TIMER
immunity_path_2 <- file.path(basic_path,"immune_checkpoint/data/immunity")
TIMER_immunity <- readr::read_tsv(file.path(immunity_path_2,"immuneEstimation.txt")) %>%
  dplyr::mutate(TIL = B_cell+CD4_Tcell+CD8_Tcell+Neutrophil+Macrophage+Dendritic)

## TCGA sample info
TCGA_cancer_info <- readr::read_tsv(file.path(tcga_path, "TCGA_RNAseq_sample_info.tsv")) 
TCGA_cancer_info %>%
  dplyr::mutate(barcode_short = substr(barcode,1,12)) %>%
  dplyr::filter(substr(barcode,14,14)=="1") -> TCGA.RNASeq.normal

TCGA_cancer_info %>%
  dplyr::mutate(barcode_short = substr(barcode,1,12)) %>%
  dplyr::filter(barcode_short %in% TCGA.RNASeq.normal$barcode_short) %>%
  dplyr::mutate(class = ifelse(substr(barcode,14,15)=="01","Tumor","Normal")) %>%
  dplyr::mutate(class = ifelse(substr(barcode,14,15) %in% c("01","11"), class, "others")) %>%
  dplyr::filter(class != "others") -> TCGA.RNAseq.TN.paired


# TIL diff between tumor and normal samples -------------------------------
TCGA.RNAseq.TN.paired %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  dplyr::select(cancer_types,barcode_short,TIL,class) %>%
  tidyr::spread(key="class",value="TIL") %>%
  dplyr::mutate(TIL.diff = Tumor-Normal) -> TCGA.RNAseq.TN.paired.TILdiff


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
cancer_colors <- readr::read_tsv(file.path(tcga_path,"02_pcc.tsv"))

# 1.distribution TIL diff  ------------------------------------------------
TCGA.RNAseq.TN.paired.TILdiff %>%
  dplyr::filter(!is.na(TIL.diff)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(TIL.diff,0.5)) %>%
  dplyr::mutate(n=dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cancer_types.label = paste(cancer_types," (n=",n,")",sep="")) %>%
  dplyr::select(cancer_types,cancer_types.label,mid) %>%
  dplyr::arrange(mid) %>%
  unique()-> cancer_rank

cancer_colors %>%
  dplyr::inner_join(TCGA.RNAseq.TN.paired.TILdiff,by="cancer_types") %>%
  dplyr::select(cancer_types,color) %>%
  unique() -> color_paired

TCGA.RNAseq.TN.paired.TILdiff %>%
  dplyr::filter(!is.na(TIL.diff)) %>%
  dplyr::inner_join(cancer_rank,by="cancer_types") %>%
  ggplot(aes(x=cancer_types.label,y=TIL.diff)) +
  # geom_jitter(size=0.5,width = 0.1) +
  scale_x_discrete(limits= cancer_rank$cancer_types.label) +
  geom_boxplot(aes(fill = cancer_types),alpha=0.5) +
  coord_flip()  +
  labs(y= "TIL difference (Tumor - Normal)", x="Cancer types") +
  scale_fill_manual(values= color_paired$color) +
  my_theme +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    legend.position = "none",
    axis.text.x = element_text(colour = "black")
  )
filename <- paste("TIMER","boxplot.TIL(T-N)",sep=".")
ggsave(file.path(result_path,paste(filename,"png",sep=".")),device = "png",width = 4,height = 6)
ggsave(file.path(result_path,paste(filename,"pdf",sep=".")),device = "pdf",width = 4,height = 6)

# 2.survival analysis using TIL diff  ------------------------------------------------
TCGA.RNAseq.TN.paired.TILdiff %>%
  dplyr::filter(!is.na(TIL.diff)) %>%
  dplyr::rename("barcode"="barcode_short") %>%
  dplyr::inner_join(clinical_all_data,by="barcode") -> TCGA_combine_TIL_survival

source("/project/huff/huff/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

## function
fn_durvival_plot <- function(data,title="Pan-cancers progress-free survival",color,filename="Pan-can.PFS.survival",out_path,legend.pos,h=3,w=4){
  
  fit <- survfit(Surv(time, status) ~ group, data = data, na.action = na.exclude)
  tmp <- 0
  for (i in 1:length(fit$strata)) {
    tmp <- fit$strata[i] + tmp
    .group <- strsplit(names(fit$strata[i]),split = "=")[[1]][2]
    tibble::tibble(mid_time=fit$time,final_surv=fit$surv) %>%
      head(tmp) %>%
      tail(1) %>%
      dplyr::mutate(group = .group)-> res.tmp
    
    if(i==1){
      tail_position <- res.tmp
    } else{
      tail_position <- rbind(tail_position,res.tmp)
    }
  }
  diff <- survdiff(Surv(time, status) ~ group, data = data, na.action = na.exclude)
  kmp <- 1 - pchisq(diff$chisq, df = length(levels(as.factor(data$group))) - 1) %>% signif(2)
  color %>%
    dplyr::inner_join(data,by="group") %>%
    dplyr::inner_join(tail_position,by="group") %>%
    dplyr::mutate(label = paste(group,", n=",n,sep="")) %>%
    dplyr::select(group,label,color,n,mid_time,final_surv) %>%
    dplyr::arrange(group) %>%
    unique() -> text_data
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        data = data,
                        # surv.median.line = "hv",
                        title = paste(title,", p=",kmp,sep=""), # change it when doing diff data
                        xlab = "Time (days)",
                        ylab = 'Probability of survival',
                        # legend.title = "Methyla group:",
                        legend= c(legend.pos),
                        # ggtheme = theme_survminer(),
                        ggtheme = theme(
                          panel.border = element_blank(), 
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), 
                          axis.line = element_line(colour = "black",size = 0.5),
                          panel.background = element_rect(fill = "white"),
                          legend.key = element_blank(),
                          legend.background = element_blank(),
                          legend.text = element_text(size = 8),
                          axis.text = element_text(size = 12,color = "black"),
                          legend.title = element_blank(),
                          axis.title = element_text(size = 12,color = "black"),
                          text = element_text(color = "black")
                        )
  )[[1]] +
    # geom_text(aes(x=time,y=surv,label=label),data=text_data) +
    geom_text_repel(aes(x=mid_time,y=final_surv,label=label,color=group),data=text_data) +
    scale_color_manual(
      values = c(text_data$color,text_data$color),
      labels = c(text_data$group,text_data$group)
    ) -> p;p
  ggsave(filename = paste(filename,kmp,"pdf",sep = "."), plot = p, path = out_path,device = "pdf",height = h,width = w)
  ggsave(filename = paste(filename,kmp,"png",sep = "."), plot = p, path = out_path,device = "png",height = h,width = w)
}
## 2.1.all sample togather, cancer_types as survival groups ----

dir.create(file.path(result_path,"TILdiff-TN_survival"))
sur_res_path <- file.path(result_path,"TILdiff-TN_survival")

### 2.1.1.OS -----
TCGA_combine_TIL_survival %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = cancer_types.x) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(n = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n>10)  -> TCGA_combine_TIL_survival.OS

cancer_colors %>%
  dplyr::rename("group" = "cancer_types") -> color_list

sur_name <- paste("OS_between_cancers")
TCGA_combine_TIL_survival.OS %>%
  fn_durvival_plot(title="Overall survival",color=color_list,filename=sur_name,out_path=sur_res_path,legend.pos="bottom",h=6,w=6)

### 2.1.1.PFS -----
TCGA_combine_TIL_survival %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = cancer_types.x) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(n = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n>10)  -> TCGA_combine_TIL_survival.PFS

cancer_colors %>%
  dplyr::rename("group" = "cancer_types") -> color_list

sur_name <- paste("PFS_between_cancers")
TCGA_combine_TIL_survival.PFS %>%
  fn_durvival_plot(title="Progression-free survival",color=color_list,filename=sur_name,out_path=sur_res_path,legend.pos="bottom",h=6,w=6)

## 2.2.all sample togather, TIL high and low as survival groups ----
## 2.2.1.PFS----
TCGA_combine_TIL_survival %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(TIL.diff>=quantile(TIL.diff,0.5),"High","Low")) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(n = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n>10)  -> TCGA_combine_TIL_survival.PFS.TILgroups

color_list <- tibble::tibble(group=c("High","Low"),
                             color=c("red","blue"))

sur_name <- paste("PFS_between_TILhigh-low")
TCGA_combine_TIL_survival.PFS.TILgroups %>%
  fn_durvival_plot(title="Progression-free survival",color=color_list,filename=sur_name,out_path=sur_res_path,legend.pos="none",h=3,w=4)

## 2.2.1.OS----
TCGA_combine_TIL_survival %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(TIL.diff>=quantile(TIL.diff,0.5),"High","Low")) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(n = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n>10)  -> TCGA_combine_TIL_survival.OS.TILgroups

color_list <- tibble::tibble(group=c("High","Low"),
                             color=c("red","blue"))

sur_name <- paste("OS_between_TILhigh-low")
TCGA_combine_TIL_survival.OS.TILgroups %>%
  fn_durvival_plot(title="Overall survival",color=color_list,filename=sur_name,out_path=sur_res_path,legend.pos="none",h=3,w=4)

## 2.3.cancer specific survival analysis, TIL high and low as survival groups ----
## 2.3.1.PFS----
TCGA_combine_TIL_survival %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::group_by(cancer_types.x) %>%
  dplyr::mutate(n = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n>10) %>%
  tidyr::nest(-cancer_types.x) %>%
  dplyr::mutate(data = purrr::map(data,.f = function(.x){
    .x %>%
      dplyr::mutate(group = ifelse(TIL.diff>=quantile(TIL.diff,0.5),"High","Low")) %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(n=dplyr::n()) %>%
      dplyr::ungroup()
  }))  -> TCGA_combine_TIL_survival.PFS.TILgroups.cancerSpecific

color_list <- tibble::tibble(group=c("High","Low"),
                             color=c("red","blue"))

sur_name <- paste("PFS_between_TILhigh-low")
TCGA_combine_TIL_survival.PFS.TILgroups.cancerSpecific %>%
  dplyr::mutate(filename = paste(cancer_types.x,sur_name,sep=".")) %>%
  dplyr::mutate(title = paste(cancer_types.x,"\n","Progression-free survival",sep="")) %>%
  dplyr::select(-cancer_types.x ) %>%
  purrr::pwalk(.f=fn_durvival_plot,color=color_list,out_path=sur_res_path,legend.pos="none",h=3,w=4)

## 2.3.2.OS----
TCGA_combine_TIL_survival %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::group_by(cancer_types.x) %>%
  dplyr::mutate(n = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n>10) %>%
  dplyr::select(-n) %>%
  tidyr::nest(-cancer_types.x) %>%
  dplyr::mutate(data = purrr::map(data,.f = function(.x){
    .x %>%
      dplyr::mutate(group = ifelse(TIL.diff>=quantile(TIL.diff,0.5),"High","Low")) %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(n=dplyr::n()) %>%
      dplyr::ungroup()
  }))  -> TCGA_combine_TIL_survival.OS.TILgroups.cancerSpecific

color_list <- tibble::tibble(group=c("High","Low"),
                             color=c("red","blue"))

sur_name <- paste("OS_between_TILhigh-low")
TCGA_combine_TIL_survival.OS.TILgroups.cancerSpecific %>%
  dplyr::mutate(filename = paste(cancer_types.x,sur_name,sep=".")) %>%
  dplyr::mutate(title = paste(cancer_types.x,"\n","Overall survival",sep="")) %>%
  dplyr::select(-cancer_types.x ) %>%
  purrr::pwalk(.f=fn_durvival_plot,color=color_list,out_path=sur_res_path,legend.pos="none",h=3,w=4)
