library(magrittr)
library(ggplot2)
library(ggpubr)

# data path 
basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")

tumor_class_by_T_N.only_paired.count <- readr::read_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.only_paired.by_geneCounts","tumor_class_by_T_N.only_paired"))
tumor_class_by_T_N.only_paired.score <- readr::read_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.only_paired.by_FCProduct","tumor_class_by_T_N.only_paired"))

# tumor_class_by_T_N.by_mean <- readr::read_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.by_mean"))

# tumor_class_by_T_N.by_peak <- readr::read_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.by_peak"))

source("/project/huff/huff/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

# get immunity score (gene percent) ---------------------------------------------------------------

tumor_class_by_T_N.only_paired.count %>%
  dplyr::select(-`<NA>`) %>%
  dplyr::mutate(hot_per = Immunity_hot/(Immunity_hot+Immunity_cold)) -> tumor_class_by_T_N.only_paired.genePercent

# tumor_class_by_T_N.by_mean %>%
#   dplyr::mutate(hot_per = Immunity_hot/(Immunity_hot+Immunity_cold)) -> tumor_class_by_T_N.by_mean.immunityScore
# 
# tumor_class_by_T_N.by_peak %>%
#   dplyr::mutate(hot_per = Immunity_hot/(Immunity_hot+Immunity_cold)) -> tumor_class_by_T_N.by_peak.immunityScore


# survival analysis -------------------------------------------------------
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
  dplyr::select(bcr_patient_barcode,PFS,PFS.time) %>%
  tidyr::drop_na() %>%
  dplyr::rename("barcode" = "bcr_patient_barcode") %>%
  unique()


# 1.only paired gene percent score analysis ----------------------------------------
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern/tumor_class_by_T_N.only_paired.by_geneCounts/")

# distribution ------------------------------------------------------------
tumor_class_by_T_N.only_paired.genePercent %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(hot_per,0.5)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,mid) %>%
  dplyr::arrange(mid) %>%
  unique()-> cancer_rank

tumor_class_by_T_N.only_paired.genePercent %>%
  dplyr::filter(!is.na(hot_per)) %>%
  ggplot(aes(x=cancer_types,y=hot_per)) +
  # geom_jitter(size=0.5,width = 0.1) +
  scale_x_discrete(limits= cancer_rank$cancer_types) +
  geom_violin(aes(fill = cancer_types),alpha=0.5) +
  rotate() +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12)
  )
ggsave(file.path(res_path,"score_distribution.png"),device = "png",width = 4,height = 6)
ggsave(file.path(res_path,"score_distribution.pdf"),device = "pdf",width = 4,height = 6)

cancer_color <- readr::read_tsv(file.path(basic_path,"data/TCGA","02_pcc.tsv"))
cancer_color %>%
  dplyr::filter(cancer_types %in% unique(tumor_class_by_T_N.only_paired.immunityScore$cancer_types)) -> cancer21_color

tumor_class_by_T_N.only_paired.genePercent %>%
  dplyr::filter(!is.na(hot_per)) %>%
  ggdensity(x="hot_per",fill = "cancer_types",alpha = 0.5,palette = "jco")+
  facet_wrap(~cancer_types) +
  # geom_vline(data=tumor_class_by_T_N.by_max.immunityScore.peak,aes(xintercept = peak.x),linetype = "dotted") +
  # geom_text(data=immune_lanscape_immunity.peak,aes(x=peak.x+0.2,y=peak.y+2,label=paste("(",signif(peak.x,2),",",signif(peak.y,2),")",sep=""))) +
  scale_fill_manual(
    values = cancer21_color$color,
    limits = cancer21_color$cancer_types
  ) +
  ggtitle("Density of ratio of ICP indicate immune hot in each cancers") +
  ylab("Density") +
  xlab("Immune Activity Score") +
  theme(
    legend.position = "none",
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(fill="white",color="black"),
    strip.text = element_text(color="black",size=12)
  )
ggsave(file.path(res_path,"score_distribution_in_cancers.png"),device = "png",height = 6,width = 8)
ggsave(file.path(res_path,"score_distribution_in_cancers.pdf"),device = "pdf",height = 6,width = 8)

##### all sample together survival -----------------------------------
# 1.only paired sample score survival -----

tumor_class_by_T_N.only_paired.genePercent %>%
  dplyr::rename("barcode" = "Participant") %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.only_paired.immunityScore.clinical

color_list <- tibble::tibble(group=c("High","Low"),
                             color=c("red","blue"))
sur_name <- paste("Score_0.5_Tumor_5year_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE",c("#CDAD00")),
                            group=c("High","Low","Mid"))

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

##### cancer specific survival -----------------------------------
# only paired sample score cacner specific survival -----

tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,group,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="group",value="n") %>%
  dplyr::filter(High >5) %>%
  .$cancer_types -> cancer_to_do_survival.only_paired

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE"),
                             group=c("High","Low"))
for(cancer in cancer_to_do_survival.only_paired){
  if(!dir.exists(file.path(res_path,"survival_cancers"))){
    dir.create(file.path(res_path,"survival_cancers"))
  }
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_PFS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_OS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
}

res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern/tumor_class_by_T_N.only_paired.by_geneCounts/")

# distribution ------------------------------------------------------------
tumor_class_by_T_N.only_paired.genePercent %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(hot_per,0.5)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,mid) %>%
  dplyr::arrange(mid) %>%
  unique()-> cancer_rank

tumor_class_by_T_N.only_paired.genePercent %>%
  dplyr::filter(!is.na(hot_per)) %>%
  ggplot(aes(x=cancer_types,y=hot_per)) +
  # geom_jitter(size=0.5,width = 0.1) +
  scale_x_discrete(limits= cancer_rank$cancer_types) +
  geom_violin(aes(fill = cancer_types),alpha=0.5) +
  rotate() +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12)
  )
ggsave(file.path(res_path,"score_distribution.png"),device = "png",width = 4,height = 6)
ggsave(file.path(res_path,"score_distribution.pdf"),device = "pdf",width = 4,height = 6)

cancer_color <- readr::read_tsv(file.path(basic_path,"data/TCGA","02_pcc.tsv"))
cancer_color %>%
  dplyr::filter(cancer_types %in% unique(tumor_class_by_T_N.only_paired.immunityScore$cancer_types)) -> cancer21_color

tumor_class_by_T_N.only_paired.genePercent %>%
  dplyr::filter(!is.na(hot_per)) %>%
  ggdensity(x="hot_per",fill = "cancer_types",alpha = 0.5,palette = "jco")+
  facet_wrap(~cancer_types) +
  # geom_vline(data=tumor_class_by_T_N.by_max.immunityScore.peak,aes(xintercept = peak.x),linetype = "dotted") +
  # geom_text(data=immune_lanscape_immunity.peak,aes(x=peak.x+0.2,y=peak.y+2,label=paste("(",signif(peak.x,2),",",signif(peak.y,2),")",sep=""))) +
  scale_fill_manual(
    values = cancer21_color$color,
    limits = cancer21_color$cancer_types
  ) +
  ggtitle("Density of ratio of ICP indicate immune hot in each cancers") +
  ylab("Density") +
  xlab("Immune Activity Score") +
  theme(
    legend.position = "none",
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(fill="white",color="black"),
    strip.text = element_text(color="black",size=12)
  )
ggsave(file.path(res_path,"score_distribution_in_cancers.png"),device = "png",height = 6,width = 8)
ggsave(file.path(res_path,"score_distribution_in_cancers.pdf"),device = "pdf",height = 6,width = 8)

##### all sample together survival -----------------------------------
# 1.only paired sample score survival -----

tumor_class_by_T_N.only_paired.genePercent %>%
  dplyr::rename("barcode" = "Participant") %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.only_paired.immunityScore.clinical

color_list <- tibble::tibble(group=c("High","Low"),
                             color=c("red","blue"))
sur_name <- paste("Score_0.5_Tumor_5year_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE",c("#CDAD00")),
                             group=c("High","Low","Mid"))

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

##### cancer specific survival -----------------------------------
# only paired sample score cacner specific survival -----

tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,group,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="group",value="n") %>%
  dplyr::filter(High >5) %>%
  .$cancer_types -> cancer_to_do_survival.only_paired

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE"),
                             group=c("High","Low"))
for(cancer in cancer_to_do_survival.only_paired){
  if(!dir.exists(file.path(res_path,"survival_cancers"))){
    dir.create(file.path(res_path,"survival_cancers"))
  }
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_PFS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_OS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
}

########## >>>>>>>>>>>>>>> ##############
# 1.1.only paired gene percent score analysis ----------------------------------------
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern/tumor_class_by_T_N.only_paired.by_FCProduct")

# distribution ------------------------------------------------------------
tumor_class_by_T_N.only_paired.score %>%
  dplyr::filter(score_mean!=0) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(score_mean,0.5)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,mid) %>%
  dplyr::arrange(mid) %>%
  unique()-> cancer_rank
library(ggbeeswarm)
tumor_class_by_T_N.only_paired.score %>%
  dplyr::filter(score_mean!=0) %>%
  ggplot(aes(x=cancer_types,y=score_mean)) +
  # geom_jitter(size=0.5,width = 0.1) +
  geom_jitter(aes(color=cancer_types),size=0.5,width=0.1) +
  scale_x_discrete(limits= cancer_rank$cancer_types) +
  geom_violin(aes(fill = cancer_types),alpha=0.5) +
  rotate() +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12)
  )
ggsave(file.path(res_path,"score_distribution.png"),device = "png",width = 4,height = 6)
ggsave(file.path(res_path,"score_distribution.pdf"),device = "pdf",width = 4,height = 6)

cancer_color <- readr::read_tsv(file.path(basic_path,"data/TCGA","02_pcc.tsv"))
cancer_color %>%
  dplyr::filter(cancer_types %in% unique(tumor_class_by_T_N.only_paired.immunityScore$cancer_types)) -> cancer21_color

tumor_class_by_T_N.only_paired.score %>%
  dplyr::filter(score_mean!=0) %>%
  ggdensity(x="score_mean",fill = "cancer_types",alpha = 0.5,palette = "jco")+
  facet_wrap(~cancer_types) +
  # geom_vline(data=tumor_class_by_T_N.by_max.immunityScore.peak,aes(xintercept = peak.x),linetype = "dotted") +
  # geom_text(data=immune_lanscape_immunity.peak,aes(x=peak.x+0.2,y=peak.y+2,label=paste("(",signif(peak.x,2),",",signif(peak.y,2),")",sep=""))) +
  scale_fill_manual(
    values = cancer21_color$color,
    limits = cancer21_color$cancer_types
  ) +
  ggtitle("Density of ICP score in each cancers") +
  ylab("Density") +
  xlab("Immune Activity Score") +
  theme(
    legend.position = "none",
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(fill="white",color="black"),
    strip.text = element_text(color="black",size=12)
  )
ggsave(file.path(res_path,"score_distribution_in_cancers.png"),device = "png",height = 6,width = 8)
ggsave(file.path(res_path,"score_distribution_in_cancers.pdf"),device = "pdf",height = 6,width = 8)

##### all sample together survival -----------------------------------
# 1.only paired sample score survival -----

tumor_class_by_T_N.only_paired.score %>%
  dplyr::rename("barcode" = "Participant") %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.only_paired.Score.clinical

color_list <- tibble::tibble(group=c("High","Low"),
                             color=c("red","blue"))
sur_name <- paste("Score_0.5_Tumor_5year_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_mean!=0) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"High","Low")) %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_mean!=0) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE",c("#CDAD00")),
                             group=c("High","Low","Mid"))

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_mean!=0) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(score_mean<quantile(score_mean,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_mean!=0) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(score_mean<quantile(score_mean,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_mean!=0) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

##### cancer specific survival -----------------------------------
# only paired sample score cacner specific survival -----

tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_mean!=0) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"High","Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,group,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="group",value="n") %>%
  dplyr::filter(High >5) %>%
  .$cancer_types -> cancer_to_do_survival.only_paired

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE"),
                             group=c("High","Low"))
for(cancer in cancer_to_do_survival.only_paired){
  if(!dir.exists(file.path(res_path,"survival_cancers"))){
    dir.create(file.path(res_path,"survival_cancers"))
  }
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_PFS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.Score.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_OS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.Score.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
}

res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern/tumor_class_by_T_N.only_paired.by_geneCounts/")

# distribution ------------------------------------------------------------
tumor_class_by_T_N.only_paired.genePercent %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(hot_per,0.5)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,mid) %>%
  dplyr::arrange(mid) %>%
  unique()-> cancer_rank

tumor_class_by_T_N.only_paired.genePercent %>%
  dplyr::filter(!is.na(hot_per)) %>%
  ggplot(aes(x=cancer_types,y=hot_per)) +
  # geom_jitter(size=0.5,width = 0.1) +
  scale_x_discrete(limits= cancer_rank$cancer_types) +
  geom_violin(aes(fill = cancer_types),alpha=0.5) +
  rotate() +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12)
  )
ggsave(file.path(res_path,"score_distribution.png"),device = "png",width = 4,height = 6)
ggsave(file.path(res_path,"score_distribution.pdf"),device = "pdf",width = 4,height = 6)

cancer_color <- readr::read_tsv(file.path(basic_path,"data/TCGA","02_pcc.tsv"))
cancer_color %>%
  dplyr::filter(cancer_types %in% unique(tumor_class_by_T_N.only_paired.immunityScore$cancer_types)) -> cancer21_color

tumor_class_by_T_N.only_paired.genePercent %>%
  dplyr::filter(!is.na(hot_per)) %>%
  ggdensity(x="hot_per",fill = "cancer_types",alpha = 0.5,palette = "jco")+
  facet_wrap(~cancer_types) +
  # geom_vline(data=tumor_class_by_T_N.by_max.immunityScore.peak,aes(xintercept = peak.x),linetype = "dotted") +
  # geom_text(data=immune_lanscape_immunity.peak,aes(x=peak.x+0.2,y=peak.y+2,label=paste("(",signif(peak.x,2),",",signif(peak.y,2),")",sep=""))) +
  scale_fill_manual(
    values = cancer21_color$color,
    limits = cancer21_color$cancer_types
  ) +
  ggtitle("Density of ratio of ICP indicate immune hot in each cancers") +
  ylab("Density") +
  xlab("Immune Activity Score") +
  theme(
    legend.position = "none",
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(fill="white",color="black"),
    strip.text = element_text(color="black",size=12)
  )
ggsave(file.path(res_path,"score_distribution_in_cancers.png"),device = "png",height = 6,width = 8)
ggsave(file.path(res_path,"score_distribution_in_cancers.pdf"),device = "pdf",height = 6,width = 8)

##### all sample together survival -----------------------------------
# 1.only paired sample score survival -----

tumor_class_by_T_N.only_paired.genePercent %>%
  dplyr::rename("barcode" = "Participant") %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.only_paired.immunityScore.clinical

color_list <- tibble::tibble(group=c("High","Low"),
                             color=c("red","blue"))
sur_name <- paste("Score_0.5_Tumor_5year_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE",c("#CDAD00")),
                             group=c("High","Low","Mid"))

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

##### cancer specific survival -----------------------------------
# only paired sample score cacner specific survival -----

tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,group,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="group",value="n") %>%
  dplyr::filter(High >5) %>%
  .$cancer_types -> cancer_to_do_survival.only_paired

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE"),
                             group=c("High","Low"))
for(cancer in cancer_to_do_survival.only_paired){
  if(!dir.exists(file.path(res_path,"survival_cancers"))){
    dir.create(file.path(res_path,"survival_cancers"))
  }
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_PFS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_OS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
}

# 2.by peak sample score analysis ----------------------------------------
# by peak means in /project/huff/huff/github/immune-cp/Exp_pattern/ICP_exp_between_tumor_normal_group_samples.R, use peak value of gene expression density in normal samples to compare with tumor sample, by which group tumor samples into immune hot and cold.
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern/score_from.tumor_class_by_T_N.by_peak")

# distribution ------------------------------------------------------------
# 2.by peak sample score
tumor_class_by_T_N.by_peak.immunityScore %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(hot_per,0.5)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,mid) %>%
  dplyr::arrange(mid) %>%
  unique()-> cancer_rank.by_peak

tumor_class_by_T_N.by_peak.immunityScore %>%
  dplyr::filter(!is.na(hot_per)) %>%
  ggplot(aes(x=cancer_types,y=hot_per)) +
  # geom_jitter(size=0.5,width = 0.1) +
  scale_x_discrete(limits= cancer_rank.by_peak$cancer_types) +
  geom_violin(aes(fill = cancer_types),alpha=0.5) +
  rotate() +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12)
  )
ggsave(file.path(res_path,"score_distribution.png"),device = "png",width = 4,height = 6)
ggsave(file.path(res_path,"score_distribution.pdf"),device = "pdf",width = 4,height = 6)

##### all sample together survival -----------------------------------
# 2.by peak sample score ----
color_list <- c( "#CD2626","#00B2EE",c("#CDAD00"))
group <- c("High","Low","Mid")


tumor_class_by_T_N.by_peak.immunityScore %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.by_peak.immunityScore.clinical

sur_name <- paste("Score_0.5_Tumor_5year_OS")
tumor_class_by_T_N.by_peak.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years OS",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_OS")
tumor_class_by_T_N.by_peak.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,sur_name,res_path,3,4,0.7,0.9)


sur_name <- paste("Score_0.75-0.25_Tumor_allyear_OS")
tumor_class_by_T_N.by_peak.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_PFS")
tumor_class_by_T_N.by_peak.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_PFS_from")
tumor_class_by_T_N.by_peak.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,sur_name,res_path,3,4,0.7,0.9)

##### cacner specific survival -----------------------------------
# 2.by peak sample score

tumor_class_by_T_N.by_peak.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,group,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="group",value="n") %>%
  dplyr::filter(High >5) %>%
  .$cancer_types -> cancer_to_do_survival.by_peak

for(cancer in cancer_to_do_survival.by_peak){
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_PFS",sep="-")
  tumor_class_by_T_N.by_peak.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    # dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
    # dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,sur_name,file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_OS",sep="-")
  tumor_class_by_T_N.by_peak.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    # dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
    # dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,sur_name,file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
}

# 3. by mean sample score analysis ----------------------------------------
# by mean means in /project/huff/huff/github/immune-cp/Exp_pattern/ICP_exp_between_tumor_normal_group_samples.R, use mean value of gene expression density in normal samples to compare with tumor sample, by which group tumor samples into immune hot and cold.
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern/score_from.tumor_class_by_T_N.by_mean")

# distribution ------------------------------------------------------------
# 2.by peak sample score
tumor_class_by_T_N.by_mean.immunityScore %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(hot_per,0.5)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,mid) %>%
  dplyr::arrange(mid) %>%
  unique()-> cancer_rank.by_mean

tumor_class_by_T_N.by_mean.immunityScore %>%
  dplyr::filter(!is.na(hot_per)) %>%
  ggplot(aes(x=cancer_types,y=hot_per)) +
  # geom_jitter(size=0.5,width = 0.1) +
  scale_x_discrete(limits= cancer_rank.by_mean$cancer_types) +
  geom_violin(aes(fill = cancer_types),alpha=0.5) +
  rotate() +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12)
  )
ggsave(file.path(res_path,"score_distribution.png"),device = "png",width = 4,height = 6)
ggsave(file.path(res_path,"score_distribution.pdf"),device = "pdf",width = 4,height = 6)

##### all sample together survival -----------------------------------
# 2.by mean sample score ----
color_list <- c( "#CD2626","#00B2EE",c("#CDAD00"))
group <- c("High","Low","Mid")


tumor_class_by_T_N.by_mean.immunityScore %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.by_mean.immunityScore.clinical

sur_name <- paste("Score_0.5_Tumor_5year_OS")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years OS",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_OS")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,sur_name,res_path,3,4,0.7,0.9)


sur_name <- paste("Score_0.75-0.25_Tumor_allyear_OS")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_PFS")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_PFS_from")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,sur_name,res_path,3,4,0.7,0.9)

##### cacner specific survival -----------------------------------
# 2.by mean sample score

tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,group,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="group",value="n") %>%
  dplyr::filter(High >5) %>%
  .$cancer_types -> cancer_to_do_survival.by_mean

for(cancer in cancer_to_do_survival.by_mean){
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_PFS",sep="-")
  tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
    dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,sur_name,file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_OS",sep="-")
  tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
    dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,sur_name,file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
}
