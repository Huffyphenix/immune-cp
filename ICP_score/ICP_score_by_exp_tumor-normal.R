library(magrittr)
library(ggplot2)
library(ggpubr)

immune_path <- "/project/huff/huff/immune_checkpoint"

tumor_class_by_T_N.only_paired <- readr::read_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.only_paired"))

tumor_class_by_T_N.by_mean <- readr::read_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.by_mean"))

tumor_class_by_T_N.by_peak <- readr::read_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.by_peak"))

source("/project/huff/huff/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

# get immunity score ---------------------------------------------------------------

tumor_class_by_T_N.only_paired %>%
  dplyr::select(-`<NA>`) %>%
  dplyr::mutate(hot_per = Immunity_hot/(Immunity_hot+Immunity_cold)) -> tumor_class_by_T_N.only_paired.immunityScore

tumor_class_by_T_N.by_mean %>%
  dplyr::mutate(hot_per = Immunity_hot/(Immunity_hot+Immunity_cold)) -> tumor_class_by_T_N.by_mean.immunityScore

tumor_class_by_T_N.by_peak %>%
  dplyr::mutate(hot_per = Immunity_hot/(Immunity_hot+Immunity_cold)) -> tumor_class_by_T_N.by_peak.immunityScore





# survival analysis -------------------------------------------------------
clinical_tcga <- readr::read_rds(file.path("/project/huff/huff/TCGA_survival/data","Pancan.Merge.clinical.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::select(-cancer_types) %>%
  unique() %>%
  dplyr::mutate(OS=as.numeric(OS),Status=as.numeric(Status),Age=as.numeric(Age)) %>%
  dplyr::group_by(barcode) %>%
  dplyr::mutate(OS= max(OS)) %>%
  dplyr::mutate(Status =  max(Status)) %>%
  dplyr::ungroup()

clinical <- readr::read_rds(file.path("/project/huff/huff/data//TCGA-survival-time/cell.2018.survival","TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::select(bcr_patient_barcode,PFS,PFS.time) %>%
  tidyr::drop_na() %>%
  dplyr::rename("barcode" = "bcr_patient_barcode") %>%
  unique()


# 1.only paired sample score analysis ----------------------------------------
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern/score_from.tumor_class_by_T_N.only_paired")

# distribution ------------------------------------------------------------
tumor_class_by_T_N.only_paired.immunityScore %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(hot_per,0.5)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,mid) %>%
  dplyr::arrange(mid) %>%
  unique()-> cancer_rank

tumor_class_by_T_N.only_paired.immunityScore %>%
  dplyr::filter(!is.na(hot_per)) %>%
  ggplot(aes(x=cancer_types,y=hot_per)) +
  geom_jitter(size=0.5,width = 0.1) +
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

##### all sample together survival -----------------------------------
# 1.only paired sample score survival -----

tumor_class_by_T_N.only_paired.immunityScore %>%
  dplyr::rename("barcode" = "Participant") %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.only_paired.immunityScore.clinical

sur_name <- paste("Score_0.5_Tumor_5year_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years OS",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,sur_name,res_path,3,4,0.7,0.9)

color_list <- c( "#CD2626","#00B2EE",c("#CDAD00"))
group <- c("High","Low","Mid")
sur_name <- paste("Score_0.75-0.25_Tumor_allyear_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,sur_name,res_path,3,4,0.7,0.9)

##### cacner specific survival -----------------------------------
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

for(cancer in cancer_to_do_survival.only_paired){
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_PFS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,sur_name,file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_OS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,sur_name,file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
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
