library(magrittr)
library(ggplot2)

immune_path <- "/project/huff/huff/immune_checkpoint"

tumor_class_by_T_N.only_paired <- readr::read_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.only_paired"))

tumor_class_by_T_N.by_mean <- readr::read_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.by_mean"))

tumor_class_by_T_N.by_peak <- readr::read_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.by_peak"))

source("/project/huff/huff/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

# group distribution in cancers -------------------------------------------

# 1. only paired distribution -----
tumor_class_by_T_N.only_paired %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::group_by(cancer_types,class) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(N=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(`%(n/N)`=100*n/N) %>%
  dplyr::select(cancer_types,class,n,`%(n/N)`) %>%
  unique() %>%
  dplyr::mutate(`%(n/N)`=ifelse(class=="Immunity_cold",-`%(n/N)`,`%(n/N)`)) %>%
  ggplot(aes(x=class,y=cancer_types,fill=`%(n/N)`)) +
  geom_tile(color="grey", size=0.5,width = 0.9) +
  scale_fill_gradient2(low="#00B2EE",high = "#CD2626",mid="white") +
  geom_text(aes(label = n)) +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 10)
  )
ggsave(file.path(result_path,"sample_distribution_in_cancers.byOnlyPaired.pdf"),device="pdf",width = 5,height = 5)
ggsave(file.path(result_path,"sample_distribution_in_cancers.byOnlyPaired.png"),device="png",width = 5,height = 5)

# 2. by peak distribution -----
tumor_class_by_T_N.by_peak %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::group_by(cancer_types,class) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(N=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(`%(n/N)`=100*n/N) %>%
  dplyr::select(cancer_types,class,n,`%(n/N)`) %>%
  unique() %>%
  dplyr::mutate(`%(n/N)`=ifelse(class=="Immunity_cold",-`%(n/N)`,`%(n/N)`)) %>%
  ggplot(aes(x=class,y=cancer_types,fill=`%(n/N)`)) +
  geom_tile(color="grey", size=0.5,width = 0.9) +
  scale_fill_gradient2(low="#00B2EE",high = "#CD2626",mid="white") +
  geom_text(aes(label = n)) +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 10)
  )
ggsave(file.path(result_path,"sample_distribution_in_cancers.by_peak.pdf"),device="pdf",width = 4,height = 5)
ggsave(file.path(result_path,"sample_distribution_in_cancers.by_peak.png"),device="png",width = 4,height = 5)

# 3. by mean distribution -----
tumor_class_by_T_N.by_mean %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::group_by(cancer_types,class) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(N=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(`%(n/N)`=100*n/N) %>%
  dplyr::select(cancer_types,class,n,`%(n/N)`) %>%
  unique() %>%
  dplyr::mutate(`%(n/N)`=ifelse(class=="Immunity_cold",-`%(n/N)`,`%(n/N)`)) %>%
  ggplot(aes(x=class,y=cancer_types,fill=`%(n/N)`)) +
  geom_tile(color="grey", size=0.5,width = 0.9) +
  scale_fill_gradient2(low="#00B2EE",high = "#CD2626",mid="white") +
  geom_text(aes(label = n)) +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 10)
  )
ggsave(file.path(result_path,"sample_distribution_in_cancers.by_mean.pdf"),device="pdf",width = 4,height = 5)
ggsave(file.path(result_path,"sample_distribution_in_cancers.by_mean.png"),device="png",width = 4,height = 5)

# survival difference between groups --------------------------------------

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

color_list <- c( "#00B2EE","#CD2626")
group <- c("Immunity_cold","Immunity_hot")

res_path <- "/project/huff/huff/immune_checkpoint/result_20171025/ICP_exp_patthern/tumor_classify_by_Exp_site/survival"

# all cancers ------
# 1. all cancers survival only paired 
tumor_class_by_T_N.only_paired %>%
  dplyr::rename("barcode" = "Participant") %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.only_paired.clinical

sur_name <- paste("Hot_Cold_Tumor_5year_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.clinical %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::filter(class != "not clear") %>%
  dplyr::rename("group" = "class","time"="OS","status"="Status") %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years OS, Hot/Cold Tumor",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Hot_Cold_Tumor_allyear_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.clinical %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::filter(class != "not clear") %>%
  dplyr::rename("group" = "class","time"="OS","status"="Status") %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS, Hot/Cold Tumor",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Hot_Cold_Tumor_5year_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.clinical %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::filter(class != "not clear") %>%
  dplyr::rename("group" = "class","time"="PFS.time","status"="PFS") %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years PFS, Hot/Cold Tumor",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Hot_Cold_Tumor_allyear_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.clinical %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::filter(class != "not clear") %>%
  dplyr::rename("group" = "class","time"="PFS.time","status"="PFS") %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS, Hot/Cold Tumor",color_list,sur_name,res_path,3,4,0.7,0.9)

# 2. all cancers survival by mean exp ----
tumor_class_by_T_N.by_mean %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.by_mean.clinical

sur_name <- paste("Hot_Cold_Tumor_5year_OS_from_by-mean")
tumor_class_by_T_N.by_mean.clinical %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::filter(class != "not clear") %>%
  dplyr::rename("group" = "class","time"="OS","status"="Status") %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years OS, Hot/Cold Tumor",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Hot_Cold_Tumor_allyear_OS_from_by-mean")
tumor_class_by_T_N.by_mean.clinical %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::filter(class != "not clear") %>%
  dplyr::rename("group" = "class","time"="OS","status"="Status") %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS, Hot/Cold Tumor",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Hot_Cold_Tumor_5year_PFS_from_by-mean")
tumor_class_by_T_N.by_mean.clinical %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::filter(class != "not clear") %>%
  dplyr::rename("group" = "class","time"="PFS.time","status"="PFS") %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years PFS, Hot/Cold Tumor",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Hot_Cold_Tumor_allyear_PFS_from_by-mean")
tumor_class_by_T_N.by_mean.clinical %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::filter(class != "not clear") %>%
  dplyr::rename("group" = "class","time"="PFS.time","status"="PFS") %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS, Hot/Cold Tumor",color_list,sur_name,res_path,3,4,0.7,0.9)

# 3. all cancers survival by peak exp -----
tumor_class_by_T_N.by_peak %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.by_peak.clinical

sur_name <- paste("Hot_Cold_Tumor_5year_OS_from_by-peak")
tumor_class_by_T_N.by_peak.clinical %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::filter(class != "not clear") %>%
  dplyr::rename("group" = "class","time"="OS","status"="Status") %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years OS, Hot/Cold Tumor",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Hot_Cold_Tumor_allyear_OS_from_by-peak")
tumor_class_by_T_N.by_peak.clinical %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::filter(class != "not clear") %>%
  dplyr::rename("group" = "class","time"="OS","status"="Status") %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS, Hot/Cold Tumor",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Hot_Cold_Tumor_5year_PFS_from_by-peak")
tumor_class_by_T_N.by_peak.clinical %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::filter(class != "not clear") %>%
  dplyr::rename("group" = "class","time"="PFS.time","status"="PFS") %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years PFS, Hot/Cold Tumor",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Hot_Cold_Tumor_allyear_PFS_from_by-peak")
tumor_class_by_T_N.by_peak.clinical %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::filter(class != "not clear") %>%
  dplyr::rename("group" = "class","time"="PFS.time","status"="PFS") %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS, Hot/Cold Tumor",color_list,sur_name,res_path,3,4,0.7,0.9)

# cancers specific ------
# 1. cancer specific survival for only paired samples ----

tumor_class_by_T_N.only_paired.clinical %>%
  dplyr::group_by(cancer_types,class) %>%
  dplyr::mutate(n =n()) %>%
  dplyr::select(cancer_types,class,n) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::filter(!is.na(class) & class != "not clear") %>%
  tidyr::spread(key="class",value="n") %>%
  dplyr::filter(Immunity_cold >= 10 & Immunity_hot>=10) %>%
  .$cancer_types -> cancer_to_do_survival.onlypaired

color_list <- c( "#00B2EE","#CD2626")
group <- c("Immunity_cold","Immunity_hot")

for (cancer in cancer_to_do_survival.onlypaired) {
  tumor_class_by_T_N.only_paired.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(class) & class != "not clear") -> .data
  
  sur_name <- paste(cancer,"Hot_Cold_Tumor_allyear_PFS_from_OnlyPaired",sep="-")
  .data %>%
    dplyr::rename("group" = "class","time"="PFS.time","status"="PFS") %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS, Hot/Cold Tumor"),color_list,sur_name,file.path(res_path,"cancer_specific_onlypaired"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Hot_Cold_Tumor_5-years_PFS_from_OnlyPaired",sep="-")
  .data %>%
    dplyr::rename("group" = "class","time"="PFS.time","status"="PFS") %>%
    dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"5-years, Hot/Cold Tumor"),color_list,sur_name,file.path(res_path,"cancer_specific_onlypaired"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Hot_Cold_Tumor_5year_OS_from_OnlyPaired",sep="-")
  .data %>%
    dplyr::rename("group" = "class","time"="OS","status"="Status") %>%
    dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"5-years OS, Hot/Cold Tumor"),color_list,sur_name,file.path(res_path,"cancer_specific_onlypaired"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Hot_Cold_Tumor_allyear_OS_from_OnlyPaired",sep="-")
  .data %>%
    dplyr::filter(!is.na(class)) %>%
    dplyr::rename("group" = "class","time"="OS","status"="Status") %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS, Hot/Cold Tumor"),color_list,sur_name,file.path(res_path,"cancer_specific_onlypaired"),3,4,0.7,0.9)
}

# 2. cancer specific survival for by peak ----

tumor_class_by_T_N.by_peak.clinical %>%
  dplyr::group_by(cancer_types,class) %>%
  dplyr::mutate(n =n()) %>%
  dplyr::select(cancer_types,class,n) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::filter(!is.na(class) & class != "not clear") %>%
  tidyr::spread(key="class",value="n") %>%
  dplyr::filter(Immunity_cold >= 10 & Immunity_hot>=10) %>%
  .$cancer_types -> cancer_to_do_survival.by_peak

color_list <- c( "#00B2EE","#CD2626")
group <- c("Immunity_cold","Immunity_hot")

for (cancer in cancer_to_do_survival.by_peak) {
  tumor_class_by_T_N.by_peak.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(class) & class != "not clear") -> .data
  
  sur_name <- paste(cancer,"Hot_Cold_Tumor_allyear_PFS_from_by-peak",sep="-")
  .data %>%
    dplyr::rename("group" = "class","time"="PFS.time","status"="PFS") %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS, Hot/Cold Tumor"),color_list,sur_name,file.path(res_path,"cancer_specific_by-peak"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Hot_Cold_Tumor_5-years_PFS_from_by-peak",sep="-")
  .data %>%
    dplyr::rename("group" = "class","time"="PFS.time","status"="PFS") %>%
    dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"5-years, Hot/Cold Tumor"),color_list,sur_name,file.path(res_path,"cancer_specific_by-peak"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Hot_Cold_Tumor_5year_OS_from_by-peak",sep="-")
  .data %>%
    dplyr::rename("group" = "class","time"="OS","status"="Status") %>%
    dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"5-years OS, Hot/Cold Tumor"),color_list,sur_name,file.path(res_path,"cancer_specific_by-peak"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Hot_Cold_Tumor_allyear_OS_from_by-peak",sep="-")
  .data %>%
    dplyr::filter(!is.na(class)) %>%
    dplyr::rename("group" = "class","time"="OS","status"="Status") %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS, Hot/Cold Tumor"),color_list,sur_name,file.path(res_path,"cancer_specific_by-peak"),3,4,0.7,0.9)
}

# 2. cancer specific survival for by mean ----

tumor_class_by_T_N.by_mean.clinical %>%
  dplyr::group_by(cancer_types,class) %>%
  dplyr::mutate(n =n()) %>%
  dplyr::select(cancer_types,class,n) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::filter(!is.na(class) & class != "not clear") %>%
  tidyr::spread(key="class",value="n") %>%
  dplyr::filter(Immunity_cold >= 10 & Immunity_hot>=10) %>%
  .$cancer_types -> cancer_to_do_survival.by_mean

color_list <- c( "#00B2EE","#CD2626")
group <- c("Immunity_cold","Immunity_hot")

for (cancer in cancer_to_do_survival.by_mean) {
  tumor_class_by_T_N.by_mean.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(class) & class != "not clear") -> .data
  
  sur_name <- paste(cancer,"Hot_Cold_Tumor_allyear_PFS_from_by-mean",sep="-")
  .data %>%
    dplyr::rename("group" = "class","time"="PFS.time","status"="PFS") %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS, Hot/Cold Tumor"),color_list,sur_name,file.path(res_path,"cancer_specific_by-mean"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Hot_Cold_Tumor_5-years_PFS_from_by-mean",sep="-")
  .data %>%
    dplyr::rename("group" = "class","time"="PFS.time","status"="PFS") %>%
    dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"5-years, Hot/Cold Tumor"),color_list,sur_name,file.path(res_path,"cancer_specific_by-mean"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Hot_Cold_Tumor_5year_OS_from_by-mean",sep="-")
  .data %>%
    dplyr::rename("group" = "class","time"="OS","status"="Status") %>%
    dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"5-years OS, Hot/Cold Tumor"),color_list,sur_name,file.path(res_path,"cancer_specific_by-mean"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Hot_Cold_Tumor_allyear_OS_from_by-mean",sep="-")
  .data %>%
    dplyr::filter(!is.na(class)) %>%
    dplyr::rename("group" = "class","time"="OS","status"="Status") %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS, Hot/Cold Tumor"),color_list,sur_name,file.path(res_path,"cancer_specific_by-mean"),3,4,0.7,0.9)
}
