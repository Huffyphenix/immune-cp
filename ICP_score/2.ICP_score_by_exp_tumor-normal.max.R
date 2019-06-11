########################################
# ICP score from immunity hot/cold groups determined by comparing gene expression in tumor and max expression in normal.
########################################
library(magrittr)
library(dplyr)

# data path 
basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern/score_from.tumor_class_by_T_N.by_max")

# load data
tumor_class_by_T_N.by_max.class <-
  readr::read_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.by_max"))

cancer_color <- readr::read_tsv(file.path("/data/shiny-data/GSCALite","02_pcc.tsv"))

# get immunity score ---------------------------------------------------------------

tumor_class_by_T_N.by_max.class %>%
  dplyr::mutate(Immunity_hot = ifelse(is.na(Immunity_hot),0,Immunity_hot)) %>%
  dplyr::mutate(hot_per = Immunity_hot/(Immunity_hot+Immunity_cold)) -> tumor_class_by_T_N.by_max.immunityScore

# density peak of scores in each cancers
fn_density_peak <-function(.name,.x){
  print(.name)
  # get peak and secondary peak
  set.seed(42)
  .x <- .x$hot_per
  if(length(unique(as.vector(.x)))<3){
    tibble::tibble(peak.x = NA, peak.y = NA)
  }else{
    d1 <- density(.x)
    d1.ym <- which.max(d1$y)
    d1.ym.x <- d1$x[d1.ym]
    d1.ym.y <- d1$y[d1.ym]
    
    tibble::tibble(peak.x = d1.ym.x, peak.y = d1.ym.y)
  }
}
tumor_class_by_T_N.by_max.immunityScore %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(peak = purrr::map2(cancer_types,data,fn_density_peak)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.by_max.immunityScore.peak
  

# distribution of scores --------------------------------------------------
cancer_color %>%
  dplyr::filter(cancer_types %in% unique(tumor_class_by_T_N.by_max.immunityScore$cancer_types)) -> cancer21_color

tumor_class_by_T_N.by_max.immunityScore %>%
  ggdensity(x="hot_per",fill = "cancer_types",alpha = 0.5,palette = "jco")+
  facet_wrap(~cancer_types) +
  geom_vline(data=tumor_class_by_T_N.by_max.immunityScore.peak,aes(xintercept = peak.x),linetype = "dotted") +
  # geom_text(data=immune_lanscape_immunity.peak,aes(x=peak.x+0.2,y=peak.y+2,label=paste("(",signif(peak.x,2),",",signif(peak.y,2),")",sep=""))) +
  scale_fill_manual(
    values = cancer21_color$color,
    limits = cancer21_color$cancer_types
  ) +
  ggtitle("Density of Immune Activity Score in each cancers") +
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

# violin plot show distribution
tumor_class_by_T_N.by_max.immunityScore %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(hot_per,0.5)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,mid) %>%
  dplyr::arrange(mid) %>%
  unique()-> cancer_rank

tumor_class_by_T_N.by_max.immunityScore %>%
  dplyr::filter(!is.na(hot_per)) %>%
  ggplot(aes(x=cancer_types,y=hot_per)) +
  geom_jitter(size=0.5,width = 0.1,alpha = 0.5) +
  scale_x_discrete(limits= cancer_rank$cancer_types) +
  geom_violin(aes(fill = cancer_types),alpha=0.5) +
  scale_fill_manual(
    values = cancer21_color$color,
    limits = cancer21_color$cancer_types
  ) +
  scale_color_manual(
    values = cancer21_color$color,
    limits = cancer21_color$cancer_types
  ) +
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
ggsave(file.path(res_path,"score_distribution_in_cancers-violin.png"),device = "png",width = 4,height = 6)
ggsave(file.path(res_path,"score_distribution_in_cancers-violin.pdf"),device = "pdf",width = 4,height = 6)

##### all sample together survival -----------------------------------
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

source("/project/huff/huff/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

tumor_class_by_T_N.by_max.immunityScore %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.by_max.immunityScore.clinical


################# by score ----
color_list <- c( "#CD2626","#00B2EE",c("#CDAD00"))
group <- c("High","Low","Mid")
sur_name <- paste("Score_0.5_Tumor_5year_OS")
tumor_class_by_T_N.by_max.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::filter(time<=1825) %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  fn_survival("5-years OS",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_OS")
tumor_class_by_T_N.by_max.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,sur_name,res_path,3,4,0.7,0.9)


sur_name <- paste("Score_0.75-0.25_Tumor_allyear_OS")
tumor_class_by_T_N.by_max.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_PFS")
tumor_class_by_T_N.by_max.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,sur_name,res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_PFS")
tumor_class_by_T_N.by_max.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,sur_name,res_path,3,4,0.7,0.9)

################# by class ----
color_list <- c( "#00B2EE","#CD2626",c("#CDAD00"))
group <- c("Immunity_cold","Immunity_hot","Mid")

sur_name <- paste("ImmunityHC_0.5_Tumor_5year_OS")
tumor_class_by_T_N.by_max.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::filter(time<=1825) %>%
  dplyr::mutate(group = ifelse(Immunity_hot>Immunity_cold,"Immunity_hot","Immunity_cold")) %>%
  fn_survival("5-years OS",color_list,sur_name,res_path,3,4,0.7,0.9)
