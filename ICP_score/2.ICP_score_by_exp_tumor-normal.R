library(magrittr)
library(ggplot2)
library(ggpubr)

# data path 
basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")

tumor_class_by_T_N.only_paired.count <- readr::read_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio","tumor_class_by_T_N.only_paired.by_geneCounts","tumor_class_by_T_N.only_paired"))
tumor_class_by_T_N.only_paired.score <- readr::read_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio","tumor_class_by_T_N.only_paired.by_FCProduct","tumor_class_by_T_N.only_paired"))

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


# survival data -------------------------------------------------------
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

fun_clinical_test <- function(expr_clinical_ready, cancer_types){
  if(nrow(expr_clinical_ready %>% dplyr::filter(status!=0)) < 5){return(tibble::tibble())}
  print(cancer_types)
  
  
  # cox p
  .cox <- survival::coxph(survival::Surv(time, status) ~ group, data = expr_clinical_ready, na.action = na.exclude)
  summary(.cox) -> .z
  
  # kmp
  kmp <- 1 - pchisq(survival::survdiff(survival::Surv(time, status) ~ group, data = expr_clinical_ready, na.action = na.exclude)$chisq,df = length(levels(as.factor(expr_clinical_ready$group))) - 1)
  
  tibble::tibble(
    n = .z$n,
    hr = .z$conf.int[1],
    hr_l = .z$conf.int[3],
    hr_h = .z$conf.int[4],
    coxp = .z$waldtest[3],
    kmp = kmp) %>%
    dplyr::mutate(status = ifelse(hr > 1, "High_risk", "Low_risk"))
}

# 1.only paired gene percent score analysis ----------------------------------------
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio/tumor_class_by_T_N.only_paired.by_geneCounts/")

# 1.1.distribution ------------------------------------------------------------
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
tumor_class_by_T_N.only_paired.genePercent %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(n=dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cancer_types = paste(cancer_types,"(n=",n,")",sep="")) -> tumor_class_by_T_N.only_paired.genePercent.plot

tumor_class_by_T_N.only_paired.genePercent.plot %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(hot_per,0.5)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,mid) %>%
  dplyr::arrange(mid) %>%
  unique()-> cancer_rank

tumor_class_by_T_N.only_paired.genePercent.plot %>%
  dplyr::filter(!is.na(hot_per)) %>%
  ggplot(aes(x=cancer_types,y=hot_per)) +
  # geom_jitter(size=0.5,width = 0.1) +
  scale_x_discrete(limits= cancer_rank$cancer_types) +
  geom_boxplot(aes(fill = cancer_types),alpha=0.5) +
  # rotate() +
  labs(y = "Ratio of ICPs indicate high immunoplasticity") +
  my_theme +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(colour = "black",vjust = 0.5,hjust = 1, angle = 90)
  )
ggsave(file.path(res_path,"score_distribution.png"),device = "png",width = 6,height = 4)
ggsave(file.path(res_path,"score_distribution.pdf"),device = "pdf",width = 6,height = 4)

cancer_color <- readr::read_tsv(file.path(basic_path,"data/TCGA","02_pcc.tsv"))
cancer_color %>%
  dplyr::filter(cancer_types %in% unique(tumor_class_by_T_N.only_paired.genePercent$cancer_types)) -> cancer21_color

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
  my_theme +
  theme(
    legend.position = "none",
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(fill="white",color="black"),
    strip.text = element_text(color="black",size=12)
  )
ggsave(file.path(res_path,"score_distribution_in_cancers.png"),device = "png",height = 6,width = 8)
ggsave(file.path(res_path,"score_distribution_in_cancers.pdf"),device = "pdf",height = 6,width = 8)

##### 1.2.all sample together survival -----------------------------------
# 1.2.1.only paired sample score survival -----

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

##### 1.2.2.cancer specific survival -----------------------------------
# 1.2.2.1.group samples into 2groups by middle score ----
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

# 1.2.2.2.group samples into 3groups by quantile 0.75 and 0.25 score ----
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.75),"Low",group)) %>%
  dplyr::group_by(cancer_types,group) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,group,n) %>%
  unique() %>%
  tidyr::spread(key="group",value="n") %>%
  dplyr::filter(High >5) %>%
  .$cancer_types -> cancer_to_do_survival.only_paired

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE",c("#CDAD00")),
                             group=c("High","Low","Mid"))
for(cancer in cancer_to_do_survival.only_paired){
  if(!dir.exists(file.path(res_path,"survival_cancers.3groups"))){
    dir.create(file.path(res_path,"survival_cancers.3groups"))
  }
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_PFS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
    dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers.3groups"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_OS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::filter(!is.na(hot_per)) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.75),"High","Mid")) %>%
    dplyr::mutate(group = ifelse(hot_per<quantile(hot_per,0.25),"Low",group)) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers.3groups"),3,4,0.7,0.9)
}

# 1.3.kmp and coxp result table -------------------
### 1.3.1.PFS -----
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"2High","1Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>=10) %>%
  dplyr::ungroup() %>%
  dplyr::rename("score" = "hot_per") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(surv_test = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.only_paired.immunityScore.PFS_res
tumor_class_by_T_N.only_paired.immunityScore.PFS_res %>%
  readr::write_tsv(file.path(res_path,"PFS_COXph_KM_HR-bygroup.tsv")) # bygroup means treat score as catergray factor instead of continus variable, change it by changing cox model in fun_clinical_test
tumor_class_by_T_N.only_paired.immunityScore.PFS_res %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::filter(hr_h<10) %>%
  dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1)-> plot_ready
plot_ready %>% 
  ggplot(aes(y = hr, x = cancer_types, ymin=hr_l,ymax=hr_h)) +
  geom_pointrange(aes(color=cox_sig)) +
  scale_color_manual(values=c("red","black")) +
  geom_hline(aes(yintercept = 1)) +
  scale_size(name = "p-value") +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),
                     labels = c("1/16","1/8","1/4","1/2",1,2,3)) +
  coord_flip() +
  ggthemes::theme_gdocs() +
  theme(
    legend.position = "none",
    axis.line.y = element_line(color="black"),
    axis.text = element_text(color = "black",size=8),
    axis.title = element_text(color = "black",size=10),
    text = element_text(color = "black")
  ) +
  labs(y = "Hazard ratio (High vs. low immunoplasticity)", x = "Cancers",title = "Progression-free survival") -> p;p
ggsave(file.path(res_path,"PFS_COXph_HR-bygroup.png"),device = "png",width = 4,height = 4)
ggsave(file.path(res_path,"PFS_COXph_HR-bygroup.pdf"),device = "pdf",width = 4,height = 4)

### 1.3.2.OS ----
tumor_class_by_T_N.only_paired.immunityScore.clinical %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::filter(!is.na(hot_per)) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(hot_per>quantile(hot_per,0.5),"2High","1Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>=10) %>%
  dplyr::ungroup() %>%
  dplyr::rename("score" = "hot_per") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(surv_test = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.only_paired.immunityScore.OS_res

tumor_class_by_T_N.only_paired.immunityScore.OS_res %>%
  readr::write_tsv(file.path(res_path,"OS_COXph_KM-bygroup.tsv"))

tumor_class_by_T_N.only_paired.immunityScore.OS_res %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::filter(hr_h<10) %>%
  dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) -> plot_ready

plot_ready %>% 
  ggplot(aes(y = hr, x = cancer_types, ymin=hr_l,ymax=hr_h)) +
  geom_pointrange(aes(color=cox_sig)) +
  scale_color_manual(values=c("red","black")) +
  geom_hline(aes(yintercept = 1)) +
  scale_size(name = "p-value") +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),
                     labels = c("1/16","1/8","1/4","1/2",1,2,3)) +
  coord_flip() +
  ggthemes::theme_gdocs() +
  theme(
    legend.position = "none",
    axis.line.y = element_line(color="black"),
    axis.text = element_text(color = "black",size=8),
    axis.title = element_text(color = "black",size=10),
    text = element_text(color = "black")
  ) +
  labs(y = "Hazard ratio (High vs. low immunoplasticity)", x = "Cancers",title = "Overall survival") -> p;p
ggsave(file.path(res_path,"OS_COXph_HR-bygroup.png"),device = "png",width = 4,height = 4)
ggsave(file.path(res_path,"OS_COXph_HR-bygroup.pdf"),device = "pdf",width = 4,height = 4)

### 1.3.3. OS and PFS res combine--------
tumor_class_by_T_N.only_paired.immunityScore.PFS_res %>%
  dplyr::mutate(class = "PFS") %>%
  rbind(tumor_class_by_T_N.only_paired.immunityScore.OS_res %>%
          dplyr::mutate(class = "OS"))  %>%
  dplyr::select(cancer_types,hr,coxp,kmp,class) %>%
  tidyr::gather(-cancer_types,-hr,-class,key="type",value="p.value") %>%
  tidyr::unite("class_type",class,type,sep="_") -> plot_ready
fn_hrclass <- function(.x){
  if(.x<=0.125){
    class <-"0" 
  }else if(.x<=0.25){
    class <-"1" 
  }else if(.x>0.25 & .x<=0.5){
    class <-"2" 
  }else if(.x>0.5 & .x<1){
    class <-"3" 
  }else if(.x>=1 & .x<=2){
    class <-"4" 
  }else if(.x>2 & .x<=3){
    class <-"5" 
  }else if(.x>3){
    class <-"6" 
  }
  class
}
plot_ready %>%
  dplyr::mutate(HR_class = purrr::map(hr,fn_hrclass)) %>%
  tidyr::unnest() %>%
  dplyr::select(HR_class) %>%
  unique() %>%
  dplyr::inner_join(tibble::tibble(values=c("#00008B","#104E8B", "#1874CD", "#87CEFF", "#FFC0CB", "#FF6A6A", "#FF0000"),
                                   HR_class=c("0","1","2","3","4","5","6"),
                                   labels = c("<=1/8","1/8-1/4","1/4-1/2","1/2-1","1-2","2-3",">3")),by="HR_class") %>%
  dplyr::arrange(HR_class)-> color
plot_ready %>%
  dplyr::mutate(HR_class = purrr::map(hr,fn_hrclass)) %>%
  tidyr::unnest() %>%
  dplyr::arrange(HR_class) %>%
  dplyr::mutate(p = ifelse(p.value>0.05,"1","2")) %>%
  ggplot(aes(x=cancer_types,y=class_type)) +
  geom_point(aes(size=p,fill=HR_class),color="grey",shape=22) +
  scale_fill_manual(
    values = color$values,
    labels = color$labels
  ) +
  scale_size_discrete(
    labels = c("p>0.05","p<=0.05"),
    name = NULL
  ) +
  my_theme +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(vjust = 1, hjust = 1, size = 10,angle = 90),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black")
  )
ggsave(file.path(res_path,"survival_merge_res-bygroup.png"),device = "png",width = 4,height = 4)
ggsave(file.path(res_path,"survival_merge_res-bygroup.pdf"),device = "pdf",width = 4,height = 4)


# 2.only paired gene FC score analysis [FC in fantom * FC in TCGA]----------------------
# mean value of each genes FC in fantom * FC in TCGA
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio/tumor_class_by_T_N.only_paired.by_FCProduct")

# 2.1.distribution ------------------------------------------------------------
tumor_class_by_T_N.only_paired.score %>%
  dplyr::filter(score_mean!=0) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(score_mean,0.5)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,mid) %>%
  dplyr::arrange(mid) %>%
  unique() -> cancer_rank
library(ggbeeswarm)
tumor_class_by_T_N.only_paired.score %>%
  dplyr::filter(score_mean!=0) %>%
  ggplot(aes(x=cancer_types,y=score_mean)) +
  # geom_jitter(size=0.5,width = 0.1) +
  geom_jitter(aes(color=cancer_types),size=0.5,width=0.1) +
  scale_x_discrete(limits= cancer_rank$cancer_types) +
  geom_violin(aes(fill = cancer_types),alpha=0.5) +
  # rotate() +
  labs(y="Immunoplasticity score") +
  my_theme +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(colour = "black",size = 12,angle = 90, vjust = 0.5)
  )
ggsave(file.path(res_path,"score_distribution.png"),device = "png",width = 4,height = 6)
ggsave(file.path(res_path,"score_distribution.pdf"),device = "pdf",width = 4,height = 6)

cancer_color <- readr::read_tsv(file.path(basic_path,"data/TCGA","02_pcc.tsv"))
cancer_color %>%
  dplyr::filter(cancer_types %in% unique(tumor_class_by_T_N.only_paired.score$cancer_types)) -> cancer21_color

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
  ggtitle("Density of immunoplasticity score") +
  ylab("Density") +
  xlab("Immunoplasticity score") +
  my_theme +
  theme(
    legend.position = "none",
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(fill="white",color="black"),
    strip.text = element_text(color="black",size=12)
  )
ggsave(file.path(res_path,"score_distribution_in_cancers.png"),device = "png",height = 6,width = 8)
ggsave(file.path(res_path,"score_distribution_in_cancers.pdf"),device = "pdf",height = 6,width = 8)

##### 2.2.all sample together survival -----------------------------------
# 2.2.1.only paired sample score survival -----

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

##### 2.2.2.cancer specific survival -----------------------------------
# 2.2.2.1.middle value classify samples into 2 groups(high and low)----
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
  if(!dir.exists(file.path(res_path,"survival_cancers."))){
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

# 2.2.2.2.quantile value classify samples into 3 groups(high and low, mid)----
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_mean!=0) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.75),"High","Mid")) %>%  
  dplyr::mutate(group = ifelse(score_mean<quantile(score_mean,0.25),"Low",group)) %>%
  dplyr::group_by(cancer_types,group) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,group,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="group",value="n") %>%
  dplyr::filter(High >5) %>%
  .$cancer_types -> cancer_to_do_survival.only_paired

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE",c("#CDAD00")),
                             group=c("High","Low","Mid"))

for(cancer in cancer_to_do_survival.only_paired){
  if(!dir.exists(file.path(res_path,"survival_cancers.3groups"))){
    dir.create(file.path(res_path,"survival_cancers.3groups"))
  }
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_PFS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.Score.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.75),"High","Mid")) %>%  
    dplyr::mutate(group = ifelse(score_mean<quantile(score_mean,0.25),"Low",group)) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers.3groups"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_OS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.Score.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.75),"High","Mid")) %>%  
    dplyr::mutate(group = ifelse(score_mean<quantile(score_mean,0.25),"Low",group)) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers.3groups"),3,4,0.7,0.9)
}

# 2.3.kmp and coxp result table -----------
### 2.3.1.PFS-----
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"2High","1Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>=10) %>%
  dplyr::ungroup() %>%
  dplyr::rename("score" = "score_mean") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(surv_test = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.only_paired.Score.PFS_res
tumor_class_by_T_N.only_paired.Score.PFS_res %>%
  readr::write_tsv(file.path(res_path,"PFS_COXph_KM_HR-bygroup.tsv"))
tumor_class_by_T_N.only_paired.Score.PFS_res %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<=0.1,"1yes","2no")) %>%
  dplyr::filter(hr_h<10) %>%
  dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort))%>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1)-> plot_ready
plot_ready %>% 
  ggplot(aes(y = hr, x = cancer_types, ymin=hr_l,ymax=hr_h)) +
  geom_pointrange(aes(color=cox_sig)) +
  scale_color_manual(values=c("red","black")) +
  geom_hline(aes(yintercept = 1)) +
  scale_size(name = "p-value") +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),
                     labels = c("1/16","1/8","1/4","1/2",1,2,3)) +
  coord_flip() +
  ggthemes::theme_gdocs() +
  theme(
    legend.position = "none",
    axis.line.y = element_line(color="black"),
    axis.text = element_text(color = "black",size=8),
    axis.title = element_text(color = "black",size=10),
    text = element_text(color = "black")
  ) +
  labs(y = "Hazard ratio (High vs. low immunoplasticity score)", x = "Cancers",title = "Progression-free survival") -> p;p
ggsave(file.path(res_path,"PFS_COXph_HR-bygroup.png"),device = "png",width = 4,height = 4)
ggsave(file.path(res_path,"PFS_COXph_HR-bygroup.pdf"),device = "pdf",width = 4,height = 4)

### 2.3.2.OS ----
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"2High","1Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>=10) %>%
  dplyr::ungroup() %>%
  dplyr::rename("score" = "score_mean") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(surv_test = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.only_paired.Score.OS_res

tumor_class_by_T_N.only_paired.Score.OS_res %>%
  readr::write_tsv(file.path(res_path,"OS_COXph_KM-bygroup.tsv"))

tumor_class_by_T_N.only_paired.Score.OS_res %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::filter(hr_h<10) %>%
  dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort))%>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1)-> plot_ready

plot_ready %>% 
  ggplot(aes(y = hr, x = cancer_types, ymin=hr_l,ymax=hr_h)) +
  geom_pointrange(aes(color=cox_sig)) +
  scale_color_manual(values=c("red","black")) +
  geom_hline(aes(yintercept = 1)) +
  scale_size(name = "p-value") +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),
                     labels = c("1/16","1/8","1/4","1/2",1,2,3)) +
  coord_flip() +
  ggthemes::theme_gdocs() +
  theme(
    legend.position = "none",
    axis.line.y = element_line(color="black"),
    axis.text = element_text(color = "black",size=8),
    axis.title = element_text(color = "black",size=10),
    text = element_text(color = "black")
  ) +
  labs(y = "Hazard ratio (High vs. low immunoplasticity score)", x = "Cancers",title = "Overall survival") -> p;p
ggsave(file.path(res_path,"OS_COXph_HR-bygroup.png"),device = "png",width = 4,height = 4)
ggsave(file.path(res_path,"OS_COXph_HR-bygroup.pdf"),device = "pdf",width = 4,height = 4)

### 2.3.3.OS and PFS res combine -----------
tumor_class_by_T_N.only_paired.Score.PFS_res %>%
  dplyr::mutate(class = "PFS") %>%
  rbind(tumor_class_by_T_N.only_paired.Score.OS_res %>%
          dplyr::mutate(class = "OS"))  %>%
  dplyr::select(cancer_types,hr,coxp,kmp,class) %>%
  tidyr::gather(-cancer_types,-hr,-class,key="type",value="p.value") %>%
  tidyr::unite("class_type",class,type,sep="_") -> plot_ready

plot_ready %>%
  dplyr::mutate(HR_class = purrr::map(hr,fn_hrclass)) %>%
  tidyr::unnest() %>%
  dplyr::select(HR_class) %>%
  unique() %>%
  dplyr::inner_join(tibble::tibble(values=c("#00008B","#104E8B", "#1874CD", "#87CEFF", "#FFC0CB", "#FF6A6A", "#FF0000"),
                                   HR_class=c("0","1","2","3","4","5","6"),
                                   labels = c("<=1/8","1/8-1/4","1/4-1/2","1/2-1","1-2","2-3",">3")),by="HR_class") %>%
  dplyr::arrange(HR_class)-> color
plot_ready %>%
  dplyr::mutate(HR_class = purrr::map(hr,fn_hrclass)) %>%
  tidyr::unnest() %>%
  dplyr::select(HR_class) %>%
  unique() %>%
  dplyr::inner_join(tibble::tibble(values=c("#104E8B", "#1874CD", "#87CEFF", "#FFC0CB", "#FF6A6A", "#FF0000"),
                                   HR_class=c("1","2","3","4","5","6"),
                                   labels = c("<0.1","0.1-0.3","0.3-1","1-3","3-10",">10")),by="HR_class") %>%
  dplyr::arrange(HR_class)-> color
plot_ready %>%
  dplyr::mutate(HR_class = purrr::map(hr,fn_hrclass)) %>%
  tidyr::unnest() %>%
  dplyr::arrange(HR_class) %>%
  dplyr::mutate(p = ifelse(p.value>0.05,"1","2")) %>%
  ggplot(aes(x=cancer_types,y=class_type)) +
  geom_point(aes(size=p,fill=HR_class),color="grey",shape=22) +
  scale_fill_manual(
    values=color$values,
    breaks = color$HR_class,
    labels = color$labels
  ) +
  scale_size_discrete(
    labels = c("p>0.05","p<=0.05"),
    name = NULL
  ) +
  my_theme +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(vjust = 1, hjust = 1, size = 10,angle = 90),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black")
  )res_path
ggsave(file.path(,"survival_merge_res-bygroup.png"),device = "png",width = 4,height = 3)
ggsave(file.path(res_path,"survival_merge_res-bygroup.pdf"),device = "pdf",width = 4,height = 3)

# 3.only paired gene FC score analysis [sum of log2FC] ----------------------------------------

res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio/tumor_class_by_T_N.only_paired.by_logFCProduct_sum/")
tumor_class_by_T_N.only_paired.by_logFCProduct_sum <- readr::read_tsv(file.path(res_path,"tumor_class_by_T_N.only_paired.by_logFCProduct_sum"))

# 3.1.distribution ------------------------------------------------------------
tumor_class_by_T_N.only_paired.by_logFCProduct_sum %>%
  dplyr::filter(score_sum!=0) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(score_sum,0.5)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,mid) %>%
  dplyr::arrange(mid) %>%
  unique()-> cancer_rank
library(ggbeeswarm)
tumor_class_by_T_N.only_paired.by_logFCProduct_sum %>%
  dplyr::filter(score_sum!=0) %>%
  ggplot(aes(x=cancer_types,y=score_sum)) +
  # geom_jitter(size=0.5,width = 0.1) +
  geom_jitter(aes(color=cancer_types),size=0.5,width=0.1) +
  scale_x_discrete(limits= cancer_rank$cancer_types) +
  geom_violin(aes(fill = cancer_types),alpha=0.5) +
  # rotate() +
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
  dplyr::filter(cancer_types %in% unique(tumor_class_by_T_N.only_paired.by_logFCProduct_sum$cancer_types)) -> cancer21_color

tumor_class_by_T_N.only_paired.by_logFCProduct_sum %>%
  dplyr::filter(score_sum!=0) %>%
  ggdensity(x="score_sum",fill = "cancer_types",alpha = 0.5,palette = "jco")+
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

##### 3.2. survival -----------------------------------
# 3.2.1.all sample together -------------

tumor_class_by_T_N.only_paired.by_logFCProduct_sum %>%
  dplyr::rename("barcode" = "Participant") %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.only_paired.Score.clinical

color_list <- tibble::tibble(group=c("High","Low"),
                             color=c("red","blue"))
sur_name <- paste("Score_0.5_Tumor_5year_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_sum!=0) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"High","Low")) %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_sum!=0) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE",c("#CDAD00")),
                             group=c("High","Low","Mid"))

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_OS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_sum!=0) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(score_sum<quantile(score_sum,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_sum!=0) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(score_sum<quantile(score_sum,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_PFS_from_OnlyPaired")
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_sum!=0) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

##### 3.2.2.cancer specific survival -----------------------------------
# 3.2.2.1.middle value classify samples into 2 groups(high and low)----
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_sum!=0) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"High","Low")) %>%
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
  if(!dir.exists(file.path(res_path,"survival_cancers."))){
    dir.create(file.path(res_path,"survival_cancers"))
  }
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_PFS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.Score.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_OS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.Score.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
}

# 3.2.2.2.quantile value classify samples into 3 groups(high and low, mid)----
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::filter(score_sum!=0) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.75),"High","Mid")) %>%  
  dplyr::mutate(group = ifelse(score_sum<quantile(score_sum,0.25),"Low",group)) %>%
  dplyr::group_by(cancer_types,group) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,group,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="group",value="n") %>%
  dplyr::filter(High >5) %>%
  .$cancer_types -> cancer_to_do_survival.only_paired

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE",c("#CDAD00")),
                             group=c("High","Low","Mid"))

for(cancer in cancer_to_do_survival.only_paired){
  if(!dir.exists(file.path(res_path,"survival_cancers.3groups"))){
    dir.create(file.path(res_path,"survival_cancers.3groups"))
  }
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_PFS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.Score.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.75),"High","Mid")) %>%  
    dplyr::mutate(group = ifelse(score_sum<quantile(score_sum,0.25),"Low",group)) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers.3groups"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_OS_from_OnlyPaired",sep="-")
  tumor_class_by_T_N.only_paired.Score.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.75),"High","Mid")) %>%  
    dplyr::mutate(group = ifelse(score_sum<quantile(score_sum,0.25),"Low",group)) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers.3groups"),3,4,0.7,0.9)
}

# 3.3.kmp and coxp result table -----------
### 3.3.1.PFS-----
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"2High","1Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>=10) %>%
  dplyr::ungroup() %>%
  dplyr::rename("score" = "score_sum") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(surv_test = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.only_paired.Score.PFS_res
tumor_class_by_T_N.only_paired.Score.PFS_res %>%
  readr::write_tsv(file.path(res_path,"PFS_COXph_KM_HR--bygroup.tsv"))
tumor_class_by_T_N.only_paired.Score.PFS_res %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::filter(hr_h<10) %>%
  dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1)-> plot_ready
plot_ready %>% 
  ggplot(aes(y = hr, x = cancer_types, ymin=hr_l,ymax=hr_h)) +
  geom_pointrange(aes(color=cox_sig)) +
  scale_color_manual(name = "P value",values=c("red","black"),labels = c("<=0.1",">0.1")) +
  geom_hline(aes(yintercept = 1)) +
  scale_size(name = "p-value") +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),
                     labels = c("1/16","1/8","1/4","1/2",1,2,3)) +
  coord_flip() +
  ggthemes::theme_gdocs() +
  theme(
    # legend.position = "none",
    axis.line.y = element_line(color="black"),
    axis.text = element_text(color = "black",size=8),
    axis.title = element_text(color = "black",size=10),
    text = element_text(color = "black")
  ) +
  labs(y = "Hazard Ratio", x = "Cancers",title = "Progression-free survival") -> p;p
ggsave(file.path(res_path,"PFS_COXph_HR-bygroup.png"),device = "png",width = 4,height = 4)
ggsave(file.path(res_path,"PFS_COXph_HR-bygroup.pdf"),device = "pdf",width = 4,height = 4)

### 3.3.2.OS----
tumor_class_by_T_N.only_paired.Score.clinical %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"2High","1Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>=10) %>%
  dplyr::ungroup() %>%
  dplyr::rename("score" = "score_sum") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(surv_test = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.only_paired.Score.OS_res

tumor_class_by_T_N.only_paired.Score.OS_res %>%
  readr::write_tsv(file.path(res_path,"OS_COXph_KM-bygroup.tsv"))

tumor_class_by_T_N.only_paired.Score.OS_res %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::filter(hr_h<10) %>%
  dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) -> plot_ready

plot_ready %>% 
  ggplot(aes(y = hr, x = cancer_types, ymin=hr_l,ymax=hr_h)) +
  geom_pointrange(aes(color=cox_sig)) +
  scale_color_manual(name = "P value",values=c("red","black"),labels = c("<=0.1",">0.1")) +
  geom_hline(aes(yintercept = 1)) +
  scale_size(name = "p-value") +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),
                     labels = c("1/16","1/8","1/4","1/2",1,2,3)) +
  coord_flip() +
  ggthemes::theme_gdocs() +
  theme(
    # legend.position = "none",
    axis.line.y = element_line(color="black"),
    axis.text = element_text(color = "black",size=8),
    axis.title = element_text(color = "black",size=10),
    text = element_text(color = "black")
  ) +
  labs(y = "Hazard Ratio", x = "Cancers",title = "Overall survival") -> p;p
ggsave(file.path(res_path,"OS_COXph_HR.png"),device = "png",width = 4,height = 4)
ggsave(file.path(res_path,"OS_COXph_HR.pdf"),device = "pdf",width = 4,height = 4)

### 3.3.3.OS and PFS res combine------
tumor_class_by_T_N.only_paired.Score.PFS_res %>%
  dplyr::mutate(class = "PFS") %>%
  rbind(tumor_class_by_T_N.only_paired.Score.OS_res %>%
          dplyr::mutate(class = "OS"))  %>%
  dplyr::select(cancer_types,hr,coxp,kmp,class) %>%
  tidyr::gather(-cancer_types,-hr,-class,key="type",value="p.value") %>%
  tidyr::unite("class_type",class,type,sep="_") -> plot_ready

plot_ready %>%
  dplyr::mutate(HR_class = purrr::map(hr,fn_hrclass)) %>%
  tidyr::unnest() %>%
  dplyr::select(HR_class) %>%
  unique() %>%
  dplyr::inner_join(tibble::tibble(values=c("#00008B","#104E8B", "#1874CD", "#87CEFF", "#FFC0CB", "#FF6A6A", "#FF0000"),
                                   HR_class=c("0","1","2","3","4","5","6"),
                                   labels = c("<=1/8","1/8-1/4","1/4-1/2","1/2-1","1-2","2-3",">3")),by="HR_class") %>%
  dplyr::arrange(HR_class)-> color

plot_ready %>%
  dplyr::mutate(HR_class = purrr::map(hr,fn_hrclass)) %>%
  tidyr::unnest() %>%
  dplyr::arrange(HR_class) %>%
  dplyr::mutate(p = ifelse(p.value>0.05,"1","2")) %>%
  ggplot(aes(x=cancer_types,y=class_type)) +
  geom_point(aes(size=p,fill=HR_class),color="grey",shape=22) +
  scale_fill_manual(
    values=color$values,
    breaks = color$HR_class,
    labels = color$labels
  ) +
  scale_size_discrete(
    labels = c("p>0.05","p<=0.05"),
    name = NULL
  ) +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(vjust = 1, hjust = 1, size = 10,angle = 90),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black")
  )
ggsave(file.path(res_path,"survival_merge_res-bygroup.png"),device = "png",width = 4,height = 3)
ggsave(file.path(res_path,"survival_merge_res-bygroup.pdf"),device = "pdf",width = 4,height = 3)

# 4.by peak sample score analysis ----------------------------------------
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
# by peak sample score ----
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

# 5. by mean sample score analysis ----------------------------------------
# by mean means in /project/huff/huff/github/immune-cp/Exp_pattern/ICP_exp_between_tumor_normal_group_samples.R, use mean value of gene expression density in normal samples to compare with tumor sample, by which group tumor samples into immune hot and cold.
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern/tumor_class_by_T_N.all-by-mean.by_FCProduct/")
tumor_class_by_T_N.by_mean.immunityScore <- readr::read_tsv(file.path(res_path,"tumor_class_by_T_N.all_tumor.by_mean"))

# distribution ------------------------------------------------------------
# 2.by peak sample score
tumor_class_by_T_N.by_mean.immunityScore %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(score_mean,0.5)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,mid) %>%
  dplyr::arrange(mid) %>%
  unique()-> cancer_rank.by_mean

tumor_class_by_T_N.by_mean.immunityScore %>%
  ggplot(aes(x=cancer_types,y=score_mean)) +
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
# by mean sample score ----

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE",c("#CDAD00")),
                             group = c("High","Low","Mid"))


tumor_class_by_T_N.by_mean.immunityScore %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.by_mean.immunityScore.clinical

sur_name <- paste("Score_0.5_Tumor_5year_OS")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"High","Low")) %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_OS")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)


sur_name <- paste("Score_0.75-0.25_Tumor_allyear_OS")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(score_mean<quantile(score_mean,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_PFS")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(score_mean<quantile(score_mean,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_PFS_from")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::filter(!is.na(score_mean)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

##### cacner specific survival -----------------------------------
# 2.by mean sample score

tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"High","Low")) %>%
  dplyr::group_by(cancer_types,group) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,group,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="group",value="n") %>%
  dplyr::filter(High >10) %>%
  .$cancer_types -> cancer_to_do_survival.by_mean

for(cancer in cancer_to_do_survival.by_mean){
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_PFS",sep="-")
  if (!dir.exists(file.path(res_path,"survival_cancers.3groups"))) {
    dir.create(file.path(res_path,"survival_cancers.3groups"))
  }
  tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.75),"High","Mid")) %>%
    dplyr::mutate(group = ifelse(score_mean<quantile(score_mean,0.25),"Low",group)) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers.3groups"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_OS",sep="-")
  tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.75),"High","Mid")) %>%
    dplyr::mutate(group = ifelse(score_mean<quantile(score_mean,0.25),"Low",group)) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers.3groups"),3,4,0.7,0.9)
}

for(cancer in cancer_to_do_survival.by_mean){
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_PFS",sep="-")
  if (!dir.exists(file.path(res_path,"survival_cancers"))) {
    dir.create(file.path(res_path,"survival_cancers"))
  }
  tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_OS",sep="-")
  tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
}

# kmp and coxp result table -----------
### PFS
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"2High","1Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>=10) %>%
  dplyr::ungroup() %>%
  dplyr::rename("score" = "score_mean") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(surv_test = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.by_mean.Score.PFS_res
tumor_class_by_T_N.by_mean.Score.PFS_res %>%
  readr::write_tsv(file.path(res_path,"PFS_COXph_KM_HR-bygroup.tsv"))
tumor_class_by_T_N.by_mean.Score.PFS_res %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::filter(hr_h<10) %>%
  dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1)-> plot_ready
plot_ready %>% 
  ggplot(aes(y = hr, x = cancer_types, ymin=hr_l,ymax=hr_h)) +
  geom_pointrange(aes(color=cox_sig)) +
  scale_color_manual(name = "P value",values=c("red","black"),labels = c("<=0.1",">0.1")) +
  geom_hline(aes(yintercept = 1)) +
  scale_size(name = "p-value") +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),
                     labels = c("1/16","1/8","1/4","1/2",1,2,3)) +
  coord_flip() +
  ggthemes::theme_gdocs() +
  theme(
    # legend.position = "none",
    axis.line.y = element_line(color="black"),
    axis.text = element_text(color = "black",size=8),
    axis.title = element_text(color = "black",size=10),
    text = element_text(color = "black")
  ) +
  labs(y = "Hazard Ratio", x = "Cancers",title = "Progression-free survival") -> p;p
ggsave(file.path(res_path,"PFS_COXph_HR-bygroup.png"),device = "png",width = 4,height = 4)
ggsave(file.path(res_path,"PFS_COXph_HR-bygroup.pdf"),device = "pdf",width = 4,height = 4)

### OS
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_mean>quantile(score_mean,0.5),"2High","1Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>=10) %>%
  dplyr::ungroup() %>%
  dplyr::rename("score" = "score_mean") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(surv_test = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.by_mean.Score.OS_res

tumor_class_by_T_N.by_mean.Score.OS_res %>%
  readr::write_tsv(file.path(res_path,"OS_COXph_KM-bygroup.tsv"))

tumor_class_by_T_N.by_mean.Score.OS_res %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::filter(hr_h<10) %>%
  dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1)-> plot_ready

plot_ready %>% 
  ggplot(aes(y = hr, x = cancer_types, ymin=hr_l,ymax=hr_h)) +
  geom_pointrange(aes(color=cox_sig)) +
  geom_hline(aes(yintercept = 1)) +
  scale_size(name = "p-value") +
  scale_color_manual(name = "P value",values=c("red","black"),labels = c("<=0.1",">0.1")) +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),
                     labels = c("1/16","1/8","1/4","1/2",1,2,3)) +
  coord_flip() +
  ggthemes::theme_gdocs() +
  theme(
    # legend.position = "none",
    axis.line.y = element_line(color="black"),
    axis.text = element_text(color = "black",size=8),
    axis.title = element_text(color = "black",size=10),
    text = element_text(color = "black")
  ) +
  labs(y = "Hazard Ratio", x = "Cancers",title = "Overall survival") -> p;p
ggsave(file.path(res_path,"OS_COXph_HR-bygroup.png"),device = "png",width = 4,height = 4)
ggsave(file.path(res_path,"OS_COXph_HR-bygroup.pdf"),device = "pdf",width = 4,height = 4)

### OS and PFS res combine
tumor_class_by_T_N.by_mean.Score.PFS_res %>%
  dplyr::mutate(class = "PFS") %>%
  rbind(tumor_class_by_T_N.by_mean.Score.OS_res %>%
          dplyr::mutate(class = "OS"))  %>%
  dplyr::select(cancer_types,hr,coxp,kmp,class) %>%
  tidyr::gather(-cancer_types,-hr,-class,key="type",value="p.value") %>%
  tidyr::unite("class_type",class,type,sep="_") -> plot_ready

plot_ready %>%
  dplyr::mutate(HR_class = purrr::map(hr,fn_hrclass)) %>%
  tidyr::unnest() %>%
  dplyr::select(HR_class) %>%
  unique() %>%
  dplyr::inner_join(tibble::tibble(values=c("#00008B","#104E8B", "#1874CD", "#87CEFF", "#FFC0CB", "#FF6A6A", "#FF0000"),
                                   HR_class=c("0","1","2","3","4","5","6"),
                                   labels = c("<=1/8","1/8-1/4","1/4-1/2","1/2-1","1-2","2-3",">3")),by="HR_class") %>%
  dplyr::arrange(HR_class)-> color
plot_ready %>%
  dplyr::mutate(HR_class = purrr::map(hr,fn_hrclass)) %>%
  tidyr::unnest() %>%
  dplyr::arrange(HR_class) %>%
  dplyr::mutate(p = ifelse(p.value>0.05,"1","2")) %>%
  ggplot(aes(x=cancer_types,y=class_type)) +
  geom_point(aes(size=p,fill=HR_class),color="grey",shape=22) +
  scale_fill_manual(
    values=color$values,
    breaks = color$HR_class,
    labels = color$labels
  ) +
  scale_size_discrete(
    labels = c("p>0.05","p<=0.05"),
    name = NULL
  ) +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(vjust = 1, hjust = 1, size = 10,angle = 90),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black")
  )
ggsave(file.path(res_path,"survival_merge_res-bygroup.png"),device = "png",width = 4,height = 3)
ggsave(file.path(res_path,"survival_merge_res-bygroup.pdf"),device = "pdf",width = 4,height = 3)


# 6. by sum sample score analysis ----------------------------------------
# by mean means in /project/huff/huff/github/immune-cp/Exp_pattern/ICP_exp_between_tumor_normal_group_samples.R, use mean value of gene expression density in normal samples to compare with tumor sample, by which group tumor samples into immune hot and cold.
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio//tumor_class_by_T_N.all-by-mean.by_log2FCProduct_sum")
tumor_class_by_T_N.immunityScore.by_mean.by_log2FCProduct_sum <- readr::read_tsv(file.path(res_path,"tumor_class_by_T_N.all_tumor.by_mean.by_log2FCProduct_sum"))

# 6.1.istribution ------------------------------------------------------------
tumor_class_by_T_N.immunityScore.by_mean.by_log2FCProduct_sum %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(mid = quantile(score_sum,0.5)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cancer_types,mid) %>%
  dplyr::arrange(mid) %>%
  unique()-> cancer_rank.by_mean

tumor_class_by_T_N.immunityScore.by_mean.by_log2FCProduct_sum %>%
  ggplot(aes(x=cancer_types,y=score_sum)) +
  # geom_jitter(size=0.5,width = 0.1) +
  scale_x_discrete(limits= cancer_rank.by_mean$cancer_types) +
  geom_violin(aes(fill = cancer_types),alpha=0.5) +
  # rotate() +
  my_theme+
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    # axis.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text.x = element_text(colour = "black",angle = 90, vjust = 0.5)
  ) +
  labs(xlab="Score",ylab="Cancer types")
ggsave(file.path(res_path,"score_distribution.png"),device = "png",width = 4,height = 6)
ggsave(file.path(res_path,"score_distribution.pdf"),device = "pdf",width = 4,height = 6)

##### 6.2. survival -----------------------------------
# 6.2.1.all sample together ----

color_list <- tibble::tibble(color=c( "#CD2626","#00B2EE",c("#CDAD00")),
                             group = c("High","Low","Mid"))


tumor_class_by_T_N.immunityScore.by_mean.by_log2FCProduct_sum %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::inner_join(clinical,by="barcode") -> tumor_class_by_T_N.by_mean.immunityScore.clinical

sur_name <- paste("Score_0.5_Tumor_5year_OS")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"High","Low")) %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("5-years OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.5_Tumor_allyear_OS")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)


sur_name <- paste("Score_0.75-0.25_Tumor_allyear_OS")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(score_sum<quantile(score_sum,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("OS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

sur_name <- paste("Score_0.75-0.25_Tumor_allyear_PFS")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.75),"High","Mid")) %>%
  dplyr::mutate(group = ifelse(score_sum<quantile(score_sum,0.25),"Low",group)) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)


sur_name <- paste("Score_0.5_Tumor_allyear_PFS_from")
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::filter(!is.na(score_sum)) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"High","Low")) %>%
  # dplyr::filter(time<=1825) %>%
  fn_survival("PFS",color_list,"group",sur_name,"Survival in days",res_path,3,4,0.7,0.9)

##### 6.2.2.cacner specific survival -----------------------------------
# 2.by mean sample score

tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"High","Low")) %>%
  dplyr::group_by(cancer_types,group) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,group,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="group",value="n") %>%
  dplyr::filter(High >10) %>%
  .$cancer_types -> cancer_to_do_survival.by_mean

for(cancer in cancer_to_do_survival.by_mean){
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_PFS",sep="-")
  if (!dir.exists(file.path(res_path,"survival_cancers.3groups"))) {
    dir.create(file.path(res_path,"survival_cancers.3groups"))
  }
  tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.75),"High","Mid")) %>%
    dplyr::mutate(group = ifelse(score_sum<quantile(score_sum,0.25),"Low",group)) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers.3groups"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_OS",sep="-")
  tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.75),"High","Mid")) %>%
    dplyr::mutate(group = ifelse(score_sum<quantile(score_sum,0.25),"Low",group)) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers.3groups"),3,4,0.7,0.9)
}

for(cancer in cancer_to_do_survival.by_mean){
  sur_name <- paste(cancer,"Score_0.5_Tumor_allyear_PFS",sep="-")
  if (!dir.exists(file.path(res_path,"survival_cancers"))) {
    dir.create(file.path(res_path,"survival_cancers"))
  }
  tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="PFS.time","status"="PFS") %>%
    dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"PFS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
  
  sur_name <- paste(cancer,"Score_0.75-0.25_Tumor_allyear_OS",sep="-")
  tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
    dplyr::filter(cancer_types == cancer) %>%
    dplyr::rename("time"="OS","status"="Status") %>%
    dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"High","Low")) %>%
    # dplyr::filter(time<=1825) %>%
    fn_survival(paste(cancer,"OS"),color_list,"group",sur_name,"Survival in days",file.path(res_path,"survival_cancers"),3,4,0.7,0.9)
}

# 6.3.kmp and coxp result table -----------
### 6.3.1.PFS ----------
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"2High","1Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>=10) %>%
  dplyr::ungroup() %>%
  dplyr::rename("score" = "score_sum") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(surv_test = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.by_mean.Score.PFS_res
tumor_class_by_T_N.by_mean.Score.PFS_res %>%
  readr::write_tsv(file.path(res_path,"PFS_COXph_KM_HR-bygroup.tsv"))
tumor_class_by_T_N.by_mean.Score.PFS_res %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::filter(hr_h<10) %>%
  dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1)-> plot_ready
plot_ready %>% 
  ggplot(aes(y = hr, x = cancer_types, ymin=hr_l,ymax=hr_h)) +
  geom_pointrange(aes(color=cox_sig)) +
  geom_hline(aes(yintercept = 1)) +
  scale_size(name = "p-value") +
  scale_color_manual(name = "P value",values=c("red","black"),labels = c("<=0.1",">0.1")) +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),
                     labels = c("1/16","1/8","1/4","1/2",1,2,3)) +
  coord_flip() +
  ggthemes::theme_gdocs() +
  theme(
    # legend.position = "none",
    axis.line.y = element_line(color="black"),
    axis.text = element_text(color = "black",size=8),
    axis.title = element_text(color = "black",size=10),
    text = element_text(color = "black")
  ) +
  labs(y = "Hazard Ratio", x = "Cancers",title = "Progression-free survival") -> p;p
ggsave(file.path(res_path,"PFS_COXph_HR-bygroup.png"),device = "png",width = 4,height = 4)
ggsave(file.path(res_path,"PFS_COXph_HR-bygroup.pdf"),device = "pdf",width = 4,height = 4)

### 6.3.2.OS------
tumor_class_by_T_N.by_mean.immunityScore.clinical %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::rename("time"="OS","status"="Status") %>%
  dplyr::mutate(group = ifelse(score_sum>quantile(score_sum,0.5),"2High","1Low")) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>=10) %>%
  dplyr::ungroup() %>%
  dplyr::rename("score" = "score_sum") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(surv_test = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.by_mean.Score.OS_res

tumor_class_by_T_N.by_mean.Score.OS_res %>%
  readr::write_tsv(file.path(res_path,"OS_COXph_KM-bygroup.tsv"))

tumor_class_by_T_N.by_mean.Score.OS_res %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::filter(hr_h<10) %>%
  dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort))%>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1)-> plot_ready

plot_ready %>% 
  ggplot(aes(y = hr, x = cancer_types, ymin=hr_l,ymax=hr_h)) +
  geom_pointrange(aes(color=cox_sig)) +
  geom_hline(aes(yintercept = 1)) +
  scale_size(name = "p-value") +
  scale_color_manual(name = "P value",values=c("red","black"),labels = c("<=0.1",">0.1")) +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),
                     labels = c("1/16","1/8","1/4","1/2",1,2,3)) +
  coord_flip() +
  ggthemes::theme_gdocs() +
  theme(
    # legend.position = "none",
    axis.line.y = element_line(color="black"),
    axis.text = element_text(color = "black",size=8),
    axis.title = element_text(color = "black",size=10),
    text = element_text(color = "black")
  ) +
  labs(y = "Hazard Ratio", x = "Cancers",title = "Overall survival") -> p;p
ggsave(file.path(res_path,"OS_COXph_HR-bygroup.png"),device = "png",width = 4,height = 4)
ggsave(file.path(res_path,"OS_COXph_HR-bygroup.pdf"),device = "pdf",width = 4,height = 4)

### 6.3.3.OS and PFS res combine------------
tumor_class_by_T_N.by_mean.Score.PFS_res %>%
  dplyr::mutate(class = "PFS") %>%
  rbind(tumor_class_by_T_N.by_mean.Score.OS_res %>%
          dplyr::mutate(class = "OS"))  %>%
  dplyr::select(cancer_types,hr,coxp,kmp,class) %>%
  tidyr::gather(-cancer_types,-hr,-class,key="type",value="p.value") %>%
  tidyr::unite("class_type",class,type,sep="_") -> plot_ready


plot_ready %>%
  dplyr::mutate(HR_class = purrr::map(hr,fn_hrclass)) %>%
  tidyr::unnest() %>%
  dplyr::select(HR_class) %>%
  unique() %>%
  dplyr::inner_join(tibble::tibble(values=c("#00008B","#104E8B", "#1874CD", "#87CEFF", "#FFC0CB", "#FF6A6A", "#FF0000"),
                                   HR_class=c("0","1","2","3","4","5","6"),
                                   labels = c("<=1/8","1/8-1/4","1/4-1/2","1/2-1","1-2","2-3",">3")),by="HR_class") %>%
  dplyr::arrange(HR_class)-> color

plot_ready %>%
  dplyr::mutate(HR_class = purrr::map(hr,fn_hrclass)) %>%
  tidyr::unnest() %>%
  dplyr::arrange(HR_class) %>%
  dplyr::mutate(p = ifelse(p.value>0.05,"1","2")) %>%
  ggplot(aes(x=cancer_types,y=class_type)) +
  geom_point(aes(size=p,fill=HR_class),color="grey",shape=22) +
  scale_fill_manual(
    values=color$values,
    breaks = color$HR_class,
    labels = color$labels
  ) +
  scale_size_discrete(
    labels = c("p>0.05","p<=0.05"),
    name = NULL
  ) +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(vjust = 1, hjust = 1, size = 10,angle = 90),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black")
  )
ggsave(file.path(res_path,"survival_merge_res.png"),device = "png",width = 4,height = 3)
ggsave(file.path(res_path,"survival_merge_res.pdf"),device = "pdf",width = 4,height = 3)

