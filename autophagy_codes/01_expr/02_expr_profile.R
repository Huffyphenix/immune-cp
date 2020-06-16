###########################################################################
# expression profile of ICP genes in cancers ------------------------------

library(magrittr)
library(ggplot2)

# data path config
basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
tcga_path <- file.path(immune_path,"/data/TCGA_data")
gene_list_path <- file.path(immune_path,"checkpoint/20171021_checkpoint")
expr_path <- file.path(basic_path,"immune_checkpoint/result_20171025/expr_rds")
out_path <- file.path(basic_path,"immune_checkpoint/result_20171025")

# load(file.path(out_path,"e_6_exp_profile","e_6_exp_profile.rdata"))

# load data
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol <- as.character(gene_list$symbol)
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))
# ICP_expr_pattern <- readr::read_tsv(file.path(out_path,"ICP_exp_patthern","manual_edit_2_ICP_exp_pattern_in_immune_tumor_cell.tsv"))

ICP_expr_pattern <- readr::read_tsv(file.path(out_path,"ICP_exp_patthern-byMeanUQ/pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(entrez_ID,symbol,Exp_site,`log2FC(I/T).mean`,`log2FC(I/T).mid`) 

ICP_family <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/checkpoint/ICP_gene_family.txt"))

ICP_ligand_receptor_pair <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/checkpoint/ICP_gene_ligand_receptor_pairs.txt")) %>%
  dplyr::mutate(Recepter_pairs = paste("Pair",pairs)) %>%
  dplyr::select(-pairs)

fn_site_color <- function(.x){
  # print(.n)
  if(.x=="Mainly_exp_on_Tumor"){
    "red"
  }else if(.x=="Mainly_exp_on_Immune"){
    "darkgreen"
  }else if(.x=="Only_exp_on_Immune"){
    "Blue"
  }else if(.x=="Both_exp_on_Tumor_Immune"){
    c("darkorange")
  }else{
    "grey"
  }
}
gene_list %>%
  dplyr::left_join(ICP_expr_pattern,by="symbol") %>%
  dplyr::mutate(Exp_site=ifelse(is.na(Exp_site),"N",Exp_site)) %>%
  dplyr::mutate(site_col = purrr::map(Exp_site,fn_site_color)) %>%
  tidyr::unnest() %>%
  dplyr::left_join(ICP_family,by="symbol") %>%
  dplyr::left_join(ICP_ligand_receptor_pair,by="symbol") %>%
  dplyr::mutate(family = ifelse(is.na(family),"Other",family)) -> gene_list

gene_list %>%
  readr::write_tsv(file.path(gene_list_path,"ICPs_all_info_class-new.tsv"))

survival_path <- "/home/huff/project/data/TCGA-survival-time/cell.2018.survival"
survival_data <- readr::read_rds(file.path(survival_path, "TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>%
  dplyr::rename("cancer_types" = "type")

clinical_tcga <- readr::read_rds(file.path(basic_path,"TCGA_survival/data","Pancan.Merge.clinical-OS-Age-stage.rds.gz")) %>%
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

# calculation of average expression of genes of samples --------------------------------
# 3 groups by function: activate, inhibit and two sides -----
# 样本中在各基因的平均表达量
MHC_class_genes <- gene_list %>%
  dplyr::filter(family %in% c("MHC class I","MHC class II")) %>%
  .$symbol

fn_average <- function(.data){
  .data %>%
    dplyr::filter(!symbol %in%MHC_class_genes) %>%
    tidyr::gather(-symbol,-entrez_id,key="barcode",value="expr") %>%
    dplyr::left_join(gene_list,by="symbol") %>%
    dplyr::filter(!is.na(expr)) %>%
    dplyr::filter(substr(barcode,14,14)==0) %>%
    dplyr::group_by(barcode,functionWithImmune) %>%
    dplyr::mutate(average_exp = mean(expr)) %>%
    dplyr::select(barcode,functionWithImmune,average_exp) %>%
    unique() %>%
    dplyr::ungroup()
}

gene_list_expr %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(ICP_mean_expr = purrr::map(filter_expr,fn_average)) %>%
  dplyr::select(-filter_expr) %>%
  dplyr::ungroup() -> ICP_mean_expr_in_cancers
ICP_mean_expr_in_cancers %>%
  readr::write_rds(file.path(out_path,"e_6_exp_profile","ICP_mean_expr_in_cancers.by_functionalRole-noMHC.rds.gz"),compress = "gz")

ICP_mean_expr_in_cancers %>%
  tidyr::unnest() %>%
  dplyr::mutate(average_exp = log2(average_exp)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(sum = quantile(average_exp,0.5)) %>%
  dplyr::select(cancer_types,sum) %>%
  unique() %>%
  dplyr::arrange(sum) -> cancer_rank
ICP_mean_expr_in_cancers %>%
  tidyr::unnest() %>%
  dplyr::mutate(average_exp = log2(average_exp)) %>%
  ggplot(aes(x=cancer_types,y=average_exp)) +
  geom_violin(aes(fill=functionWithImmune),alpha=0.5) +
  facet_wrap(~functionWithImmune) +
  coord_flip() +
  scale_x_discrete(limits= cancer_rank$cancer_types) +
  scale_fill_manual(values = c( "#1C86EE", "#EE3B3B","#EE7600")) +
  theme_bw() +
  xlab("Cancer Types") +
  ylab("log2 (Expression)") +
  theme(
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(size = 12,color = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black"),
    panel.grid = element_line(linetype = "dashed")
  )
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_cancers.pdf"),device = "pdf", width = 6,height = 6)  
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_cancers.png"),device = "png",width = 6,height = 6)

# class by cancers
ICP_mean_expr_in_cancers %>%
  tidyr::unnest() %>%
  dplyr::mutate(average_exp = log2(average_exp)) %>%
  ggplot(aes(x=functionWithImmune,y=average_exp)) +
  geom_violin(aes(fill=cancer_types),alpha=0.5, size = 0.25) +
  facet_wrap(~cancer_types) +
  coord_flip() +
  # scale_x_discrete(limits= cancer_rank$cancer_types) +
  # scale_fill_manual(values = c( "#1C86EE", "#EE3B3B","#EE7600")) +
  theme_bw() +
  xlab("Funtional class of ICPs") +
  ylab("log2 (Expression)") +
  theme(
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(size = 12,color = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black"),
    panel.grid = element_blank()
  )
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_functionRoles-bycancers-noMHC.pdf"),device = "pdf", width = 8,height = 6)
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_functionRoles-bycancer-noMHC.png"),device = "png",width = 8,height =6)

# survival analysis
fun_clinical_test <- function(expr_clinical_ready, cancer_types){
  if(nrow(expr_clinical_ready %>% dplyr::filter(status!=0)) < 5){return(tibble::tibble())}
  print(cancer_types)
  
  # cox p
  .cox <- survival::coxph(survival::Surv(time, status) ~ group, data = expr_clinical_ready, na.action = na.exclude)
  summary(.cox) -> .z
  
  # cox p
  .cox <- survival::coxph(survival::Surv(time, status) ~ average_exp, data = expr_clinical_ready, na.action = na.exclude)
  summary(.cox) -> .z.continus
  
  # kmp
  kmp <- 1 - pchisq(survival::survdiff(survival::Surv(time, status) ~ group, data = expr_clinical_ready, na.action = na.exclude)$chisq,df = length(levels(as.factor(expr_clinical_ready$group))) - 1)
  
  tibble::tibble(
    n = .z$n,
    hr = .z$conf.int[1],
    hr_l = .z$conf.int[3],
    hr_h = .z$conf.int[4],
    coxp = .z$waldtest[3],
    hr.c = .z.continus$conf.int[1],
    hr.c.l = .z.continus$conf.int[3],
    hr.c.h = .z.continus$conf.int[4],
    hr.c.coxp = .z.continus$waldtest[3],
    kmp = kmp) %>%
    dplyr::mutate(status = ifelse(hr > 1, "High_risk", "Low_risk"))
}
## PFS
ICP_mean_expr_in_cancers %>%
  tidyr::unnest() %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical,by="barcode") %>%
  dplyr::rename("time" = "PFS.time", "status" = "PFS") %>%
  dplyr::group_by(cancer_types,functionWithImmune) %>%
  dplyr::mutate(group = ifelse(average_exp >= quantile(average_exp,0.5),"2High","1Low")) %>%
  dplyr::ungroup() %>%
  tidyr::nest(-cancer_types,-functionWithImmune) %>%
  dplyr::mutate(survival = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ICP_mean_expr_in_cancers.byfunction.PFS

ICP_mean_expr_in_cancers.byfunction.PFS  %>%
  readr::write_tsv(file.path(out_path,"e_6_exp_profile","survival","PFS_survival.ICP_mean_expr_in_cancers.byfunction-noMHC.tsv"))
## OS
ICP_mean_expr_in_cancers %>%
  tidyr::unnest() %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::rename("time" = "OS", "status" = "Status") %>%
  dplyr::group_by(cancer_types,functionWithImmune) %>%
  dplyr::mutate(group = ifelse(average_exp >= quantile(average_exp,0.5),"2High","1Low")) %>%
  dplyr::ungroup() %>%
  tidyr::nest(-cancer_types,-functionWithImmune) %>%
  dplyr::mutate(survival = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ICP_mean_expr_in_cancers.byfunction.OS
ICP_mean_expr_in_cancers.byfunction.OS  %>%
  readr::write_tsv(file.path(out_path,"e_6_exp_profile","survival","OS_survival.ICP_mean_expr_in_cancers.byfunction-noMHC.tsv"))

## survival plot, cox
# plot function
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
fn_cox_plot <- function(functionWithImmune,data,filename,hr,hr_l,hr_h,title,facet, dir,w=4,h=4){
  data %>% 
    dplyr::mutate(functionWithImmune=functionWithImmune) %>%
    dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
    ggplot(aes(y = hr, x = cancer_types, ymin=hr_l,ymax=hr_h)) +
    geom_pointrange(aes(color=cox_sig),size=0.5) +
    scale_color_manual(values=c("red","black")) +
    geom_hline(aes(yintercept = 1), linetype =2) +
    scale_size(name = "p-value") +
    scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),
                       labels = c("1/16","1/8","1/4","1/2",1,2,3)) +
    facet_wrap(as.formula(facet)) +
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

fn_cox_plot.all <- function(data,filename,hr_sum,hr,hr_l,hr_h,title,facet, dir,w=4,h=4){
  # data %>%
  #   dplyr::arrange(hr_sum) %>%
  #   .$cancer_types %>% unique() -> cancer_rank.cox
  # data <- within(data,cancer_types <- factor(cancer_types,labels = cancer_rank.cox))
  data %>% 
    dplyr::mutate(functionWithImmune=functionWithImmune) %>%
    dplyr::mutate(cancer_types = reorder(cancer_types,hr_sum,mean)) %>%
    ggplot(aes(y = hr, x = cancer_types, ymin=hr_l,ymax=hr_h)) +
    geom_pointrange(aes(color=cox_sig),size=0.5) +
    scale_color_manual(values=c("red","black")) +
    geom_hline(aes(yintercept = 1), linetype =2) +
    scale_size(name = "p-value") +
    scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),
                       labels = c("1/16","1/8","1/4","1/2",1,2,3)) +
    facet_wrap(as.formula(facet)) +
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
# PFS, cox, group
ICP_mean_expr_in_cancers.byfunction.PFS %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  tidyr::nest(-functionWithImmune) %>%
  dplyr::mutate(data = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort))
  })) %>%
  dplyr::mutate(filename = paste("Meanexp.COX_PFS.by-functionRole-noMHC",functionWithImmune,sep=".")) %>%
  purrr::pwalk(.f=fn_cox_plot,hr="hr",hr_l="hr_l",hr_h="hr_h",title="Progression-free survival",facet="~ functionWithImmune",dir = "survival",w = 4, h = 6)

ICP_mean_expr_in_cancers.byfunction.PFS %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(hr_sum = sum(hr)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  # tidyr::nest(-functionWithImmune) %>%
  # dplyr::mutate(data = purrr::map(data,.f=function(.x){
  #   .x %>%
  #     dplyr::mutate(cancer_types = reorder(cancer_types,hr_sum,mean))
  # })) %>%
  fn_cox_plot.all(filename="Meanexp.COX_PFS.by-functionRole.all-noMHC",hr="hr",hr_l="hr_l",hr_h="hr_h",title="Progression-free survival",facet="~ functionWithImmune",dir = "survival",w = 4, h = 3)

# PFS, cox, continus
ICP_mean_expr_in_cancers.byfunction.PFS %>% 
  dplyr::mutate(cox_sig = ifelse(hr.c.coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr.c=log2(hr.c)+1,hr.c.l=log2(hr.c.l)+1,hr.c.h=log2(hr.c.h)+1) %>%
  tidyr::nest(-functionWithImmune) %>%
  dplyr::mutate(data = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::mutate(cancer_types = reorder(cancer_types,hr.c,sort))
  })) %>%
  dplyr::mutate(filename = paste("Meanexp.COX_PFS.by-functionRole.continus",functionWithImmune,sep=".")) %>%
  purrr::pwalk(.f=fn_cox_plot,hr="hr.c",hr_l="hr.c.l",hr_h="hr.c.h",title="Progression-free survival",facet="~ functionWithImmune",dir = "survival",w = 4, h = 6)

ICP_mean_expr_in_cancers.byfunction.PFS %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr.c=log2(hr.c)+1,hr.c.l=log2(hr.c.l)+1,hr.c.h=log2(hr.c.h)+1) %>%
  fn_cox_plot.all(filename="Meanexp.COX_PFS.by-functionRole.all.continus",hr="hr.c",hr_l="hr.c.l",hr_h="hr.c.h",title="Progression-free survival",facet="~ functionWithImmune",dir = "survival",w = 8, h = 6)

# OS, cox, group
ICP_mean_expr_in_cancers.byfunction.OS %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  tidyr::nest(-functionWithImmune) %>%
  dplyr::mutate(data = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort))
  })) %>%
  dplyr::mutate(filename = paste("Meanexp.COX_OS.by-functionRole",functionWithImmune,sep=".")) %>%
  purrr::pwalk(.f=fn_cox_plot,hr="hr",hr_l="hr_l",hr_h="hr_h",title="Overall survival",facet="~ functionWithImmune",dir = "survival",w = 4, h = 6)

ICP_mean_expr_in_cancers.byfunction.OS %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(hr_sum = sum(hr)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  fn_cox_plot.all(filename="Meanexp.COX_OS.by-functionRole.all",hr="hr",hr_l="hr_l",hr_h="hr_h",title="Overall survival",facet="~ functionWithImmune",dir = "survival",w = 8, h = 6)

# OS, cox, continus
ICP_mean_expr_in_cancers.byfunction.OS %>% 
  dplyr::mutate(cox_sig = ifelse(hr.c.coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr.c=log2(hr.c)+1,hr.c.l=log2(hr.c.l)+1,hr.c.h=log2(hr.c.h)+1) %>%
  tidyr::nest(-functionWithImmune) %>%
  dplyr::mutate(data = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::mutate(cancer_types = reorder(cancer_types,hr.c,sort))
  })) %>%
  dplyr::mutate(filename = paste("Meanexp.COX_OS.by-functionRole.continus",functionWithImmune,sep=".")) %>%
  purrr::pwalk(.f=fn_cox_plot,hr="hr.c",hr_l="hr.c.l",hr_h="hr.c.h",title="Overall survival",facet="~ functionWithImmune",dir = "survival",w= 4, h = 6)

ICP_mean_expr_in_cancers.byfunction.OS %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr.c=log2(hr.c)+1,hr.c.l=log2(hr.c.l)+1,hr.c.h=log2(hr.c.h)+1) %>%
  fn_cox_plot.all(filename="Meanexp.COX_OS.by-functionRole.all.continus",hr="hr.c",hr_l="hr.c.l",hr_h="hr.c.h",title="Overall survival",facet="~ functionWithImmune",dir = "survival",w = 8, h = 6)

# expression in 3 groups by function: exp site Mainly_Tumor, Mainly_Immune, both -----
# 样本中在各基因的平均表达量

fn_average.by_expsite <- function(.data){
  .data %>%
    dplyr::filter(!symbol %in% MHC_class_genes) %>%
    tidyr::gather(-symbol,-entrez_id,key="barcode",value="expr") %>%
    dplyr::left_join(gene_list,by="symbol") %>%
    dplyr::filter(!is.na(expr)) %>%
    dplyr::filter(substr(barcode,14,14)==0) %>%
    dplyr::filter(Exp_site %in% c("Immune and tumor cell almost" ,"Immune cell dominate","Tumor cell dominate" )) %>%
    # dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Only_exp_on_Immune","Mainly_exp_on_Immune"),"Mainly_exp_on_Immune","Mainly_exp_on_Tumor")) %>%
    # dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Both_exp_on_Tumor_Immune"),"Both_exp_on_Tumor_Immune",Exp_site.1)) %>%
    # dplyr::mutate(Exp_site = Exp_site.1) %>%
    dplyr::group_by(barcode,Exp_site) %>%
    dplyr::mutate(average_exp = mean(expr)) %>%
    dplyr::select(barcode,Exp_site,average_exp) %>%
    unique() %>%
    dplyr::ungroup()
}

gene_list_expr %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(ICP_mean_expr = purrr::map(filter_expr,fn_average.by_expsite)) %>%
  dplyr::select(-filter_expr) %>%
  dplyr::ungroup() -> ICP_mean_expr_in_cancers.by_expsite

ICP_mean_expr_in_cancers.by_expsite %>%
  readr::write_rds(file.path(out_path,"e_6_exp_profile","ICP_mean_expr_in_cancers.by_expsite-noMHC.rds.gz"),compress = "gz") 

ICP_mean_expr_in_cancers.by_expsite %>%
  tidyr::unnest() %>%
  dplyr::mutate(average_exp = log2(average_exp)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(sum = quantile(average_exp,0.5)) %>%
  dplyr::select(cancer_types,sum) %>%
  unique() %>%
  dplyr::arrange(sum) -> cancer_rank.by_expsite

ICP_mean_expr_in_cancers.by_expsite %>%
  tidyr::unnest() %>%
  dplyr::filter(Exp_site!="N") %>%
  dplyr::mutate(average_exp = log2(average_exp)) %>%
  ggplot(aes(x=cancer_types,y=average_exp)) +
  geom_violin(aes(fill=Exp_site),alpha=0.5) +
  facet_wrap(~Exp_site) +
  coord_flip() +
  scale_x_discrete(limits= cancer_rank$cancer_types) +
  scale_fill_manual(values = c("#9A32CD", "Blue" ,"red","pink","orange","grey")) +
  theme_bw() +
  xlab("Cancer Types") +
  ylab("log2 (Expression)") +
  theme(
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(size = 12,color = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black"),
    panel.grid = element_line(linetype = "dashed")
  )
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_cancers-by_expsite.pdf"),device = "pdf", width = 6,height = 6)  
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_cancers-by_expsite.png"),device = "png",width = 6,height = 6)  

# class by cancers
ICP_mean_expr_in_cancers.by_expsite %>%
  tidyr::unnest() %>%
  dplyr::filter(!Exp_site %in% c("N","Not_sure")) -> plot_ready

plot_ready %>%
  dplyr::group_by(cancer_types, Exp_site) %>%
  dplyr::mutate(mean = mean(average_exp)) %>%
  dplyr::select(cancer_types, Exp_site,mean)%>%
  unique() %>%
  tidyr::spread(key="Exp_site",value="mean") %>%
  dplyr::mutate(`log2FC(I/T)`=log2(`Immune cell dominate`/`Tumor cell dominate`)) %>%
  dplyr::mutate(title = paste(cancer_types,", log2FC=",signif(`log2FC(I/T)`,2),sep="")) %>%
  dplyr::select(cancer_types,title,`log2FC(I/T)`) -> label
plot_ready %>%
  dplyr::inner_join(label,by="cancer_types") -> plot_ready

plot_ready %>%
  dplyr::select(Exp_site) %>%
  dplyr::inner_join(data.frame(Exp_site =c("Immune and tumor cell almost" ,"Immune cell dominate","Tumor cell dominate" ),
                               rank = c(3,2,1)), by = "Exp_site") %>%
  dplyr::arrange(rank) %>%
  .$Exp_site -> Exp_site.rank
plot_ready <- within(plot_ready,Exp_site <- factor(Exp_site,levels = unique(Exp_site.rank)))
with(plot_ready, levels(Exp_site))

plot_ready %>%
  dplyr::arrange(`log2FC(I/T)`) %>%
  .$title -> title.rank
plot_ready <- within(plot_ready,title <- factor(title,levels = unique(title.rank)))
with(plot_ready, levels(title))

plot_ready %>%
  dplyr::mutate(average_exp = log2(average_exp)) %>%
  ggplot(aes(x=Exp_site,y=average_exp)) +
  geom_violin(aes(fill=cancer_types),alpha=0.5, size = 0.5, color = "grey") +
  facet_wrap(~title) +
  coord_flip() +
  theme_bw() +
  xlab("Expression pattern of ICPs") +
  ylab("log2 (Average expression)") +
  ggpubr::stat_compare_means(comparisons = list(c("Tumor cell dominate","Immune cell dominate")),label = "p.signif") +
  theme(
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(size = 10,color = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black"),
    panel.grid = element_blank()
  )
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_expsite-by-cancers-noMHC.pdf"),device = "pdf", width = 12,height = 7)  
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_expsite-by-cancers-noMHC.png"),device = "png",width = 12,height = 7)  

# compare exp pattern among, with normal and tumor ---------
# expression in 3 groups by function: exp site Mainly_Tumor, Mainly_Immune, both -----
# 样本中在各基因的平均表达量

fn_average.by_expsite.TN <- function(.data){
  .data %>%
    dplyr::filter(!symbol %in% MHC_class_genes) %>%
    tidyr::gather(-symbol,-entrez_id,key="barcode",value="expr") %>%
    dplyr::left_join(gene_list,by="symbol") %>%
    dplyr::filter(!is.na(expr)) %>%
    dplyr::mutate(group = ifelse(substr(barcode,14,14)==0, "Tumor", "Normal")) %>%
    # dplyr::filter(Exp_site %in% c("Immune and tumor cell almost" ,"Immune cell dominate","Tumor cell dominate" )) %>%
    dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Only_exp_on_Immune","Mainly_exp_on_Immune"),"Mainly_exp_on_Immune","Mainly_exp_on_Tumor")) %>%
    dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Both_exp_on_Tumor_Immune"),"Both_exp_on_Tumor_Immune",Exp_site.1)) %>%
    dplyr::mutate(Exp_site = Exp_site.1) %>%
    dplyr::group_by(barcode,Exp_site, group) %>%
    dplyr::mutate(average_exp = mean(expr)) %>%
    dplyr::select(barcode,Exp_site,group,average_exp) %>%
    unique() %>%
    dplyr::ungroup()
}

gene_list_expr %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(ICP_mean_expr = purrr::map(filter_expr,fn_average.by_expsite.TN)) %>%
  dplyr::select(-filter_expr) %>%
  dplyr::ungroup() -> ICP_mean_expr_in_cancers.by_expsite.TN

ICP_mean_expr_in_cancers.by_expsite.TN %>%
  readr::write_rds(file.path(out_path,"e_6_exp_profile","ICP_mean_expr_in_TN.by_expsite-noMHC.rds.gz"),compress = "gz") 

# ICP_mean_expr_in_cancers.by_expsite.TN <- readr::read_rds(file.path(out_path,"e_6_exp_profile","ICP_mean_expr_in_TN.by_expsite-noMHC.rds.gz"))

ICP_mean_expr_in_cancers.by_expsite.TN %>%
  tidyr::unnest() %>%
  dplyr::mutate(average_exp = log2(average_exp)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(sum = quantile(average_exp,0.5)) %>%
  dplyr::select(cancer_types,sum) %>%
  unique() %>%
  dplyr::arrange(sum) -> cancer_rank.by_expsite.TN

# ICP_mean_expr_in_cancers.by_expsite %>%
#   tidyr::unnest() %>%
#   dplyr::filter(Exp_site!="N") %>%
#   dplyr::mutate(average_exp = log2(average_exp)) %>%
#   ggplot(aes(x=cancer_types,y=average_exp)) +
#   geom_violin(aes(fill=Exp_site),alpha=0.5) +
#   facet_wrap(~Exp_site) +
#   coord_flip() +
#   scale_x_discrete(limits= cancer_rank$cancer_types) +
#   scale_fill_manual(values = c("#9A32CD", "Blue" ,"red","pink","orange","grey")) +
#   theme_bw() +
#   xlab("Cancer Types") +
#   ylab("log2 (Expression)") +
#   theme(
#     strip.background = element_rect(colour = "black", fill = "white"),
#     strip.text = element_text(size = 12,color = "black"),
#     axis.text = element_text(size = 10, colour = "black"),
#     legend.position = "none",
#     panel.background = element_blank(),
#     panel.border = element_rect(fill='transparent',colour = "black"),
#     panel.grid = element_line(linetype = "dashed")
#   )
# ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_cancers-by_expsite.pdf"),device = "pdf", width = 6,height = 6)  
# ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_cancers-by_expsite.png"),device = "png",width = 6,height = 6)  

# class by cancers
ICP_mean_expr_in_cancers.by_expsite.TN %>%
  tidyr::unnest() %>%
  dplyr::filter(!Exp_site %in% c("N","Not_sure")) -> plot_ready

plot_ready %>%
  dplyr::group_by(cancer_types, Exp_site, group) %>%
  dplyr::mutate(mean = mean(average_exp)) %>%
  dplyr::select(cancer_types, Exp_site,mean)%>%
  unique() %>%
  tidyr::spread(key="Exp_site",value="mean") %>%
  # dplyr::mutate(`log2FC(I/T)`=log2(`Immune cell dominate`/`Tumor cell dominate`)) %>%
  dplyr::mutate(`log2FC(I/T)`=log2(Mainly_exp_on_Immune/Mainly_exp_on_Tumor)) %>%
  dplyr::mutate(title = paste(cancer_types,", log2FC=",signif(`log2FC(I/T)`,2),sep="")) %>%
  dplyr::select(cancer_types,title,`log2FC(I/T)`) -> label
plot_ready %>%
  dplyr::inner_join(label,by=c("cancer_types","group")) -> plot_ready

plot_ready %>%
  dplyr::select(Exp_site) %>%
  dplyr::inner_join(data.frame(Exp_site =c("Both_exp_on_Tumor_Immune" ,"Mainly_exp_on_Immune","Mainly_exp_on_Tumor" ), rank = c(3,2,1), Exp_site.s = c("TIC-ICGs","IC-ICGs","TC-ICGs")), by = "Exp_site") %>%
  dplyr::arrange(rank) %>%
  dplyr::select(-rank)%>%
  unique() -> Exp_site.rank
plot_ready %>%
  dplyr::inner_join(Exp_site.rank, by="Exp_site") %>%
  dplyr::mutate(Exp_site.s = paste(Exp_site.s,group,sep = ".")) -> plot_ready
# plot_ready <- within(plot_ready,Exp_site.s <- factor(Exp_site.s,levels = unique(Exp_site.rank$Exp_site.s)))
# with(plot_ready, levels(Exp_site.s))

plot_ready %>%
  dplyr::filter(group == "Tumor") %>%
  dplyr::arrange(`log2FC(I/T)`) %>%
  .$cancer_types -> c.rank
plot_ready <- within(plot_ready,cancer_types <- factor(cancer_types,levels = unique(c.rank)))
with(plot_ready, levels(cancer_types))


plot_ready %>%
  dplyr::mutate(average_exp = log2(average_exp)) %>%
  ggplot(aes(x=Exp_site.s,y=average_exp)) +
  geom_violin(aes(fill=group), size = 0.5, color = "black") +
  facet_wrap(~cancer_types) +
  coord_flip() +
  theme_bw() +
  xlab("Expression pattern of ICGs") +
  ylab("Average expression (log2)") +
  scale_fill_manual(values=c("#1E90FF","#CD9B1D")) +
  ggpubr::stat_compare_means(comparisons = list(c("IC-ICGs.Normal","IC-ICGs.Tumor"),
                                                c("TC-ICGs.Normal","TC-ICGs.Tumor"),
                                                c("TIC-ICGs.Normal","TIC-ICGs.Tumor"),
                                                c("IC-ICGs.Tumor","TC-ICGs.Tumor"),
                                                c("IC-ICGs.Normal","TC-ICGs.Normal")),
                             label = "p.signif") +
  theme(
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(size = 10,color = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    # legend.position = "none",
    legend.title = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black"),
    panel.grid = element_blank()
  )
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_expsite-by-TN-noMHC.pdf"),device = "pdf", width = 12,height = 10)

## added by chunjie
# split figure into 2 figure
# 1 compare the tumor between ICG
# 2 compare the ICG between normal and tumor

gene_color <- c('IC-ICGs' = '#bbf4b0', 'TC-ICGs' = '#f4c2c2', 'TIC-ICGs' = '#e7ddac')
plot_ready_tumor <- plot_ready %>% 
  dplyr::mutate(cancer_types = as.character(cancer_types)) %>% 
  dplyr::filter(group == 'Tumor') %>% 
  dplyr::mutate(average_exp = log2(average_exp)) %>% 
  dplyr::mutate(average_exp = ifelse(average_exp > 12, 12, average_exp)) %>% 
  dplyr::mutate(average_exp = ifelse(average_exp < 5.5, 5.5, average_exp))

plot_ready_tumor %>% 
  dplyr::filter(Exp_site.s == 'IC-ICGs.Tumor') %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(m = mean(average_exp)) %>% 
  dplyr::arrange(-m) %>% 
  dplyr::pull(cancer_types) ->
  plot_ready_tumor_tumor_rank

plot_ready_tumor %>% 
  ggplot(aes(x = cancer_types, y = average_exp, color = Exp_site.s)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  scale_color_manual(values = unname(gene_color), name = 'ICG Type', label = c('IC-ICG', 'TC-ICG', 'TIC-ICG')) +
  scale_x_discrete(limit = plot_ready_tumor_tumor_rank) +
  labs(x = 'Cancer Types', y = 'Expression(log2)') +
  theme(
    plot.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.line.x.bottom = element_line(color = 'black'),
    axis.line.y.left = element_line(color = 'black'),
    legend.position = 'top',
    legend.background = element_blank(),
    legend.key = element_blank()
  ) -> plot_within_tumor

ggsave(filename = 'plot_within_tumor.pdf', plot = plot_within_tumor, device = 'pdf', path = '/project/huff/huff/bib-newoutput', height = 3, width = 11)

# plot2
# 2 compare the ICG between normal and tumor
plot_ready %>% 
  dplyr::mutate(cancer_types = as.character(cancer_types)) %>% 
  dplyr::group_by(cancer_types) %>% 
  tidyr::nest() %>% 
  dplyr::mutate(data = purrr::map(.x = data, .f = function(.x) {
    .x %>% 
      dplyr::select(barcode) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(patient = substr(x = barcode, start = 1, stop = 12)) %>% 
      dplyr::mutate(code = substr(x = barcode, start = 14, stop = 15)) %>% 
      dplyr::filter(code %in% c('01', '11')) %>% 
      dplyr::group_by(patient) %>% 
      dplyr::mutate(n = dplyr::n()) %>% 
      dplyr::filter(n > 1) %>% dplyr::ungroup() %>% dplyr::pull(barcode) -> 
      .pair_barcode
    .x %>% dplyr::filter(barcode %in% .pair_barcode)
  })) %>% 
  dplyr::filter(purrr::map_lgl(.x = data, .f = function(.x){nrow(.x)> 60})) %>% 
  tidyr::unnest(cols = data) %>% 
  dplyr::ungroup() ->
  plot_ready_tumor_normal

plot_ready_tumor_normal %>% 
  dplyr::mutate(average_exp = log2(average_exp)) %>% 
  dplyr::mutate(average_exp = ifelse(average_exp > 11, 11, average_exp)) %>% 
  dplyr::mutate(average_exp = ifelse(average_exp < 6, 6, average_exp)) %>% 
  dplyr::mutate(icg_type = gsub(pattern = 's.Normal|s.Tumor', replacement = '', x = Exp_site.s)) %>% 
  dplyr::mutate(group = factor(x = group, levels = c('Tumor', 'Normal'))) %>% 
  ggplot(aes(x = cancer_types, y = average_exp, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  scale_color_manual(values = c('#bfa500', '#6775ff'), name = '') +
  facet_wrap(~icg_type, ncol = 1, strip.position = 'right') +
  labs(x = 'Cancer Types', y = 'Expression(log2)') +
  theme(
    plot.title = element_blank(),
    panel.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.line.x.bottom = element_line(color = 'black'),
    axis.line.y.left = element_line(color = 'black'),
    legend.position = 'top',
    legend.background = element_blank(),
    legend.key = element_blank(),
    strip.background = element_blank()
  ) -> plot_compare_tumor_normal
ggsave(filename = 'plot_compare_tumor_normal.pdf', plot = plot_compare_tumor_normal, device = 'pdf', path = '/project/huff/huff/bib-newoutput', height = 4, width = 11)

load(file = '/project/huff/huff/bib-newoutput/cj_02_expr_proile.rda')

save.image(file = '/project/huff/huff/bib-newoutput/cj_02_expr_proile.rda')


# survival analysis [univariable survival analysis] ----------------------------
## PFS
ICP_mean_expr_in_cancers.by_expsite %>%
  tidyr::unnest() %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical,by="barcode") %>%
  dplyr::rename("time" = "PFS.time", "status" = "PFS") %>%
  dplyr::group_by(cancer_types,Exp_site) %>%
  dplyr::mutate(group = ifelse(average_exp >= quantile(average_exp,0.5),"2High","1Low")) %>%
  dplyr::ungroup() %>%
  tidyr::nest(-cancer_types,-Exp_site) %>%
  dplyr::mutate(survival = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ICP_mean_expr_in_cancers.byexpsite.PFS
ICP_mean_expr_in_cancers.byexpsite.PFS  %>%
  readr::write_tsv(file.path(out_path,"e_6_exp_profile","survival_byExpsite","PFS_survival.ICP_mean_expr_in_cancers.byexpsite-noMHC.tsv"))

## OS
ICP_mean_expr_in_cancers.by_expsite %>%
  tidyr::unnest() %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::rename("time" = "OS", "status" = "Status") %>%
  dplyr::group_by(cancer_types,Exp_site) %>%
  dplyr::mutate(group = ifelse(average_exp >= quantile(average_exp,0.5),"2High","1Low")) %>%
  dplyr::ungroup() %>%
  tidyr::nest(-cancer_types,-Exp_site) %>%
  dplyr::mutate(survival = purrr::map2(data,cancer_types,fun_clinical_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ICP_mean_expr_in_cancers.byexpsite.OS

ICP_mean_expr_in_cancers.byexpsite.OS  %>%
  readr::write_tsv(file.path(out_path,"e_6_exp_profile","survival_byExpsite","OS_survival.ICP_mean_expr_in_cancers.byexpsite-noMHC.tsv"))

## survival plot, cox [univariable]
# PFS, cox, group
ICP_mean_expr_in_cancers.byexpsite.PFS %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  tidyr::nest(-Exp_site) %>%
  dplyr::mutate(data = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort))
  })) %>%
  dplyr::mutate(filename = paste("Meanexp.COX_PFS.by-Expsite-noMHC",Exp_site,sep=".")) %>%
  dplyr::rename("functionWithImmune"="Exp_site") %>%
  purrr::pwalk(.f=fn_cox_plot,hr="hr",hr_l="hr_l",hr_h="hr_h",title="Progression-free survival",facet="~ functionWithImmune",dir = "survival_byExpsite",w = 4, h = 6)

ICP_mean_expr_in_cancers.byexpsite.PFS %>% 
  dplyr::filter(Exp_site != "Immune and tumor cell almost") %>%
  dplyr::rename("functionWithImmune"="Exp_site") %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(hr_sum = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::select(hr,functionWithImmune) %>%
      tidyr::spread(key="functionWithImmune",value="hr") %>%
      dplyr::mutate(hr_sum=`Tumor cell dominate`-`Immune cell dominate`) %>%
      dplyr::select(hr_sum)-> tmp
    rbind(tmp,tmp)
  })) %>%
  tidyr::unnest() %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  fn_cox_plot.all(filename="Meanexp.COX_PFS.by-Expsite.all-noMHC-reorder",hr="hr",hr_l="hr_l",hr_h="hr_h",title="Progression-free survival",facet="~ functionWithImmune",dir = "survival_byExpsite",w = 6, h = 6)

# PFS, cox, continus
ICP_mean_expr_in_cancers.byexpsite.PFS %>% 
  dplyr::rename("functionWithImmune"="Exp_site") %>%
  dplyr::mutate(cox_sig = ifelse(hr.c.coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr.c=log2(hr.c)+1,hr.c.l=log2(hr.c.l)+1,hr.c.h=log2(hr.c.h)+1) %>%
  tidyr::nest(-functionWithImmune) %>%
  dplyr::mutate(data = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::mutate(cancer_types = reorder(cancer_types,hr.c,sort))
  })) %>%
  dplyr::mutate(filename = paste("Meanexp.COX_PFS.by-Expsite.continus",functionWithImmune,sep=".")) %>%
  purrr::pwalk(.f=fn_cox_plot,hr="hr.c",hr_l="hr.c.l",hr_h="hr.c.h",title="Progression-free survival",facet="~ functionWithImmune",dir = "survival_byExpsite",w = 4, h = 6)

ICP_mean_expr_in_cancers.byexpsite.PFS %>% 
  dplyr::rename("functionWithImmune"="Exp_site") %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr.c=log2(hr.c)+1,hr.c.l=log2(hr.c.l)+1,hr.c.h=log2(hr.c.h)+1) %>%
  fn_cox_plot.all(filename="Meanexp.COX_PFS.by-expsite.all.continus",hr="hr.c",hr_l="hr.c.l",hr_h="hr.c.h",title="Progression-free survival",facet="~ functionWithImmune",dir = "survival_byExpsite",w = 8, h = 6)

# OS, cox, group
ICP_mean_expr_in_cancers.byexpsite.OS %>% 
  dplyr::rename("functionWithImmune"="Exp_site") %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  tidyr::nest(-functionWithImmune) %>%
  dplyr::mutate(data = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort))
  })) %>%
  dplyr::mutate(filename = paste("Meanexp.COX_OS.by-Expsite",functionWithImmune,sep=".")) %>%
  purrr::pwalk(.f=fn_cox_plot,hr="hr",hr_l="hr_l",hr_h="hr_h",title="Overall survival",facet="~ functionWithImmune",dir = "survival_byExpsite",w = 4, h = 6)

ICP_mean_expr_in_cancers.byexpsite.OS %>% 
  dplyr::filter(Exp_site != "Immune and tumor cell almost") %>%
  dplyr::rename("functionWithImmune"="Exp_site") %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(hr_sum = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::select(hr,functionWithImmune) %>%
      tidyr::spread(key="functionWithImmune",value="hr") %>%
      dplyr::mutate(hr_sum=`Immune cell dominate`-`Tumor cell dominate`) %>%
      dplyr::select(hr_sum)-> tmp
    rbind(tmp,tmp)
  })) %>%
  tidyr::unnest() %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  fn_cox_plot.all(filename="Meanexp.COX_OS.by-Expsite.all-noMHC",hr="hr",hr_l="hr_l",hr_h="hr_h",title="Overall survival",facet="~ functionWithImmune",dir = "survival_byExpsite",w = 6, h = 6)

# OS, cox, continus
ICP_mean_expr_in_cancers.byexpsite.OS %>% 
  dplyr::rename("functionWithImmune"="Exp_site") %>%
  dplyr::mutate(cox_sig = ifelse(hr.c.coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr.c=log2(hr.c)+1,hr.c.l=log2(hr.c.l)+1,hr.c.h=log2(hr.c.h)+1) %>%
  tidyr::nest(-functionWithImmune) %>%
  dplyr::mutate(data = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::mutate(cancer_types = reorder(cancer_types,hr.c,sort))
  })) %>%
  dplyr::mutate(filename = paste("Meanexp.COX_OS.by-Expsite.continus",functionWithImmune,sep=".")) %>%
  purrr::pwalk(.f=fn_cox_plot,hr="hr.c",hr_l="hr.c.l",hr_h="hr.c.h",title="Overall survival",facet="~ functionWithImmune",dir = "survival_byExpsite",w = 4, h = 6)

ICP_mean_expr_in_cancers.byexpsite.OS %>% 
  dplyr::rename("functionWithImmune"="Exp_site") %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr.c=log2(hr.c)+1,hr.c.l=log2(hr.c.l)+1,hr.c.h=log2(hr.c.h)+1) %>%
  fn_cox_plot.all(filename="Meanexp.COX_OS.by-expsite.all.continus",hr="hr.c",hr_l="hr.c.l",hr_h="hr.c.h",title="Overall survival",facet="~ functionWithImmune",dir = "survival_byExpsite",w = 8, h = 6)

#### survival analysis [multi-variable survival analysis] ----------------------------

# calculation of average expression of samples of genes ----
# 基因在各样本中的平均表达量
fn_average_a <- function(.data){
  .data %>%
    tidyr::gather(-symbol,-entrez_id,key="barcode",value="expr") %>%
    dplyr::filter(!is.na(expr)) %>%
    dplyr::filter(substr(barcode,14,14)==0) %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(average_exp = mean(expr)) %>%
    dplyr::select(symbol,average_exp) %>%
    unique() %>%
    dplyr::ungroup()
}

gene_list_expr %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(ICP_mean_expr = purrr::map(filter_expr,fn_average_a)) %>%
  dplyr::select(-filter_expr) %>%
  dplyr::ungroup() -> mean_exp_of_ICP_in_cancers_samples

# heatmap -----------------
mean_exp_of_ICP_in_cancers_samples %>%
  tidyr::unnest() %>%
  tidyr::spread(key="symbol",value="average_exp") %>%
  as.data.frame() -> mean_exp_of_ICP_in_cancers_samples.df

mean_exp_of_ICP_in_cancers_samples %>%
  tidyr::unnest() %>%
  dplyr::select(symbol) %>%
  dplyr::inner_join(gene_list,by="symbol") %>%
  dplyr::arrange(symbol) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("N","Not_sure"),"Not_sure",Exp_site)) %>%
  dplyr::mutate(Exp_site = Exp_site.1) %>%
  dplyr::mutate(Exp_site = ifelse(Exp_site=="N","Not_sure",Exp_site)) %>%
  dplyr::select(symbol,type,functionWithImmune,Exp_site,family) %>%
  unique() %>%
  dplyr::rename("Function type"="type","Immunity"="functionWithImmune","Exp. pattern" = "Exp_site","Gene family" = "family") %>%
  as.data.frame() -> symbol_anno
rownames(symbol_anno) <- symbol_anno[,1]
symbol_anno <- symbol_anno[,-1]
library(ComplexHeatmap)
# gene row annotation
gene_anno <- HeatmapAnnotation(df=symbol_anno,
                               col = list(Immunity=c("Inhibit" = "#FF82AB",
                                                     "Activate" = "#1E90FF",
                                                     "TwoSide" = "#FFB90F"),
                                          `Exp. pattern`=c("Tumor cell dominate" = "#FFC1C1",
                                                     "Immune and tumor cell almost" = "#EED8AE",
                                                     "Immune cell dominate" = "#B4EEB4",
                                                     "Not_sure" = "grey"),
                                          `Function type` = c("Receptor" = "#9400D3",
                                                   "Ligand" = "#556B2F",
                                                   "Ligand&Receptor" = "#76EE00"),
                                          `Gene family` = c("BTN" = "#838B8B", 
                                                     "KIR_activate"= "#000000",
                                                     "KIR_inhibit"  ="#104E8B",
                                                     "MHC class I"= "#8B2323",
                                                     "MHC class II"="#CDAA7D",
                                                     "Other"="#8EE5EE")),
                               height = unit(0.2, "cm"),
                               show_annotation_name =T)

draw(gene_anno,1:69)

rownames(mean_exp_of_ICP_in_cancers_samples.df) <- mean_exp_of_ICP_in_cancers_samples.df[,1]
mean_exp_of_ICP_in_cancers_samples.df <- mean_exp_of_ICP_in_cancers_samples.df[,-1]

mean_exp_of_ICP_in_cancers_samples.df <- as.matrix(mean_exp_of_ICP_in_cancers_samples.df) # complexheatmap need matrix data form

mean_exp_of_ICP_in_cancers_samples.df.scaled <- apply(mean_exp_of_ICP_in_cancers_samples.df,1,scale) %>% t()
colnames(mean_exp_of_ICP_in_cancers_samples.df.scaled) <- colnames(mean_exp_of_ICP_in_cancers_samples.df)

library(circlize)
library(dendextend)

col_dend = hclust(dist(t(mean_exp_of_ICP_in_cancers_samples.df.scaled))) # column clustering
row_dend = hclust(dist((mean_exp_of_ICP_in_cancers_samples.df.scaled))) # row clustering 
pdf(file.path(out_path,"e_6_exp_profile","ICP_cluster_tree.pdf"),width = 12,height = 4)
plot(col_dend)
dev.off()

plot(row_dend)
he = Heatmap(mean_exp_of_ICP_in_cancers_samples.df.scaled,
             col = colorRamp2(c(-1, 0, 5), c("blue", "white", "red")),
             row_names_gp = gpar(fontsize = 8),
             show_row_names = T, 
             show_column_names = FALSE,
             show_row_dend = F, # whether show row clusters.
             top_annotation = gene_anno,
             row_names_side = c("left"),
             cluster_columns = color_branches(col_dend, k = 6),   # add color on the column tree branches
             # cluster_rows = color_branches(row_dend, k = 4), 
             heatmap_legend_param = list(title = c("Mean exp.")))
pdf(file.path(out_path,"e_6_exp_profile","ICP_exp_in_cancers.pdf"),width = 8,height = 6)
he
dev.off()

tiff(file.path(out_path,"e_6_exp_profile","ICP_exp_in_cancers.tiff"),width = 800,height = 500)
he
dev.off()

# GSVA score in 3 groups by function: exp site Mainly_Tumor, Mainly_Immune, both -----
# GSVA score of exp sites
tcga_gsvascore <- readr::read_rds(file.path(immune_res_path,"TCGA_GSVAScore","TCGA_cancer_specific.allsamples(T-N)_GSVA.score_ICPs_features.rds.gz")) 

tcga_gsvascore %>%
  tidyr::unnest() %>%
  dplyr::select(cancer_types,barcode,Both_exp_on_Tumor_Immune,Mainly_exp_on_Immune,Mainly_exp_on_Tumor) %>%
  tidyr::gather(-cancer_types,-barcode,key="Exp_site",value="GSVA score") -> plot_ready.gsva
plot_ready.gsva %>%
  dplyr::group_by(cancer_types, Exp_site) %>%
  dplyr::mutate(mean = mean(`GSVA score`)) %>%
  dplyr::select(cancer_types, Exp_site,mean)%>%
  unique() %>%
  tidyr::spread(key="Exp_site",value="mean") %>%
  dplyr::mutate(diff_IT=Mainly_exp_on_Immune-Mainly_exp_on_Tumor) %>%
  dplyr::select(cancer_types,diff_IT) -> label
plot_ready.gsva %>%
  dplyr::inner_join(label,by="cancer_types") -> plot_ready.gsva

plot_ready.gsva %>%
  dplyr::select(Exp_site) %>%
  dplyr::inner_join(data.frame(Exp_site = c("Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor"),
                               rank = c(3,2,1)), by = "Exp_site") %>%
  dplyr::arrange(rank) %>%
  .$Exp_site -> Exp_site.rank
plot_ready.gsva <- within(plot_ready.gsva,Exp_site <- factor(Exp_site,levels = unique(Exp_site.rank)))
with(plot_ready.gsva, levels(Exp_site))

plot_ready.gsva %>%
  dplyr::arrange(diff_IT) %>%
  .$cancer_types -> cancer.rank
plot_ready.gsva <- within(plot_ready.gsva,cancer_types <- factor(cancer_types,levels = unique(cancer.rank)))
with(plot_ready.gsva, levels(cancer_types))

plot_ready.gsva %>%
  ggplot(aes(x=Exp_site,y=`GSVA score`)) +
  geom_violin(aes(fill=cancer_types),alpha=0.5) +
  facet_wrap(~cancer_types) +
  coord_flip() +
  theme_bw() +
  xlab("Expression pattern of ICPs") +
  ylab("GSVA score") +
  ggpubr::stat_compare_means(method = "anova",label = "p.signif") +
  theme(
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(size = 10,color = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black"),
    panel.grid = element_line(linetype = "dashed")
  )
ggsave(file.path(out_path,"e_6_exp_profile","ICP_GSVA-score_in_expsite-by-cancers.pdf"),device = "pdf", width = 12,height = 7)  
ggsave(file.path(out_path,"e_6_exp_profile","ICP_GSVA-score_in_expsite-by-cancers.png"),device = "png",width = 12,height = 7)
# save image --------------------------------------------------------------

save.image(file.path(out_path,"e_6_exp_profile","e_6_exp_profile.rdata"))

