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

# load data
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol <- as.character(gene_list$symbol)
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))
# ICP_expr_pattern <- readr::read_tsv(file.path(out_path,"ICP_exp_patthern","manual_edit_2_ICP_exp_pattern_in_immune_tumor_cell.tsv"))
ICP_expr_pattern <- readr::read_tsv(file.path(out_path,"ICP_exp_patthern-byratio/pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(entrez_ID,symbol,Exp_site,`log2FC(I/T)`) 

fn_site_color <- function(.x){
  # print(.n)
  if(.x=="Mainly_exp_on_Tumor"){
    "red"
  }else if(.x=="Mainly_exp_on_Immune"){
    "Blue"
  }else if(.x=="Only_exp_on_Immune"){
    "Blue"
  }else if(.x=="Both_exp_on_Tumor_Immune"){
    c("#9A32CD")
  }else{
    "grey"
  }
}
gene_list %>%
  dplyr::left_join(ICP_expr_pattern,by="symbol") %>%
  dplyr::mutate(Exp_site=ifelse(is.na(Exp_site),"N",Exp_site)) %>%
  dplyr::mutate(site_col = purrr::map(Exp_site,fn_site_color)) -> gene_list

survival_path <- "/home/huff/project/data/TCGA-survival-time/cell.2018.survival"
survival_data <- readr::read_rds(file.path(survival_path, "TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>%
  dplyr::rename("cancer_types" = "type")

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

# calculation of average expression of genes of samples --------------------------------
# 3 groups by function: activate, inhibit and two sides -----
# 样本中在各基因的平均表达量

fn_average <- function(.data){
  .data %>%
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
  readr::write_rds(file.path(out_path,"e_6_exp_profile","ICP_mean_expr_in_cancers.by_functionalRole.rds.gz"),compress = "gz")

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
  geom_violin(aes(fill=cancer_types),alpha=0.5) +
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
    panel.grid = element_line(linetype = "dashed")
  )
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_functionRoles-bycancers.pdf"),device = "pdf", width = 6,height = 5)
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_functionRoles-bycancer.png"),device = "png",width = 6,height =5)

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
  readr::write_tsv(file.path(out_path,"e_6_exp_profile","survival","PFS_survival.ICP_mean_expr_in_cancers.byfunction.tsv"))
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
  readr::write_tsv(file.path(out_path,"e_6_exp_profile","survival","OS_survival.ICP_mean_expr_in_cancers.byfunction.tsv"))

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

# PFS, cox, group
ICP_mean_expr_in_cancers.byfunction.PFS %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  tidyr::nest(-functionWithImmune) %>%
  dplyr::mutate(data = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort))
  })) %>%
  dplyr::mutate(filename = paste("Meanexp.COX_PFS.by-functionRole",functionWithImmune,sep=".")) %>%
  purrr::pwalk(.f=fn_cox_plot,hr="hr",hr_l="hr_l",hr_h="hr_h",title="Progression-free survival",facet="~ functionWithImmune",dir = "survival",w = 4, h = 6)

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

# 3 groups by function: exp site Mainly_Tumor, Mainly_Immune, both -----
# 样本中在各基因的平均表达量

fn_average.by_expsite <- function(.data){
  .data %>%
    tidyr::gather(-symbol,-entrez_id,key="barcode",value="expr") %>%
    dplyr::left_join(gene_list,by="symbol") %>%
    dplyr::filter(!is.na(expr)) %>%
    dplyr::filter(substr(barcode,14,14)==0) %>%
    dplyr::filter(Exp_site %in% c("Only_exp_on_Immune","Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor","Only_exp_on_Tumor")) %>%
    dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Only_exp_on_Immune","Mainly_exp_on_Immune"),"Mainly_exp_on_Immune","Mainly_exp_on_Tumor")) %>%
    dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Both_exp_on_Tumor_Immune"),"Both_exp_on_Tumor_Immune",Exp_site.1)) %>%
    dplyr::mutate(Exp_site = Exp_site.1) %>%
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
  readr::write_rds(file.path(out_path,"e_6_exp_profile","ICP_mean_expr_in_cancers.by_expsite.rds.gz"),compress = "gz") 

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
  dplyr::mutate(`log2FC(I/T)`=log2(Mainly_exp_on_Immune/Mainly_exp_on_Tumor)) %>%
  dplyr::mutate(title = paste(cancer_types,", log2FC=",signif(`log2FC(I/T)`,2),sep="")) %>%
  dplyr::select(cancer_types,title,`log2FC(I/T)`) -> label
plot_ready %>%
  dplyr::inner_join(label,by="cancer_types") -> plot_ready

plot_ready %>%
  dplyr::select(Exp_site) %>%
  dplyr::inner_join(data.frame(Exp_site = c("Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor"),
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
  geom_violin(aes(fill=cancer_types),alpha=0.5) +
  facet_wrap(~title) +
  coord_flip() +
  theme_bw() +
  xlab("Expression site of ICPs") +
  ylab("log2 (Expression)") +
  ggpubr::stat_compare_means(comparisons = list(c("Mainly_exp_on_Tumor","Mainly_exp_on_Immune")),label = "p.signif") +
  theme(
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(size = 8,color = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black"),
    panel.grid = element_line(linetype = "dashed")
  )
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_expsite-by-cancers.pdf"),device = "pdf", width = 8,height = 5)  
ggsave(file.path(out_path,"e_6_exp_profile","ICP_average_exp_in_expsite-by-cancers.png"),device = "png",width = 8,height = 5)  

# survival analysis
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
  readr::write_tsv(file.path(out_path,"e_6_exp_profile","survival_byExpsite","PFS_survival.ICP_mean_expr_in_cancers.byexpsite.tsv"))

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
  readr::write_tsv(file.path(out_path,"e_6_exp_profile","survival_byExpsite","OS_survival.ICP_mean_expr_in_cancers.byexpsite.tsv"))
## survival plot, cox

# PFS, cox, group
ICP_mean_expr_in_cancers.byexpsite.PFS %>% 
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  tidyr::nest(-Exp_site) %>%
  dplyr::mutate(data = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort))
  })) %>%
  dplyr::mutate(filename = paste("Meanexp.COX_PFS.by-Expsite",Exp_site,sep=".")) %>%
  dplyr::rename("functionWithImmune"="Exp_site") %>%
  purrr::pwalk(.f=fn_cox_plot,hr="hr",hr_l="hr_l",hr_h="hr_h",title="Progression-free survival",facet="~ functionWithImmune",dir = "survival_byExpsite",w = 4, h = 6)

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

# heatmap
mean_exp_of_ICP_in_cancers_samples %>%
  tidyr::unnest() %>%
  tidyr::spread(key="symbol",value="average_exp") %>%
  as.data.frame() -> mean_exp_of_ICP_in_cancers_samples.df

mean_exp_of_ICP_in_cancers_samples %>%
  tidyr::unnest() %>%
  dplyr::select(symbol) %>%
  dplyr::inner_join(gene_list,by="symbol") %>%
  dplyr::arrange(symbol) %>%
  dplyr::select(symbol,type,functionWithImmune,Exp_site) %>%
  unique() %>%
  dplyr::mutate(Exp_site = ifelse(Exp_site=="N","Not_sure",Exp_site)) %>%
  as.data.frame() -> symbol_anno
rownames(symbol_anno) <- symbol_anno[,1]
symbol_anno <- symbol_anno[,-1]
library(ComplexHeatmap)
# gene row annotation
gene_anno <- HeatmapAnnotation(df=symbol_anno,
                               col = list(functionWithImmune=c("Inhibit" = "#EE3B3B",
                                                               "Activate" = "#1C86EE",
                                                               "TwoSide" = "#EE7600"),
                                          Exp_site=c("Mainly_exp_on_Tumor" = "pink",
                                                     "Only_exp_on_Tumor" = "red",
                                                     "Both_exp_on_Tumor_Immune" = "#9A32CD",
                                                     "Mainly_exp_on_Immune" = "green",
                                                     "Only_exp_on_Immune" = "blue",
                                                     "Not_sure" = "grey"),
                                          type = c("Receptor" = "black",
                                                   "Ligand" = "red",
                                                   "Ligand&Receptor" = "purple")),
                               width = unit(0.2, "cm"),
                               name = c("Immunity","Type"))

draw(gene_anno,1:20)

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
             col = colorRamp2(c(-2, 0, 4), c("blue", "white", "red")),
             row_names_gp = gpar(fontsize = 8),
             show_row_names = T, 
             show_column_names = FALSE,
             show_row_dend = T, # whether show row clusters.
             top_annotation = gene_anno,
             row_names_side = c("left"),
             cluster_columns = color_branches(col_dend, k = 6),   # add color on the column tree branches
             cluster_rows = color_branches(row_dend, k = 4), 
             heatmap_legend_param = list(title = c("Scaled Exp.")))
pdf(file.path(out_path,"e_6_exp_profile","ICP_exp_in_cancers.pdf"),width = 6,height = 4)
he
dev.off()

tiff(file.path(out_path,"e_6_exp_profile","ICP_exp_in_cancers.tiff"),width = 6,height = 4)
he
dev.off()

# save image --------------------------------------------------------------

save.image(file.path(out_path,"e_6_exp_profile","e_6_exp_profile.rdata"))
