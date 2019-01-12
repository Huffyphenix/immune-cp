###########################################################################
# survival analysis for ICP mutation --------------------------------------
# ICP mutated - non_ICP mutated

library(magrittr)
library(ggplot2)

out_path <- "/project/huff/huff/immune_checkpoint/result_20171025"
snv_path <- file.path(out_path, "m_2_snv")
gene_list_snv.hypermutation_class <- readr::read_rds(file.path(snv_path,"gene_list_snv.hypermutation_class"))

gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)

clinical <- readr::read_rds(file.path("/project/huff/huff/data/TCGA-survival-time/cell.2018.survival","TCGA_pancan_cancer_cell_survival_time.rds.gz"))
cancers_except <- c("COADREAD","GBMLGG","KIPAN","STES")

clinical %>%
  dplyr::filter(! type %in% cancers_except) %>%
  # dplyr::mutate(cli_snv_merge = purrr::map2(.x=clinical,type,fn_merge)) %>%
  # dplyr::select(-clinical) %>%
  tidyr::unnest() %>%
  dplyr::select(bcr_patient_barcode,PFS,PFS.time) %>%
  tidyr::drop_na() %>%
  dplyr::rename("barcode" = "bcr_patient_barcode")-> time_status

# survival function
fn_survival <- function(cancer,.data){
  print(cancer)
  .data %>%
    dplyr::inner_join(gene_list,by="symbol") %>%
    dplyr::filter(! class == "hyper") -> .tmp
  .tmp %>%
    dplyr::filter(functionWithImmune == "Inhibit") %>%
    dplyr::select(barcode,symbol,count) %>%
    dplyr::group_by(barcode) %>% 
    dplyr::mutate(n = sum(count)) %>%
    dplyr::select(barcode,n) %>%
    unique() %>%
    dplyr::ungroup() -> inhibt_mut
  .tmp %>%
    dplyr::filter(functionWithImmune == "Activate") %>%
    dplyr::select(barcode,symbol,count) %>%
    dplyr::group_by(barcode) %>% 
    dplyr::mutate(n = sum(count)) %>%
    dplyr::select(barcode,n) %>%
    unique() %>%
    dplyr::ungroup() -> Activate_mut
  .tmp %>%
    dplyr::filter(functionWithImmune == "TwoSide") %>%
    dplyr::select(barcode,symbol,count) %>%
    dplyr::group_by(barcode) %>% 
    dplyr::mutate(n = sum(count)) %>%
    dplyr::select(barcode,n) %>%
    unique() %>%
    dplyr::ungroup() -> TwoSide_mut
  .tmp %>%
    dplyr::select(barcode,symbol,count) %>%
    dplyr::group_by(barcode) %>% 
    dplyr::mutate(n = sum(count)) %>%
    dplyr::select(barcode,n) %>%
    unique() %>%
    dplyr::ungroup() -> ICP_mut
  
  # ICP mutation survival----
  ICP_mut %>%
    dplyr::inner_join(time_status,by="barcode") %>%
    dplyr::rename("time" = "PFS.time", "status" = "PFS") %>%
    dplyr::mutate(group = ifelse(n>0,"ICP Mutated","Wild Type")) -> ICP_sur_data
  
  if (ICP_sur_data$group %>% table() %>% length() >1) {
    title <- paste(cancer,"ICPs-mutated")
    color_list <- c("red","blue")
    sur_name <- paste(cancer,"ICPs-mutated.survival.pdf",sep=".")
    ICP.res.kmp <- fn_draw_survival(ICP_sur_data,title,color_list,sur_name,file.path(snv_path,"snv_survival/ICP_class.survival"),3,4)
  }else{
    ICP.res.kmp <- 1
  }
  # Inhibited ICP mutation survival----
  inhibt_mut %>%
    dplyr::inner_join(time_status,by="barcode") %>%
    dplyr::rename("time" = "PFS.time", "status" = "PFS") %>%
    dplyr::mutate(group = ifelse(n>0,"I-ICP Mutated","Wild Type")) -> Inhibted_ICP_sur_data
  
  if (Inhibted_ICP_sur_data$group %>% table() %>% length() >1) {
    title <- paste(cancer,"Inhibted-ICPs-mutated")
    color_list <- c("red","blue")
    sur_name <- paste(cancer,"Inhibted-ICPs-mutated.survival.pdf",sep=".")
    Inhibted_ICP.res.kmp <- fn_draw_survival(Inhibted_ICP_sur_data,title,color_list,sur_name,file.path(snv_path,"snv_survival/ICP_class.survival"),3,4)
  }else{
    Inhibted_ICP.res.kmp <- 1
  }
  # Activate_ICP mutation survival----
  Activate_mut %>%
    dplyr::inner_join(time_status,by="barcode") %>%
    dplyr::rename("time" = "PFS.time", "status" = "PFS") %>%
    dplyr::mutate(group = ifelse(n>0,"A-ICP Mutated","Wild Type")) -> Activate_ICP_sur_data
  
  if (Activate_ICP_sur_data$group %>% table() %>% length() >1) {
    title <- paste(cancer,"Activated-ICPs-mutated")
    color_list <- c("red","blue")
    sur_name <- paste(cancer,"Activated-ICPs-mutated.survival.pdf",sep=".")
    Activated_ICP.res.kmp <- fn_draw_survival(Activate_ICP_sur_data,title,color_list,sur_name,file.path(snv_path,"snv_survival/ICP_class.survival"),3,4)
  }else{
    Activated_ICP.res.kmp <- 1
  }
  # ICP mutation survival----
  TwoSide_mut %>%
    dplyr::inner_join(time_status,by="barcode") %>%
    dplyr::rename("time" = "PFS.time", "status" = "PFS") %>%
    dplyr::mutate(group = ifelse(n>0,"T-ICPs Mutated","Wild Type")) -> TwoSide_ICP_sur_data
  
  if (TwoSide_ICP_sur_data$group %>% table() %>% length() >1) {
    title <- paste(cancer,"TwoSide-ICPs-mutated")
    color_list <- c("red","blue")
    sur_name <- paste(cancer,"TwoSide-ICPs-mutated.survival.pdf",sep=".")
    TwoSide_ICP.res.kmp <- fn_draw_survival(TwoSide_ICP_sur_data,title,color_list,sur_name,file.path(snv_path,"snv_survival/ICP_class.survival"),3,4)
  }else{
    TwoSide_ICP.res.kmp <- 1
  }
  tibble::tibble(group = c("ICP_mut","Inhibit_mut","Activate_mut","TwoSide_mut"),kmp = c(ICP.res.kmp,Inhibted_ICP.res.kmp,Activated_ICP.res.kmp,TwoSide_ICP.res.kmp))
}

fn_draw_survival <- function(data,title,color,sur_name,result_path,h,w){
  library(survival)
  library(survminer)
  fit <- survfit(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)
  diff <- survdiff(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)
  kmp <- 1 - pchisq(diff$chisq, df = length(levels(as.factor(data$group))) - 1)
  legend <- data.frame(group=paste(sort(unique(data$group)),sep=""),n=fit$n)
  legend %>%
    dplyr::mutate(
      label = purrr::map2(
        .x = group,
        .y = n,
        .f = function(.x,.y){
          latex2exp::TeX(glue::glue("<<.x>>, n = <<.y>>", .open = "<<", .close = ">>"))
        }
      )
    ) -> legend
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        data = data,
                        surv.median.line = "hv",
                        title = paste(title,", p =", signif(kmp, 2)), # change it when doing diff data
                        xlab = "Survival in days",
                        ylab = 'Probability of survival',
                        # legend.title = "Methyla group:",
                        legend= c(0.75,0.85),
                        # ggtheme = theme_survminer(),
                        ggtheme = theme(
                          panel.border = element_blank(), panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", 
                                                                                       size = 0.5), 
                          panel.background = element_rect(fill = "white"),
                          legend.key = element_blank(),
                          legend.background = element_blank(),
                          text = element_text(size = 10, colour = "black"),
                          legend.text = element_text(size = 8, colour = "black"),
                          axis.text = element_text(size = 10, colour = "black"),
                          legend.title = element_blank(),
                          axis.title = element_text(size = 10,color = "black")
                        )
  ) +
    scale_color_manual(
      values = color,
      labels = legend$label
    )
  ggsave(filename =sur_name, path = result_path,device = "pdf",height = h,width = w)
  return(kmp)
}

# calculation of the pvalue and output the survival plots
gene_list_snv.hypermutation_class %>%
  dplyr::mutate(survival = purrr::map2(cancer_types,res,fn_survival)) %>%
  dplyr::select(-res) -> cancers_ICP_mutation_survival

# filter
cancers_ICP_mutation_survival %>%
  tidyr::unnest() %>%
  dplyr::filter(kmp <= 0.05)

cancers_ICP_mutation_survival %>%
  tidyr::unnest() %>%
  readr::write_tsv(file.path(snv_path,"cancers_ICP_mutation_survival.tsv"))

save.image(file = file.path(snv_path, ".rda_02_snv_a_gene_list_sruvival_class.rda"))
load(file = file.path(snv_path, ".rda_02_snv_a_gene_list_sruvival_class.rda"))
