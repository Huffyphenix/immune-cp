
# verify the immune ration using published data ---------------------------
library(magrittr)

# data path -----
xcell_path <- "/project/huff/huff/data/TCGA/immune_infiltration/xCell"
purity_path <- "/project/huff/huff/data/TCGA_tumor_purity_from_ncomms9971"
immune_path <- "/project/huff/huff/immune_checkpoint"
immunity_path <- "/project/huff/huff/immune_checkpoint/data/immunity"

# load data -----
# immune ratio we calculated
xCell_TCGA_RSEM.immune_stroma.ratio <- readr::read_rds(file.path(immune_path,"genelist_data","Pancan21.tcga.xCell.immune_stroma.ratio.rds.gz")) %>%
  dplyr::select(-barcode) %>%
  dplyr::rename("barcode" = "Sample ID")

# TIMER data
TIMER_immunity <- readr::read_tsv(file.path(immunity_path,"immuneEstimation.txt")) %>%
  dplyr::mutate(TIL = B_cell+CD4_Tcell+CD8_Tcell+Neutrophil+Macrophage+Dendritic)

# purity data
TCGA_purity <- readr::read_tsv(file.path(purity_path,"ncomms9971-s2.txt")) %>%
  dplyr::select(`Cancer type`,`Sample ID`,`CPE`) %>%
  dplyr::mutate(barcode=substr(`Sample ID`,1,15))

# data combination ----

TIMER_immunity %>%
  dplyr::inner_join(xCell_TCGA_RSEM.immune_stroma.ratio,by="barcode") -> immunity_combine_xCell.immnue.ratio


# distribution ------------------------------------------------------------
xCell_TCGA_RSEM.immune_stroma.ratio %>%
  dplyr::filter(!is.na(Immune_ratio)) %>%
  dplyr::group_by(`Cancer type`) %>%
  dplyr::mutate(mid = quantile(Immune_ratio,0.5)) %>%
  dplyr::select(`Cancer type`,mid) %>%
  unique() %>%
  dplyr::arrange(mid) -> cancer_rank
xCell_TCGA_RSEM.immune_stroma.ratio %>%
  dplyr::filter(!is.na(Immune_ratio)) %>%
  dplyr::select(barcode,`Cancer type`,CPE,Immune_ratio,Stromal_ratio,other_Stromal_ratio) %>%
  tidyr::gather(-barcode,-`Cancer type`,key="Cell Type",value="Ratio") %>%
  ggplot(aes(x=`Cancer type`,y=Ratio)) +
  geom_violin(aes(fill = `Cancer type`),alpha = 0.5) +
  facet_grid(~`Cell Type`,scales = "free") +
  scale_x_discrete(limit = cancer_rank$`Cancer type`) +
  coord_flip() +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(vjust = 1, hjust = 1, size = 10),
    axis.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12,angle = 90),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(size = 20),
    axis.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "white",colour = "black")
  )
  
ggsave(file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","ImmuneRatio_distribution_in_cancers.pdf"),device = "pdf",height = 4,width = 6)
ggsave(file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","ImmuneRatio_distribution_in_cancers.png"),device = "png",height = 4,width = 6)

# calculate correlation ---------------------------------------------------

fn_correlation <- function(.x){
  broom::tidy(cor.test(.x$Immune_ratio,.x$TIL,method = "pearson")) %>%
    dplyr::mutate(x="Immune_ratio",y="Total TIL") -> .tmp1
  broom::tidy(cor.test(.x$Immune_ratio,.x$B_cell,method = "pearson")) %>%
    dplyr::mutate(x="Immune_ratio",y="B_cell") -> .tmp2
  broom::tidy(cor.test(.x$Immune_ratio,.x$CD4_Tcell,method = "pearson")) %>%
    dplyr::mutate(x="Immune_ratio",y="CD4_Tcell") -> .tmp3
  broom::tidy(cor.test(.x$Immune_ratio,.x$CD8_Tcell,method = "pearson")) %>%
    dplyr::mutate(x="Immune_ratio",y="CD8_Tcell") -> .tmp4
  broom::tidy(cor.test(.x$Immune_ratio,.x$Neutrophil,method = "pearson")) %>%
    dplyr::mutate(x="Immune_ratio",y="Neutrophil") -> .tmp5
  broom::tidy(cor.test(.x$Immune_ratio,.x$Macrophage,method = "pearson")) %>%
    dplyr::mutate(x="Immune_ratio",y="Macrophage") -> .tmp6
  broom::tidy(cor.test(.x$Immune_ratio,.x$Dendritic,method = "pearson")) %>%
    dplyr::mutate(x="Immune_ratio",y="Dendritic") -> .tmp7
  rbind(.tmp1,.tmp2,.tmp3,.tmp4,.tmp5,.tmp6,.tmp7)
}

immunity_combine_xCell.immnue.ratio %>%
  tidyr::nest(-`Cancer type`) %>%
  dplyr::group_by(`Cancer type`) %>%
  dplyr::mutate(cor = purrr::map(data,fn_correlation)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::ungroup() -> correlation_TIMER.our_results

immunity_combine_xCell.immnue.ratio %>%
  dplyr::mutate(`Cancer type` = "Total") %>%
  tidyr::nest(-`Cancer type`) %>%
  dplyr::group_by(`Cancer type`) %>%
  dplyr::mutate(cor = purrr::map(data,fn_correlation)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::ungroup() -> Total.correlation_TIMER.our_results
# draw plot -----
library(ggplot2)
CPCOLS <- c("red", "white", "#1C86EE")
Total.correlation_TIMER.our_results %>%
  rbind(correlation_TIMER.our_results) %>%
  ggplot(aes(y=`Cancer type`,x=y)) +
  geom_tile(aes(fill=estimate),colour = "white") +
  scale_fill_gradient2(
    name = "Correlation", #"Methylation diff (T - N)",
    low = CPCOLS[3],
    high = CPCOLS[1],
    mid = CPCOLS[2],
    breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)
  ) +
  guides(fill=guide_colorbar(title.position="left",
                             title.theme = element_text(angle = 90))) +
  scale_y_discrete(limit = c(sort(correlation_TIMER.our_results$`Cancer type` %>% unique()),"Total")) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(vjust = 1, hjust = 1, angle = 90, size = 10),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12,angle = 90),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(size = 20),
    axis.text = element_text(colour = "black")
  ) -> p
p
ggsave(file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","TIMER","heatmap_correlation.TIMER.ourImmuneRatio.pdf"),device = "pdf",height = 6,width = 4)
ggsave(file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","TIMER","heatmap_correlation.TIMER.ourImmuneRatio.png"),device = "png",height = 6,width = 4)

Total.correlation_TIMER.our_results %>%
  rbind(correlation_TIMER.our_results) %>%
  readr::write_tsv(file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","TIMER","person.correlation_TIMER.our_results"))

# survival and immune ratio -----------------------------------------------
# sigle factor survival analysis ,OS
clinical <- readr::read_rds(file.path("/project/huff/huff/data/TCGA-survival-time/cell.2018.survival","TCGA_pancan_cancer_cell_survival_time.rds.gz"))
clinical_tcga <- readr::read_rds(file.path("/project/huff/huff/TCGA_survival/data","Pancan.Merge.clinical.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::select(-cancer_types) %>%
  unique() %>%
  dplyr::mutate(OS=as.numeric(OS),Status=as.numeric(Status),Age=as.numeric(Age))
clinical %>%
  tidyr::unnest() %>%
  dplyr::select(bcr_patient_barcode,PFS,PFS.time) %>%
  tidyr::drop_na() %>%
  dplyr::rename("barcode" = "bcr_patient_barcode")-> time_status

fn_survival <-function(.cancer,.x){
  .up=75
  .low=25
  .x %>%
    dplyr::filter(!is.na(Immune_ratio)) %>%
    dplyr::mutate(group = ifelse(Immune_ratio>quantile(Immune_ratio,.up/100),"High",NA)) %>%
    dplyr::mutate(group = ifelse(Immune_ratio<quantile(Immune_ratio,.low/100),"Low",group))-> .data
  .d_diff <- try(survival::survdiff(survival::Surv(OS, Status) ~ group, data = .data),silent = TRUE)
  if ('try-error' %in% class(.d_diff)) {return(data.frame(KMP=1,Coxp=1))}else{
    # print(.d_diff)
    kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.data$group))) - 1)
    # print(kmp)
    coxp <-  broom::tidy(survival::coxph(survival::Surv(OS, Status) ~ Immune_ratio, data = .data, na.action = na.exclude))
    # print(coxp)
    
    .fit <- survival::survfit(survival::Surv(OS, Status) ~ group, data = .data, na.action = na.exclude) 
    low_n=.fit$n[2]
    high_n=.fit$n[1]
    legend <- data.frame(group=c("Low","High"),n=c(low_n,high_n))
    legend %>%
      dplyr::mutate(
        label = purrr::map2(
          .x = group,
          .y = n,
          .f = function(.x,.y){
            latex2exp::TeX(glue::glue("IR^{<<.x>>}, n = <<.y>>", .open = "<<", .close = ">>"))
          }
        )
      ) %>%
      dplyr::arrange(group)-> legend
    survminer::ggsurvplot(.fit,pval=F, #pval.method = T,
                          surv.median.line = "hv",
                          title = paste(paste(.cancer, sep = "-"), "p =", signif(kmp, 2)),
                          xlab = "Time (days)",
                          ylab = 'Survival',
                          legend.title = "Expression group:",
                          legend= c(0.8,0.8),
                          ggtheme = theme(
                            panel.border = element_blank(), panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.5), 
                            panel.background = element_rect(fill = "white"),
                            legend.key = element_blank(),
                            legend.background = element_blank(),
                            legend.text = element_text(size = 12, colour = "black"),
                            axis.text = element_text(size = 12, colour = "black"),
                            legend.title = element_blank(),
                            axis.title = element_text(size = 12,color = "black")
                          ),
                          font.main = c(16),
                          font.x = c(14),
                          font.y = c(14),
                          font.tickslab = c(12)
    ) +
      scale_color_manual(
        values = c( "#EE6363","#1C86EE"),
        labels = legend$label
      )
    # dev.off()
    fig_name <- paste(.cancer,"PFS",signif(kmp, 2),"pdf", sep = ".")
    fig_name_png <- paste(.cancer,"PFS",signif(kmp, 2),"png", sep = ".")
    # ggsave(filename = fig_name, device = "pdf", path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","survival","75-25"), width = 4, height = 3)
    # ggsave(filename = fig_name_png, device = "png", path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","survival","75-25"), width = 4, height = 3)
    ggsave(filename = fig_name, device = "pdf", path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","survival","OS","75-25"), width = 4, height = 3)
    ggsave(filename = fig_name_png, device = "png", path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","survival","OS","75-25"), width = 4, height = 3)
    # print(.fit)
    # fn_sur_draw(.gene,kmp,.fit)
    data.frame(KMP=kmp,Coxp=coxp) -> .out
    return(.out)
  }
}

xCell_TCGA_RSEM.immune_stroma.ratio %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  tidyr::nest(-`Cancer type`) %>%
  dplyr::mutate(sur = purrr::map2(`Cancer type`,data,fn_survival)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() -> OS_50_ImmuneRatio_survival
OS_50_ImmuneRatio_survival %>%
  dplyr::filter(KMP <= 0.05)
OS_50_ImmuneRatio_survival %>% readr::write_tsv(file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","survival","OS","50","OS_50_ImmuneRatio_survival.kmp.tsv"))

xCell_TCGA_RSEM.immune_stroma.ratio %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  tidyr::nest(-`Cancer type`) %>%
  dplyr::mutate(sur = purrr::map2(`Cancer type`,data,fn_survival)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() -> OS_75.25_ImmuneRatio_survival
OS_75.25_ImmuneRatio_survival %>%
  dplyr::filter(KMP <= 0.05)
OS_75.25_ImmuneRatio_survival %>% readr::write_tsv(file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","survival","OS","75-25","OS_75.25_ImmuneRatio_survival.kmp.tsv"))


# IR between relapse or not, other form to represent PFS survival ---------
# logistic regression, refer to TIMER method

# prepare data
time_status %>%
  dplyr::mutate(Relapse = ifelse(PFS==0,"TUMOR FREE","WITH TUMOR")) %>%
  dplyr::inner_join(clinical_tcga,by="barcode") -> pfs_relapse

# prepare function
fn_logic_regression <- function(.x,.cancer){
  print(.cancer)
  if (length(unique(.x$Stage))==1) {
    model <- glm(Immune_ratio~PFS + Age_class,data=.x, family = binomial) #control = age+stage,
  }else if(length(unique(.x$Age))==1){
    model <- glm(Immune_ratio~PFS + Stage,data=.x, family = binomial) #control = age+stage,
  }else if (length(unique(.x$Stage))==1 && length(unique(.x$Age))==1){
    model <- glm(Immune_ratio~PFS,data=.x, family = binomial)
  }else{
    model <- glm(Immune_ratio~PFS + Age_class + Stage,data=.x, family = binomial) #control = age+stage,
  }
  coef <- summary(model)$coef
  call <- as.character(summary(model)$call)[2]
  p <- as.vector(coef[,"Pr(>|z|)"])
  q <- qvalue::qvalue(p,lambda=0)
  qval <- q$qvalues[1]
  fdr <- q$lfdr[1]
  tibble::tibble(call=call,pvalue=p[1],qvalue=qval,fdr=fdr)
}

# get results
xCell_TCGA_RSEM.immune_stroma.ratio %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(pfs_relapse,by="barcode") %>%
  dplyr::filter(!is.na(Stage)) %>%
  dplyr::filter(Age != "Not_applicable") %>%
  dplyr::mutate(Age_class = ifelse(Age<30,"Youth","Old")) %>%
  dplyr::mutate(Age_class = ifelse(Age>=30 & Age<50,"Middle",Age_class)) %>%
  tidyr::nest(-`Cancer type`) %>%
  dplyr::mutate(logic_immuneratio_relapse = purrr::map2(data, `Cancer type`,fn_logic_regression)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> immuneratio_relapse_logistic_regression

immuneratio_relapse_logistic_regression %>%
  readr::write_tsv(file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","relapse","immuneratio_relapse_logistic_regression.tsv"))

# draw plot
immuneratio_relapse_logistic_regression %>%
  dplyr::filter(fdr<0.15) -> immuneratio_relapse_logistic_regression.sig

xCell_TCGA_RSEM.immune_stroma.ratio %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(pfs_relapse,by="barcode") %>%
  dplyr::filter(`Cancer type` %in% immuneratio_relapse_logistic_regression.sig$`Cancer type`) %>%
  dplyr::filter(!is.na(Immune_ratio)) %>%
  dplyr::group_by(`Cancer type`) %>%
  dplyr::mutate(y=max(Immune_ratio)-0.15) %>%
  dplyr::mutate(x=2.4) %>%
  dplyr::select(`Cancer type`,x,y) %>%
  unique() %>%
  dplyr::inner_join(immuneratio_relapse_logistic_regression.sig,by="Cancer type") %>%
  dplyr::mutate(qvalue=paste("q=",signif(qvalue,2))) -> logistic_relapse_anno

xCell_TCGA_RSEM.immune_stroma.ratio %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(pfs_relapse,by="barcode") %>%
  dplyr::filter(`Cancer type` %in% immuneratio_relapse_logistic_regression.sig$`Cancer type`) %>%
  dplyr::filter(!is.na(Immune_ratio)) %>%
  ggplot(aes(x=Relapse,y=Immune_ratio)) +
  geom_violin(aes(fill = Relapse),alpha = 0.5) +
  geom_jitter(aes(color = Relapse),alpha = 0.2,width = 0.1) +
  geom_boxplot(width = 0.1,fill="white",alpha = 0) +
  geom_text(data=logistic_relapse_anno,aes(x=x,y=y,label=qvalue),size = 4, color="red") +
  facet_wrap(~`Cancer type`,scales = "free_x") +
  coord_flip() +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(vjust = 1, hjust = 1, size = 10),
    axis.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12,angle = 90),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(size = 20),
    axis.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "white",colour = "black"),
    strip.text = element_text(size=12)
  )
ggsave(filename = "immuneratio_with_relapse(PFS status)_logistic_regression.sig.png",device = "png", path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","relapse"), width = 6, height = 5)
ggsave(filename = "immuneratio_with_relapse(PFS status)_logistic_regression.sig.pdf",device = "pdf", path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","relapse"), width = 6, height = 5)

# cox model ---------------------------------------------------------------
# refer to TIMER method
fn_cox <- function(.x){
  covariates <- c("Immune_ratio","CPE", "Stromal_ratio",  "other_Stromal_ratio", "Age", "Stage_n")
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(OS, Status)~', x)))
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = .x)})
  
  # Extract data 
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR <- HR
                           HR.confint  <- paste0(" (", 
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR,HR.confint, wald.test, p.value)
                           names(res)<-c("beta", "HR","(95% CI for HR)", "wald.test", 
                                         "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  res.df <- as.data.frame(res)
  res.df %>%
    dplyr::as.tbl() %>%
    dplyr::mutate(covariates=covariates,n=nrow(.x)) %>%
    dplyr::select(covariates,beta:n)
}

xCell_TCGA_RSEM.immune_stroma.ratio %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::filter(!is.na(Stage)) %>%
  dplyr::filter(Stage!= "Not_applicable") %>%
  dplyr::filter(Age != "Not_applicable") %>%
  dplyr::filter(!is.na(Immune_ratio)) %>%
  dplyr::mutate(Stage_n = as.numeric(as.factor(Stage))) %>%
  # dplyr::mutate(Age_class = ifelse(Age<30,"Youth","Old")) %>%
  # dplyr::mutate(Age_class = ifelse(Age>=30 & Age<50,"Middle",Age_class)) %>%
  tidyr::nest(-`Cancer type`) %>%
  dplyr::mutate(cox_immuneratio = purrr::map(data, fn_cox)) %>%
  dplyr::select(-data) -> immuneratio_cox_model
fn_hrclass <- function(.x){
  if(.x<0.1){
    class <-"1" 
  }else if(.x>=0.1 & .x<0.3){
    class <-"2" 
  }else if(.x>=0.3 & .x<1){
    class <-"3" 
  }else if(.x>=1 & .x<3){
    class <-"4" 
  }else if(.x>=3 & .x<10){
    class <-"5" 
  }else if(.x>10){
    class <-"6" 
  }
  class
}

immuneratio_cox_model %>%
  tidyr::unnest() %>%
  dplyr::mutate(HR=as.numeric(HR),p.value=as.numeric(p.value)) %>%
  dplyr::mutate(HR_class = purrr::map(HR,fn_hrclass)) %>%
  tidyr::unnest() %>%
  dplyr::mutate(p = ifelse(p.value>0.05,"1","2")) %>%
  ggplot(aes(x=`Cancer type`,y=covariates)) +
  geom_point(aes(size=p,fill=HR_class),color="grey",shape=22) +
  scale_fill_manual(
    values=c("#104E8B", "#1874CD", "#87CEFF", "#FFC0CB", "#FF6A6A", "#FF0000"),
    labels = c("<0.1","0.1-0.3","0.3-1","1-3","3-10",">10")
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
ggsave(filename = "Cox_model_of_immuneratio.png",device = "png",path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","cox"),width = 5,height = 3)
ggsave(filename = "Cox_model_of_immuneratio.pdf",device = "pdf",path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","cox"),width = 5,height = 3)

immuneratio_cox_model %>%
  tidyr::unnest() %>%
  readr::write_tsv(file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","cox","Cox_model_of_immuneratio.tsv"))
