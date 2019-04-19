############################
# GSVA score
############################
library(GSVA)
library(magrittr)

# data path ---------------------------------------------------------------
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"ICP_score")

# load data ---------------------------------------------------------------
exp_data <- readr::read_rds(file.path(TCGA_path,"pancan33_expr.rds.gz"))
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header = T)
gene_list$symbol <- as.character(gene_list$symbol)
survival_data <- readr::read_rds(file.path("/home/huff/project/data/TCGA-survival-time/cell.2018.survival","TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>%
  dplyr::rename("cancer_types" = "type")

source("/home/huff/project/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

# calculation of GSVA score -----------------------------------------------

fn_GSVA <- function(cancer_types, exp){
  genelist <- list(ICPs = gene_list$GeneID)
  exp <- as.matrix(exp)
  rownames.exp <- exp[,2]
  exp <- exp[,-c(1:2)]
  exp <- apply(exp, 2, as.numeric)
  rownames(exp) <- rownames.exp
  res.gsva <- gsva(exp,genelist, mx.diff = FALSE, 
                   method = c("gsva"),
                   kcdf = c("Gaussian"), 
                   verbose = FALSE, parallel.sz = 1)
  res.gsva %>%
    as.data.frame() %>%
    dplyr::as.tbl() %>%
    tidyr::gather(key = "barcode",value = "GSVA_score") -> gsva.score
  
  fn_compare_TN(gsva.score, cancer_types, result_path = file.path(res_path,"TN_compare"))
  
  fn_survival(gsva.score)
}

fn_compare_TN <- function(gsva.score, cancer_types, result_path){
  gsva.score %>%
    dplyr::mutate(group = ifelse(substr(barcode,14,15) == "01","Tumor","Not_primary")) %>%
    dplyr::mutate(group = ifelse(substr(barcode,14,15) == "11","Normal",group)) %>%
    dplyr::filter(group != "Not_primary") -> gsva.score.TN
  color_paired <- tibble::tibble(group = c("Normal", "Tumor"),
                          color = c("#00B2EE", "#CD2626")) %>%
    dplyr::inner_join(gsva.score.TN, by = "group") %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(n = n(), y = min(GSVA_score) - 0.05) %>%
    dplyr::select(-barcode,-GSVA_score) %>%
    unique()
    
  comp_list <- list(c("Tumor", "Normal"))
  if (nrow(color_paired) == 2) {
    gsva.score.TN %>%
      ggpubr::ggboxplot(x = "group", y = "GSVA_score", fill = "white",alpha = 0,width = 0.1,
                        color = "group" #add = "jitter",#, palette = "npg"
      ) +
      geom_violin(aes(fill = group),alpha = 0.5) +
      # geom_jitter(aes(color = group),alpha = 0.2,width = 0.1,size=0.1) +
      geom_boxplot(fill = "white",alpha = 0,width = 0.1) +
      geom_text(aes(x = group,y = y,label = paste("n=",n)), data = color_paired) +
      scale_color_manual(
        values = color_paired$color
      ) +
      scale_fill_manual(
        values = color_paired$color
      ) +
      # ylim(4,12) +
      labs(y = "GSVA score", title  = cancer_types) +
      theme(legend.position = "none",
            axis.title = element_text(colour = "black"),
            axis.title.x = element_blank(),
            strip.background = element_rect(fill = "white",colour = "white"),
            text = element_text(size = 8, colour = "black")) +
      # ggpubr::stat_compare_means(method = "t.test") +
      ggpubr::stat_compare_means(comparisons = comp_list,method = "t.test",label = "p.signif")
    fig_name <- paste(cancer_types,"GSVA_score-T_N")
    ggsave(filename = paste(fig_name,"png",sep = "."), path = result_path,device = "png",width = w,height = h)
    ggsave(filename = paste(fig_name,"pdf",sep = "."), path = result_path,device = "pdf",width = w,height = h)
  }
}

fn_survival.calculate <- function(gsva.score, cancer_types, result_path){
  gsva.score %>%
    dplyr::filter(substr(barcode,14,15) == "01") %>%
    dplyr::mutate(barcode =  substr(barcode,1,12)) -> gsva.score.T
  
  ct <- cancer_types
  survival_data %>%
    dplyr::filter(cancer_types == ct) %>%
    tidyr::unnest() %>%
    dplyr::select(bcr_patient_barcode, PFI.1, PFI.time.1) %>%
    dplyr::rename("barcode" = "bcr_patient_barcode","status" = "PFI.1","time" = "PFI.time.1") %>%
    dplyr::mutate(time = time/365) %>%
    dplyr::inner_join(gsva.score.T, by = "barcode") %>%
    dplyr::mutate(group = ifelse(GSVA_score > quantile(GSVA_score,0.5), "High", "Low")) -> .survival
    # dplyr::filter(time <= 1825) %>%
    # 
  
  # do survival diff analysis, get pvalue 
  fit <- survfit(survival::Surv(time, status) ~ group, data = .survival, na.action = na.exclude)
  diff <- survdiff(survival::Surv(time, status) ~ group, data = .survival, na.action = na.exclude)
  kmp <- 1 - pchisq(diff$chisq, df = length(levels(as.factor(.survival$group))) - 1)
  
  if (kmp <= 0.05) {
    title <- paste(cancer_types)
    color_list <- tibble::tibble(group = c("Low", "High"),
                                 color = c("#00B2EE", "#CD2626"))
    sur_name <- paste(cancer_types, "GSVA_middle-survival", sep = "_")
    result_path <- file.path(res_path,"survival")
    fn_survival(.survival,title,color_list,"group",sur_name,xlab = "Time (years)",result_path,3,4,lx = 0.8,ly = 0.8)
  }
}

exp_data %>%
  head(1) %>%
  dplyr::mutate(GSVA = purrr::map2(cancer_types, expr, fn_GSVA)) %>%
  dplyr::select(-expr) -> GSVA.score

GSVA.score %>%
  readr::write_rds(file.path(res_path, "ICP_GSVA_score.rds.gz"))
