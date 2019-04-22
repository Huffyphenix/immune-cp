############################
# GSVA score
# all 79 ICP genes as one feature
############################
library(GSVA)
library(magrittr)
library(ggplot2)
library(survival)

# data path ---------------------------------------------------------------
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"ICP_score/2.GSVA-ICPs_exp_site_5_feature")

# load data ---------------------------------------------------------------
exp_data <- readr::read_rds(file.path(TCGA_path,"pancan33_expr.rds.gz"))
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header = T)
gene_list$symbol <- as.character(gene_list$symbol)
survival_data <- readr::read_rds(file.path("/home/huff/project/data/TCGA-survival-time/cell.2018.survival","TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>%
  dplyr::rename("cancer_types" = "type")

ICPs_exp_in_TI.cancer.site_groups <- 
  readr::read_tsv(file.path(res_path, "ICPs_exp_in_TI.cancer.site_groups.tsv")) %>%
  dplyr::inner_join(gene_list,by = "symbol") %>%
  dplyr::select(symbol, Tissues, `log2FC(I/T)`, TCGA_Cancer, ICP_exp_group, GeneID, type, functionWithImmune)

source("/home/huff/project/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

# calculation of GSVA score -----------------------------------------------
fn_get_gene_list_feature   <- function(cancer_types){
  ICPs_exp_in_TI.cancer.site_groups %>%
    dplyr::filter(TCGA_Cancer %in% cancer_types) -> .data
  
  .group <- unique(.data$ICP_exp_group)
  
  .genelist <- list()
  for (group in .group){
    .data %>%
      dplyr::filter(ICP_exp_group %in% group) %>%
      .$GeneID -> .geneid
    .genelist[[group]] <- .geneid
  }
  .genelist
}

fn_GSVA <- function(cancer_types, exp){
  print(paste("GSVA",cancer_types))
  
  genelist <- fn_get_gene_list_feature(cancer_types)
  if (length(genelist)>0) {
    index <- which(substr(colnames(exp),14,14) == "0")
    exp <- exp[,c(2,index)]
    
    exp <- as.matrix(exp)
    rownames.exp <- exp[,1]
    exp <- exp[,-1]
    exp <- apply(exp, 2, as.numeric)
    rownames(exp) <- rownames.exp
    res.gsva <- gsva(exp,genelist, mx.diff = FALSE, 
                     method = c("gsva"),
                     kcdf = c("Gaussian"), 
                     verbose = FALSE, parallel.sz = 1)
    res.gsva %>%
      as.data.frame() %>%
      dplyr::as.tbl() %>%
      dplyr::mutate(feature = rownames(res.gsva)) %>%
      tidyr::gather(-feature, key = "barcode",value = "GSVA_score") %>%
      tidyr::spread(key = "feature", value = "GSVA_score") -> gsva.score
  } else{
    gsva.score <- tibble::tibble()
  }
  
  
  # fn_compare_TN(gsva.score, cancer_types, result_path = file.path(res_path,"TN_compare"))
  
  # fn_survival.calculate(gsva.score, cancer_types)
  
  # fn_heatmap(res.gsva)
  gsva.score
}

library(ComplexHeatmap)
library(pheatmap)

fn_heatmap <- function(res.gsva){
  pheatmap(res.gsva,scale = "row",
           clustering_distance_rows = "correlation",
           color = colorRampPalette(c("seagreen","white","red"))(50),
           border_color = NA,cutree_rows = 2,
           cutree_cols = 2,main = "DE_TF_expression")
  
  he = Heatmap(res.gsva,
               show_row_names = FALSE, 
               show_column_names = FALSE,
               cluster_columns = FALSE,
               # top_annotation = sample_anno,
               show_row_dend = FALSE # whether show row clusters.
               )
}
fn_compare_TN <- function(gsva.score, cancer_types, result_path){
  print(paste("Compare TN",cancer_types))
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
    ggsave(filename = paste(fig_name,"png",sep = "."), path = result_path,device = "png",width = 4,height = 3)
    ggsave(filename = paste(fig_name,"pdf",sep = "."), path = result_path,device = "pdf",width = 4,height = 3)
  }
}

fn_survival.calculate <- function(gsva.score, cancer_types){
  print(paste("Survival",cancer_types))
  if (cancer_types == "LAML") {
    gsva.score %>%
      dplyr::filter(substr(barcode,14,15) == "03") %>%
      dplyr::mutate(barcode =  substr(barcode,1,12)) -> gsva.score.T
  } else{
    gsva.score %>%
      dplyr::filter(substr(barcode,14,15) == "01") %>%
      dplyr::mutate(barcode =  substr(barcode,1,12)) -> gsva.score.T
  }
  
  ct <- cancer_types
  survival_data %>%
    dplyr::filter(cancer_types == ct) %>%
    tidyr::unnest() %>%
    dplyr::select(bcr_patient_barcode, PFI.1, PFI.time.1) %>%
    dplyr::rename("barcode" = "bcr_patient_barcode","status" = "PFI.1","time" = "PFI.time.1") %>%
    dplyr::mutate(time = time/365) %>%
    dplyr::inner_join(gsva.score.T, by = "barcode") %>%
    dplyr::mutate(group = ifelse(GSVA_score > quantile(GSVA_score,0.5), "High", "Low")) %>%
    dplyr::filter(!is.na(time)) -> .survival
  # dplyr::filter(time <= 1825) %>%
  # 
  
  # do survival diff analysis, get pvalue 
  if (length(unique(.survival$group)) == 2) {
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
}

exp_data %>%
  # head(1) %>%
  dplyr::mutate(GSVA = purrr::map2(cancer_types, expr, fn_GSVA)) %>%
  dplyr::select(-expr) -> GSVA.score

GSVA.score %>%
  readr::write_rds(file.path(res_path, "ICP_GSVA_score.rds.gz"), compress = "gz")
