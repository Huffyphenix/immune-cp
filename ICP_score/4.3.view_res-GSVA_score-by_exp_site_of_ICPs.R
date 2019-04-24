############################
# GSVA score
# all 79 ICP genes as one feature
# result display
############################

library(ComplexHeatmap)
library(pheatmap)
library(ConsensusClusterPlus)

# data path ---------------------------------------------------------------
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"ICP_score/2.GSVA-ICPs_exp_site_5_feature")

GSVA.score <- readr::read_rds(file.path(res_path, "ICP_GSVA_score.rds.gz"), compress = "gz")

# survival data of PFI
survival_data <- readr::read_rds(file.path("/home/huff/project/data/TCGA-survival-time/cell.2018.survival","TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>%
  dplyr::rename("cancer_types" = "type")

# data from TIMER
immunity_path_2 <- "/home/huff/project/immune_checkpoint/data/immunity"
TIMER_immunity_onlyTumor <- readr::read_tsv(file.path(immunity_path_2,"immuneEstimation.txt")) %>%
  dplyr::filter(substr(barcode,14,14) == 0) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  tidyr::gather(-barcode, key = "Cell_type", value = "TIL")

mutation_burden_class <- readr::read_rds(file.path("/home/huff/project/data/TCGA","classfication_of_26_cancers_by_mutation_burden192.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::rename("cancer_types" = "Cancer_Types") %>%

# draw pic p---------------------------------------------------------------
fn_heatmap <- function(gsva.score, cancer_types){
  res.gsva <- as.data.frame(gsva.score) %>% t() 
  colnames.res.gsva <- res.gsva[1,]
  rownames.res.gsva <- rownames(res.gsva)
  res.gsva <- res.gsva[-1,]
  res.gsva <- apply(res.gsva, 2, as.numeric)
  rownames(res.gsva) <- rownames.res.gsva
  colnames(res.gsva) <- colnames.res.gsva
  
  p <- pheatmap(res.gsva, #scale = "row",
           clustering_distance_rows = "correlation",
           color = colorRampPalette(c("#00B2EE","white","red"))(50),
           border_color = NA,cutree_rows = 2,
           show_colnames = F, treeheight_col = 30, treeheight_row = 20,
           cutree_cols = 2,main = paste(cancer_types," GSVA score", sep = ","))
  pdf(file.path(res_path,cancer_types,paste(cancer_types,"GSVAscore", "heatmap.pdf", sep = "_")),
      height = 3, width = 6)
  p
  dev.off()
  
  tiff(file.path(res_path,cancer_types,paste(cancer_types,"GSVAscore", "heatmap.tiff", sep = "_")),
      height = 300, width = 600,compression = c("lzw"))
  p
  dev.off()
  results = ConsensusClusterPlus(res.gsva,maxK = 5,reps = 100,pItem = 0.8,
                                 pFeature = 1, 
                                 title = file.path(res_path, cancer_types),
                                 distance = "pearson", plot = "pdf",
                                 clusterAlg = "hc",seed = 1262118388.71279)
  
  results %>%
    readr::write_rds(file.path(res_path, cancer_types, paste(cancer_types, "CC_res.rds.gz")), compress = "gz")
  for (i in 2:5) {
    C <- i
    group <- results[[C]][['consensusClass']]
    data.frame(sample = names(group),group = group) %>%
      readr::write_tsv(file.path(res_path, cancer_types, paste(C,"clusters_group_info.tsv",sep = "_"))) 
    
    data.frame(barcode = names(group),group = group) %>%
      dplyr::as.tbl() %>%
      dplyr::mutate(barcode = substr(barcode, 1, 12)) -> group
    
    fn_compare(group = group, cancer_types = cancer_types, data = TIMER_immunity_onlyTumor, 
               value = "sm_count", facet = "~cancer_types", ylab = "log2 (Mutation burden)",
               title  = "", result_path = file.path(res_path, cancer_types))
    fn_compare(group, cancer_types, mutation_burden_class, "sm_count",file.path(res_path, cancer_types))
    
    fn_survival.calculate(group, cancer_types)
  }

}

fn_compare <- function(group, cancer_types, data, value, facet, ylab, title,  result_path){

  print(paste("Compare TN",cancer_types))
  group %>%
    dplyr::inner_join(data, by = "barcode") %>%
    dplyr::rename("value" = value) -> plot_ready
  color_paired <- tibble::tibble(group = c(1:5),
                                 color = c("pink1", "skyblue1", "darkseagreen1", "darkgreen", "dodgerblue4")) %>%
    dplyr::inner_join(plot_ready, by = "group") %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(n = n(), y = min(value) - max(value)*0.05) %>%
    dplyr::select(group, color, n, y) %>%
    unique()
  
  # comp_list <- list(c("Tumor", "Normal"))
  if (nrow(color_paired) >= 2) {
    plot_ready %>%
      dplyr::arrange(group) %>%
      ggpubr::ggboxplot(x = "group", y = "value", fill = "white",alpha = 0,width = 0.1,
                        color = "group" #add = "jitter",#, palette = "npg"
      ) +
      geom_violin(aes(fill = group),alpha = 0.5) +
      # geom_jitter(aes(color = group),alpha = 0.2,width = 0.1,size=0.1) +
      geom_boxplot(fill = "white",alpha = 0,width = 0.1) +
      geom_text(aes(x = group,y = y,label = paste("n=",n)), data = color_paired) +
      facet_wrap(as.formula(facet), strip.position = "bottom", scales = "free") +
      scale_color_manual(
        values = color_paired$color
      ) +
      scale_fill_manual(
        values = color_paired$color
      ) +
      # ylim(4,12) +
      labs(y = ylab, title  = title) +
      theme(legend.position = "none",
            axis.title = element_text(colour = "black"),
            axis.title.x = element_blank(),
            strip.background = element_rect(fill = "white",colour = "white"),
            text = element_text(size = 8, colour = "black")) +
      # ggpubr::stat_compare_means(method = "t.test") +
      ggpubr::stat_compare_means(method = "t.test",label = "p.signif")
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