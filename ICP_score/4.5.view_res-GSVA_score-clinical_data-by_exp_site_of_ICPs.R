############################
# GSVA score
# 79 ICP genes classified into groups by exp site, and these groups are used as features to do GSVA
# result display
############################

library(ComplexHeatmap)
library(pheatmap)
library(ConsensusClusterPlus)
library(magrittr)
library(ggplot2)
library(survival)
library(survminer)

# data path ---------------------------------------------------------------
# server 1 basic path
basic_path <- file.path("/home/huff/project")

# E zhou basi path
basic_path <- file.path("F:/我的坚果云")

immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
# TCGA_path <- file.path("/data/TCGA/TCGA_data")
# gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"ICP_score/2.1.Clinical_validation-GSVA-ICPs_exp_site_5_feature")

GSVA.score <- readr::read_rds(file.path(res_path, "ICP_GSVA_score.rds.gz"))

# survival data of PFI
# server 1
survival_path <- "/home/huff/project/data/TCGA-survival-time/cell.2018.survival"
# e zhou
survival_path <- file.path(basic_path,"immune_checkpoint/clinical_response_data")

survival_data <- readr::read_tsv(file.path(survival_path, "RNAseq-sample_info_complete.tsv")) %>%
  dplyr::filter( Library_strategy == "RNA-Seq") %>%
  tidyr::nest(-`pubmed ID`, .key = "clinical")

# data from TIMER
# server 1
# immunity_path_2 <- "/home/huff/project/immune_checkpoint/data/immunity"
# #  e zhou
# immunity_path_2 <- res_path
# TIMER_immunity_onlyTumor <- readr::read_tsv(file.path(immunity_path_2,"immuneEstimation.txt")) %>%
#   dplyr::filter(substr(barcode,14,14) == 0) %>%
#   dplyr::mutate(barcode = substr(barcode,1,12)) %>%
#   tidyr::gather(-barcode, key = "Cell_type", value = "TIL")
# 
# # data from Mutation burden
# # server 1
# MB_path_2 <- "/home/huff/project/data/TCGA"
# #  e zhou
# MB_path_2 <- res_path
# mutation_burden_class <- readr::read_rds(file.path(MB_path_2,"classfication_of_26_cancers_by_mutation_burden192.rds.gz")) %>%
#   tidyr::unnest() %>%
#   dplyr::rename("cancer_types" = "Cancer_Types")

source("F:/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

# draw pic p---------------------------------------------------------------
fn_heatmap <- function(pubmed_ID, Cancer.y, tissue, GSVA, clinical){
  cancer_types <- Cancer.y
  if(nrow(GSVA)>0){
    res.gsva <- as.data.frame(GSVA) %>% t() 
    colnames.res.gsva <- res.gsva[1,]
    rownames.res.gsva <- rownames(res.gsva)[-1]
    res.gsva <- res.gsva[-1,]
    res.gsva <- apply(res.gsva, 2, as.numeric)
    rownames(res.gsva) <- rownames.res.gsva
    colnames(res.gsva) <- colnames.res.gsva
    
    ## pheatmap
    pdf(file.path(res_path,paste(pubmed_ID, cancer_types,"GSVAscore", "heatmap.pdf", sep = "_")),
        height = 3, width = 6)
    pheatmap(res.gsva, #scale = "row",
             clustering_distance_rows = "correlation",
             color = colorRampPalette(c("#00B2EE","white","red"))(50),
             border_color = NA,cutree_rows = 2,
             show_colnames = F, treeheight_col = 30, treeheight_row = 20,
             cutree_cols = 2,main = paste(cancer_types," GSVA score", sep = ","))
    dev.off()

    tiff(file.path(res_path,paste(pubmed_ID, cancer_types,"GSVAscore", "heatmap.tiff", sep = "_")),
         height = 300, width = 600,compression = c("lzw"))
    pheatmap(res.gsva, #scale = "row",
             clustering_distance_rows = "correlation",
             color = colorRampPalette(c("#00B2EE","white","red"))(50),
             border_color = NA,cutree_rows = 2,
             show_colnames = F, treeheight_col = 30, treeheight_row = 20,
             cutree_cols = 2,main = paste(paste("PMID:",pubmed_ID),cancer_types," GSVA score", sep = ","))
    dev.off()
    
    ## do consensus cluster
    # results = ConsensusClusterPlus(res.gsva,maxK = 5,reps = 100,pItem = 0.8,
    #                                pFeature = 1,  innerLinkage="complete", finalLinkage="complete",
    #                                title = file.path(res_path, cancer_types),
    #                                distance = "euclidean", plot = "pdf",
    #                                clusterAlg = "hc",seed = 1262118388.71279)
    # 
    # results %>%
    #   readr::write_rds(file.path(res_path, pubmed_ID, paste(cancer_types, "CC_res.rds.gz")), compress = "gz")
    
    ## hclust
    # Compute the dissimilarity matrix
    gsva.dist <- dist(t(res.gsva), method = "euclidean")
    gsva.hc <- hclust(d = gsva.dist, method = "complete")
    # library("factoextra")
    # fviz_dend(gsva.hc, cex = 0.5)
    
    # Compute cophentic distance
    # gsva.coph <- cophenetic(gsva.hc)
    
    # Correlation between cophenetic distance and
    # the original distance
    # cor(gsva.dist, gsva.coph)
    
    # Cut tree into 4 groups
    # group <- cutree(gsva.hc, k = 4)
    # head(grp, n = 4)
    
    
    
    ## draw pic 
    for (i in 2:6) {
      C <- i
      group <- cutree(gsva.hc, k = C)
      # data.frame(Run = names(group),group = group) %>%
      #   readr::write_tsv(file.path(res_path, pubmed_ID, paste(C,"clusters_group_info.tsv",sep = "_"))) 
      
      data.frame(Run = names(group),group = group) %>%
        dplyr::as.tbl()  -> group
      
      ###### complex heatmap #####
      clinical %>%
        dplyr::inner_join(group, by = "Run") %>%
        dplyr::select(Response, group) %>%
        dplyr::mutate(group = as.character(group)) %>%
        as.data.frame() -> sample_info
      rownames(sample_info) <- clinical$Run
      
      sample_anno <- HeatmapAnnotation(df = sample_info,
                                       col = list(Response=c("CR" = "#228B22",
                                                             "PR" =  "#C1FFC1",
                                                             "SD" = "#969696",
                                                             "PD" = "black"),
                                                  group=color_list),
                                       width = unit(0.5, "cm"),
                                       name = c("Response", "CC group"))
      # draw(sample_anno,1:45)
      
      # heatmap
      png(file.path(res_path,"complex_heatmap",paste(C, pubmed_ID, "GSVAscore", "heatmap.png", sep = "_")),
          height = 300, width = 600)
      Heatmap(res.gsva,
              show_row_names = T,
              show_column_names = FALSE,
              cluster_columns = T,
              # clustering_distance_columns = "euclidean",
              # clustering_method_columns = "complete",
              top_annotation = sample_anno,
              show_row_dend = FALSE, # whether show row clusters.
              heatmap_legend_param = list(title = c("GSVA score")))
      dev.off()
      
      #### PCA analysis #####
      fn_pca_scatter(clinical=clinical,group=group,gsva.dist,cancer_types,pubmed_ID=pubmed_ID,C=C)
      
      # fn_compare(group = group, cancer_types = cancer_types,
      #            data = mutation_burden_class %>%
      #              dplyr::rename("Cell_type" = "cancer_types") %>%
      #              dplyr::mutate(sm_count = ifelse(sm_count==0,0.1,sm_count)) %>%
      #              dplyr::mutate(sm_count = log2(sm_count)),
      #            value = "sm_count", facet = "Cell_type", ylab = "log2 (Mutation burden)",
      #            title  = "", data_type = "MutationBurden", C = C,
      #            result_path = file.path(res_path, cancer_types), h = 3, w = 4)
      # fn_compare(group = group, cancer_types = cancer_types, data = TIMER_immunity_onlyTumor,
      #            value = "TIL", facet = "Cell_type", ylab = "TIL", C = C,
      #            title  = cancer_types, data_type = "TIMER",result_path = file.path(res_path, cancer_types), h = 6, w = 8)
      #### survival analysis #####
      fn_survival.calculate(group, clinical, pubmed_ID, C,cancer_types, result_path = res_path, h = 3, w = 4)
    }
  }
}
color_list <- c("1" = "pink1",
                "2" = "skyblue1",
                "3" = "darkseagreen1",
                "4" = "darkgreen",
                "5" = "dodgerblue4",
                "6" = "tan2")

library(ade4)
library(factoextra)
fn_pca_scatter <- function(clinical, group, gsva.dist, cancer_types,pubmed_ID,C){
  ##### scatter plot######
  
  group %>%
    dplyr::inner_join(clinical, by = "Run") %>%
    dplyr::select(Run, group, Response) -> group.info
  # group.info %>% dplyr::group_by(group, Response) %>% dplyr::mutate(n=n()) %>% dplyr::select(-Run) %>% unique() %>% dplyr::arrange(group)
  
  #### do PCA for distance data
  mds <- gsva.dist %>%          
    cmdscale() %>%
    tibble::as.tibble() %>%
    dplyr::rename("Dim.1" = "V1", "Dim.2" = "V2") %>%
    dplyr::mutate(Run = colnames(as.matrix(gsva.dist))) %>%
    dplyr::inner_join(group.info, by = "Run") %>%
    dplyr::mutate(group = as.factor(group))
  ggscatter(mds, x = "Dim.1", y = "Dim.2", 
            color = "group",
            shape = "Response",
            palette = "jco",
            size = 1, 
            ellipse = TRUE,
            ellipse.type = "convex",
            repel = TRUE) +
    theme(text = element_text(size = 10))
  
  ggsave(file.path(res_path,"PCA",paste(C, pubmed_ID, cancer_types,"GSVAscore", "PCA.png", sep = "_")), device = "png", height = 3, width = 4)
  ggsave(file.path(res_path,"PCA",paste(C, pubmed_ID, cancer_types,"GSVAscore", "PCA.pdf", sep = "_")), device = "pdf", height = 3, width = 4)
}

fn_compare <- function(group, cancer_types, data, value, facet, ylab, title, data_type, result_path,C, h, w){
  
  print(paste("Compare",data_type,cancer_types))
  
  group %>%
    dplyr::inner_join(data, by = "barcode") %>%
    dplyr::rename("value" = value) -> plot_ready
  
  assign("group1",facet)
  color_paired <- tibble::tibble(group = c(1:6),
                                 color = c("pink1", "skyblue1", "darkseagreen1", "darkgreen", "dodgerblue4", "tan2")) %>%
    dplyr::inner_join(plot_ready, by = "group") %>%
    dplyr::select(group, color) %>%
    unique()
  anno_text <- plot_ready %>%
    dplyr::group_by(group,Cell_type) %>%
    dplyr::mutate(n = n(), y = min(value) - max(value)*0.05) %>%
    dplyr::select(group,Cell_type,n, y) %>%
    unique() %>%
    dplyr::ungroup() 
  
  # combn_matrix <- combn(sort(unique(plot_ready$group)),2)
  # comp_list <- list()
  # for(i in 1:ncol(combn_matrix)){
  #   comp_list[[i]] <- combn_matrix[,i]
  # }
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
      geom_text(aes(x = group,y = y,label = paste("n=",n)), data = anno_text, size = 2) +
      facet_wrap(as.formula(paste("~",facet)), strip.position = "bottom", scales = "free") +
      scale_color_manual(
        values = color_paired$color
      ) +
      scale_fill_manual(
        values = color_paired$color
      ) +
      # ylim(4,12) +
      labs(y = ylab, x = "Cluster", title  = title) +
      theme(legend.position = "none",
            axis.title = element_text(colour = "black"),
            strip.background = element_rect(fill = "white",colour = "white"),
            text = element_text(size = 12, colour = "black")) +
      ggpubr::stat_compare_means(method = "anova") 
    # ggpubr::stat_compare_means(comparisons = comp_list, method = "t.test",label = "p.signif")
    fig_name <- paste("C",C,cancer_types, data_type, sep = "_")
    result_path <- file.path(result_path, data_type)
    ggsave(filename = paste(fig_name,"tiff",sep = "."), path = result_path,device = "tiff",width = w,height = h)
    ggsave(filename = paste(fig_name,"pdf",sep = "."), path = result_path,device = "pdf",width = w,height = h)
  }
}

fn_survival.calculate <- function(group, clinical, pubmed_ID, C,cancer_types, result_path, h, w){
  print(paste("Compare","survival",cancer_types))
  
  ct <- cancer_types
  group %>%
    dplyr::inner_join(clinical, by = "Run") %>%
    dplyr::select(Run, group, Survival_time, Survival_status) %>%
    dplyr::filter(! is.na(Survival_time) & ! is.na(Survival_status)) %>%
    dplyr::rename("status" = "Survival_status","time" = "Survival_time") %>%
    dplyr::mutate(status = ifelse(status == "Alive", 0, 1))-> .survival
  
  
  # do survival diff analysis, get pvalue 
  if (length(unique(.survival$group)) >= 2) {
    .survival %>%
      dplyr::mutate(time = time/365) -> .survival
    fit <- survfit(survival::Surv(time, status) ~ group, data = .survival, na.action = na.exclude)
    diff <- survdiff(survival::Surv(time, status) ~ group, data = .survival, na.action = na.exclude)
    kmp <- 1 - pchisq(diff$chisq, df = length(levels(as.factor(.survival$group))) - 1)
    print(kmp)
    if (kmp <= 0.1) {
      title <- paste(cancer_types)
      color_list <- tibble::tibble(group = c(1:6),
                                   color = c("pink1", "skyblue1", "darkseagreen1", "darkgreen", "dodgerblue4", "tan2"))
      sur_name <- paste("C", C, pubmed_ID, "GSVA_middle-survival-5y", sep = "_")
      result_path <- file.path(res_path,"survival")
      fn_survival(.survival,title,color_list,"group",sur_name,xlab = "Time (years)",result_path,h=h,w=w,lx = 0.8,ly = 0.8)
    }
  }
}


GSVA.score %>%
  dplyr::inner_join(survival_data, by="pubmed ID") %>%
  dplyr::rename("pubmed_ID"="pubmed ID") %>%
  # tail(8) %>%
  purrr::pmap(.f=fn_heatmap)

i = 5
  GSVA.score %>%
    dplyr::inner_join(survival_data, by="pubmed ID") %>%
    dplyr::rename("pubmed_ID"="pubmed ID") %>% .$GSVA %>% .[[i]] -> GSVA
  GSVA.score %>%
    dplyr::inner_join(survival_data, by="pubmed ID") %>%
    dplyr::rename("pubmed_ID"="pubmed ID") %>% .$clinical %>% .[[i]] -> clinical
  GSVA.score %>%
    dplyr::inner_join(survival_data, by="pubmed ID") %>%
    dplyr::rename("pubmed_ID"="pubmed ID") %>% .$pubmed_ID %>% .[[i]] -> pubmed_ID
  
  res.gsva <- as.data.frame(GSVA) %>% t() 
  colnames.res.gsva <- res.gsva[1,]
  rownames.res.gsva <- rownames(res.gsva)[-1]
  res.gsva <- res.gsva[-1,]
  res.gsva <- apply(res.gsva, 2, as.numeric)
  rownames(res.gsva) <- rownames.res.gsva
  colnames(res.gsva) <- colnames.res.gsva
  
  # Compute the dissimilarity matrix
  gsva.dist <- dist(t(res.gsva), method = "euclidean")
  gsva.hc <- hclust(d = gsva.dist, method = "complete")
  
j=2
    C <- j
    group <- cutree(gsva.hc, k = C)
    data.frame(Run = names(group),group = group) %>%
      readr::write_tsv(file.path(res_path, pubmed_ID, paste(C,"clusters_group_info.tsv",sep = "_"))) 
    
    data.frame(Run = names(group),group = group) %>%
      dplyr::as.tbl()  -> group
    
    ###### complex heatmap #####
    clinical %>%
      dplyr::inner_join(group, by = "Run") %>%
      dplyr::select(Response, group) %>%
      dplyr::mutate(group = as.character(group)) %>%
      as.data.frame() -> sample_info
    rownames(sample_info) <- clinical$Run
    
    sample_anno <- HeatmapAnnotation(df = sample_info,
                                     col = list(Response=c("CR" = "#228B22",
                                                           "PR" =  "#C1FFC1",
                                                           "SD" = "#969696",
                                                           "PD" = "black",
                                                           "NE" = "white",
                                                           "X" = "white",
                                                           "PRCR" = "#C1FFC1",
                                                           "NR" = "black",
                                                           "R" = "#228B22"),
                                                group=color_list),
                                     width = unit(0.5, "cm"),
                                     name = c("Response", "CC group"))
    # draw(sample_anno,1:45)
    
    # heatmap
    png(file.path(res_path,"complex_heatmap",paste(C, pubmed_ID, "GSVAscore", "heatmap.png", sep = "_")),
        height = 300, width = 600)
    Heatmap(res.gsva,
            show_row_names = T,
            show_column_names = FALSE,
            cluster_columns = T,
            # clustering_distance_columns = "euclidean",
            # clustering_method_columns = "complete",
            top_annotation = sample_anno,
            show_row_dend = FALSE, # whether show row clusters.
            heatmap_legend_param = list(title = c("GSVA score")))
    dev.off()


# all data in one ---------------------------------------------------------

    GSVA.score %>%
      dplyr::inner_join(survival_data, by="pubmed ID") %>%
      dplyr::rename("pubmed_ID"="pubmed ID") %>%
      # tail(8) %>%
      purrr::pmap(.f=fn_heatmap)
    