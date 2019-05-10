##############################################
# apply logistic model get from clinical data on TCGA samples
##############################################
library(magrittr)
library(tidyverse)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
theme_set(theme_classic())

# data path ---------------------------------------------------------------
# server 1 basic path
basic_path <- file.path("/home/huff/project")

# E zhou basi path
basic_path <- file.path("F:/我的坚果云")

# Hust basi path
basic_path <- file.path("S:/坚果云/我的坚果云")

immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
TCGA_path <- file.path("H:/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"ICP_score/4.logistics-model-on-TCGA")

# load data ---------------------------------------------------------------
## TCGA expression data from stad and skcm
exp_data <- readr::read_rds(file.path(TCGA_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types %in% c("SKCM","STAD"))
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header = T)
gene_list$symbol <- as.character(gene_list$symbol)

## clinical info of TCGA samples
survival_data <- readr::read_rds(file.path(TCGA_path,"TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>%
  dplyr::rename("cancer_types" = "type") %>%
  dplyr::rename("survival" = "data")

## clinical expression and response data
clinical_data_for_logistic <- readr::read_rds(file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical","genelist_exp_for_logistic.rds.gz"))

## features filtered by logistic regression
final_feature.5 <- readr::read_tsv(file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical","by_cancer_targets_30test.5","final_feature.5.tsv")) %>%
  tidyr::nest(features,.key="features")

## stepAIC features
final_feature.5.1 <- readr::read_tsv(file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical","by_cancer_targets_30test.5.stepAIC","gastric cancer_ALL30_by-feature-from_anti–PD-1_features.tsv")) %>%
  dplyr::mutate(Cancer.y = "gastric cancer",blockade="anti-PD-1")

final_feature.5.2 <- readr::read_tsv(file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical","by_cancer_targets_30test.5.stepAIC","melanoma_anti-CTLA-4_by-feature-from_anti-CTLA-4_features.tsv")) %>%
  dplyr::mutate(Cancer.y = "melanoma",blockade="anti-CTLA-4")

final_feature.5.3 <- readr::read_tsv(file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical","by_cancer_targets_30test.5.stepAIC","melanoma_anti-CTLA-4_by-feature-from_anti–PD-1_features.tsv")) %>%
  dplyr::mutate(Cancer.y = "melanoma",blockade="anti-PD-1")

final_feature.5.4 <- readr::read_tsv(file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical","by_cancer_targets_30test.5.stepAIC","melanoma_anti–PD-1、CTLA-4_by-feature-from_anti-CTLA-4_features.tsv")) %>%
  dplyr::mutate(Cancer.y = "melanoma",blockade="anti-CTLA-4")

final_feature.5.5 <- readr::read_tsv(file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical","by_cancer_targets_30test.5.stepAIC","melanoma_anti–PD-1、CTLA-4_by-feature-from_anti–PD-1_features.tsv")) %>%
  dplyr::mutate(Cancer.y = "melanoma",blockade="anti-PD-1")

final_feature.5.6 <- readr::read_tsv(file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical","by_cancer_targets_30test.5.stepAIC","melanoma_anti–PD-1_by-feature-from_anti-CTLA-4_features.tsv")) %>%
  dplyr::mutate(Cancer.y = "melanoma",blockade="anti-CTLA-4")

final_feature.5.7 <- readr::read_tsv(file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical","by_cancer_targets_30test.5.stepAIC","melanoma_anti–PD-1_by-feature-from_anti–PD-1_features.tsv")) %>%
  dplyr::mutate(Cancer.y = "melanoma",blockade="anti-PD-1")

final_feature.5.1 %>%
  rbind(final_feature.5.2) %>%
  rbind(final_feature.5.3) %>%
  rbind(final_feature.5.4) %>%
  rbind(final_feature.5.5) %>%
  rbind(final_feature.5.6) %>%
  rbind(final_feature.5.7) %>%
  # dplyr::select(-blockade) %>%
  unique() %>%
  tidyr::nest(features,.key="features")-> final_feature.5.x
## timer TIL data
immunity_path_2 <- file.path(immune_res_path,"ICP_score/2.1.GSVA-ICPs_exp_site_5_feature")
TIMER_immunity_onlyTumor <- readr::read_tsv(file.path(immunity_path_2,"immuneEstimation.txt")) %>%
  dplyr::filter(substr(barcode,14,14) == 0) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  tidyr::gather(-barcode, key = "Cell_type", value = "TIL")

# mutation burden data
mutation_burden_class <- readr::read_rds(file.path(TCGA_path,"classfication_of_26_cancers_by_mutation_burden192.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::rename("cancer_types"="Cancer_Types") %>%
  dplyr::filter(cancer_types %in% c("SKCM","STAD"))

# filter genelist expression ----------------------------------------------
exp_data %>%
  dplyr::mutate(filter_exp = purrr::map(expr,.f = function(.x){
    .x %>% 
      dplyr::filter(symbol %in% gene_list$symbol)
  })) %>%
  dplyr::select(-expr) -> genelist_exp


# combine data for analysis------------------------------------------------

final_feature.5 %>%
  dplyr::mutate(cancer_types = ifelse(Cancer.y == "gastric cancer","STAD","SKCM")) %>%
  dplyr::inner_join(genelist_exp,by="cancer_types") %>% 
  dplyr::inner_join(survival_data,by="cancer_types") %>%
  purrr::pwalk(.f=fn_run)

final_feature.5.x %>%
  dplyr::mutate(cancer_types = ifelse(Cancer.y == "gastric cancer","STAD","SKCM")) %>%
  dplyr::inner_join(genelist_exp,by="cancer_types") %>% 
  dplyr::inner_join(survival_data,by="cancer_types") %>%
  purrr::pwalk(.f=fn_run,blockade="stepAIC_combined")

final_feature.5.x %>%
  dplyr::mutate(cancer_types = ifelse(Cancer.y == "gastric cancer","STAD","SKCM")) %>%
  dplyr::inner_join(genelist_exp,by="cancer_types") %>% 
  dplyr::inner_join(survival_data,by="cancer_types") %>%
  purrr::pwalk(.f=fn_run)
# functions ---------------------------------------------------------------
color_list <- c("1" = "pink1",
                "2" = "skyblue1",
                "3" = "darkseagreen1",
                "4" = "darkgreen",
                "5" = "dodgerblue4",
                "6" = "tan2")

fn_run <- function(Cancer.y,blockade,features,cancer_types,filter_exp,survival){
  survival %>%
    dplyr::rename("barcode" = "bcr_patient_barcode") %>%
    dplyr::select(barcode,PFS,PFS.time) %>%
    dplyr::rename("status" = "PFS", "time" = "PFS.time") -> survival
  # ## get the logistic model
  model <- fn_logistic_model(data_spread,response,features)
  model.summary <- summary(model)
  features.sig.coef <- model.summary$coefficients %>%
    as.data.frame() %>%
    dplyr::as.tbl() %>%
    dplyr::mutate(symbol = rownames(model.summary$coefficients)) %>%
    dplyr::filter(abs(Estimate)>0.2) %>%
    dplyr::select(symbol,Estimate)
  
  ## use features to class samples
  filter_exp %>%
    dplyr::select(-entrez_id) %>%
    dplyr::mutate(symbol = gsub("-",".",symbol)) %>%
    dplyr::inner_join(features.sig.coef,by="symbol") %>%
    tidyr::gather(-symbol,-Estimate,key="barcode",value="exp") %>%
    dplyr::filter(substr(barcode,14,15)=="01") %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(group = ifelse(exp > quantile(exp,0.5),"High", "Low")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(barcode)
    dplyr::mutate(Estimate)
    tidyr::spread(key="symbol",value="exp") %>%
    dplyr::mutate(barcode = substr(barcode,1,12)) %>%
    tidyr::gather(-barcode,key="symbol",value="exp") %>%
    tidyr::spread(key="barcode",value="exp") -> filter_exp.ready
  
  ## cluster analysis
  filter_exp.ready.m <- as.matrix(as.data.frame(filter_exp.ready[,-1]))
  rownames(filter_exp.ready.m) <- filter_exp.ready$symbol
  # results = ConsensusClusterPlus(filter_exp.ready.m,maxK = 5,reps = 100,pItem = 0.8,
  #                                pFeature = 1, 
  #                                title = file.path(res_path, cancer_types),
  #                                distance = "binary", plot = "pdf",
  #                                clusterAlg = "hc",seed = 1262118388.71279)
  
  # Compute the dissimilarity matrix
  dist <- dist(t(filter_exp.ready.m), method = "euclidean")
  hc <- hclust(d = dist, method = "complete")
  
  ## use model to predict TCGA data
  # probabilities <- model %>% predict(filter_exp.ready,type = "response")
  # probabilities %>%
  #   as.data.frame() %>%
  #   dplyr::mutate(barcode = filter_exp.ready$barcode) %>%
  #   dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  #   dplyr::inner_join(survival,by="barcode") %>%
  #   dplyr::rename("probabilities" = ".") %>%
  #   dplyr::mutate(group = ifelse(probabilities> 0.5, "yes", "no")) %>%
  #   dplyr::arrange(group)-> data.for.survival
  ## draw pic 
  for (i in 2:5) {
    C <- i
    group <- cutree(hc, k = C)
    # group <- results[[i]]$consensusClass
    # data.frame(Run = names(group),group = group) %>%
    #   readr::write_tsv(file.path(res_path, pubmed_ID, paste(C,"clusters_group_info.tsv",sep = "_"))) 
    
    data.frame(barcode = names(group),group = group) %>%
      dplyr::as.tbl()  -> group
    
    ###### complex heatmap #####
    survival %>%
      dplyr::inner_join(group, by = "barcode") %>%
      # dplyr::select(Response, group) %>%
      dplyr::mutate(group = as.character(group)) %>%
      as.data.frame() -> sample_info
    # rownames(sample_info) <- sample_info$barcode
    
    # sample_anno <- HeatmapAnnotation(df = sample_info,
    #                                  col = list(group=color_list),
    #                                  width = unit(0.5, "cm"),
    #                                  name = c("CC group"))
    # draw(sample_anno,1:45)
    
    # heatmap
    # png(file.path(res_path,"complex_heatmap",paste(C, pubmed_ID, "GSVAscore", "heatmap.png", sep = "_")),
    #     height = 300, width = 600)
    # Heatmap(filter_exp.ready.m,
    #         show_row_names = T,
    #         show_column_names = FALSE,
    #         cluster_columns = T,
    #         clustering_distance_columns = "euclidean",
    #         clustering_method_columns = "complete",
    #         top_annotation = sample_anno,
    #         show_row_dend = FALSE, # whether show row clusters.
    #         heatmap_legend_param = list(title = c("GSVA score")))
    # dev.off()
    
    #### PCA analysis #####
    # fn_pca_scatter(clinical=clinical,group=group,gsva.dist,cancer_types,pubmed_ID=pubmed_ID,C=C)
    
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
    color <- data.frame(group = c("1","2","3","4","5","6"),
                        color = c("pink1","skyblue1","darkseagreen1","darkgreen","dodgerblue4","tan2"))
    fn_survival(data=sample_info,title=paste(cancer_types,blockade,sep="_"),color=color,
                "group",sur_name=paste(C,cancer_types,"from_features",Cancer.y,blockade,sep="_"),
                xlab = "Time (days)",result_path = file.path(res_path,"survival.cluster"),h=3,w=4,lx=0.8,ly=0.6)
    
    fn_compare(group = sample_info %>% dplyr::mutate(group = as.integer(group)),
               cancer_types = paste(cancer_types,Cancer.y,blockade), data = TIMER_immunity_onlyTumor,
               value = "TIL", facet = "Cell_type", ylab = "TIL", C = C,
               title  = cancer_types, data_type = "TIMER",result_path = file.path(res_path, "TIL"), h = 6, w = 8)
    
    fn_mutation_burden_all(data = sample_info %>% 
                             dplyr::mutate(group = as.integer(group)) %>% 
                             dplyr::inner_join(mutation_burden_class,by="barcode") %>%
                             dplyr::mutate(sm_count = log2(sm_count)),
                           group = "group",value = "sm_count", xlab = "log2(Mutation burden)", 
                           m_a_name = paste(C,cancer_types,Cancer.y,blockade,sep="_"),
                           result_path = file.path(res_path, "MB"), h = 6, w = 8)
    
    apply(filter_exp.ready.m,1,scale)
    tiff(file.path(res_path,paste(pubmed_ID, cancer_types,"GSVAscore", "heatmap.tiff", sep = "_")),
         height = 300, width = 600,compression = c("lzw"))
    pheatmap(res.gsva, #scale = "row",
             clustering_distance_rows = "correlation",
             color = colorRampPalette(c("#00B2EE","white","red"))(50),
             border_color = NA,cutree_rows = 2,
             show_colnames = F, treeheight_col = 30, treeheight_row = 20,
             cutree_cols = 2,main = paste(paste("PMID:",pubmed_ID),cancer_types," GSVA score", sep = ","))
    dev.off()
    # fn_survival.calculate(group, clinical, pubmed_ID, C,cancer_types, result_path = res_path, h = 3, w = 4)
  }
}
## do survival analysis
# color <- tibble::tibble(group = c("no","yes"),
#                         color = c("red","blue"))
# fn_survival(data=data.for.survival,title=paste(cancer_types,"from_features",Cancer.y,blockade,sep="_"),color=color,
#             "group",sur_name=paste(cancer_types,"from_features",Cancer.y,blockade,sep="_"),
#             xlab = "Time (days)",result_path = file.path(res_path,"survival"),h=3,w=4,lx=0.8,ly=0.6)

## construct logistic model
fn_logistic_model <- function(data_spread,response,features){
  data_spread %>%
    dplyr::inner_join(response, by = "Run") %>%
    dplyr::select(-Run) %>%
    dplyr::mutate(Response=as.factor(Response)) -> data.ready
  data.ready <- na.omit(data.ready)
  colnames(data.ready) <- gsub("-",".",colnames(data.ready))
  for(i in 1:nrow(features)){
    if(i == 1){
      formula <- paste("Response", features$features[i], sep = "~")
    } else{
      formula <- paste(formula, features$features[i], sep = "+")
    }
  }
  ## do logistic regression
  model <- glm( as.formula(formula), data = data.ready, family = binomial)
  model
}

## function to draw survival
fn_survival <- function(data,title,color,group,sur_name,xlab,result_path,h,w,lx=0.8,ly=0.6){
  library(survival)
  library(survminer)
  fit <- survfit(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)
  diff <- survdiff(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)
  kmp <- 1 - pchisq(diff$chisq, df = length(levels(as.factor(data$group))) - 1)
  # legend <- data.frame(group=paste("C",sort(unique(data$group)),sep=""),n=fit$n)
  # legend %>%
  #   dplyr::mutate(
  #     label = purrr::map2(
  #       .x = group,
  #       .y = n,
  #       .f = function(.x,.y){
  #         latex2exp::TeX(glue::glue("<<.x>>, n = <<.y>>", .open = "<<", .close = ">>"))
  #       }
  #     )
  #   ) -> legend
  color %>%
    dplyr::inner_join(data,by="group") %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(group = paste(group,", n=",n,sep="")) %>%
    dplyr::select(group,color) %>%
    unique() -> color_paired
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        data = data,
                        surv.median.line = "hv",
                        title = paste(title,", p =", signif(kmp, 2)), # change it when doing diff data
                        xlab = xlab,
                        ylab = 'Probability of survival',
                        # legend.title = "Methyla group:",
                        legend= c(lx,ly),
                        # ggtheme = theme_survminer(),
                        ggtheme = theme(
                          panel.border = element_blank(), panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", 
                                                                                       size = 0.5), 
                          panel.background = element_rect(fill = "white"),
                          legend.key = element_blank(),
                          legend.background = element_blank(),
                          legend.text = element_text(size = 8),
                          axis.text = element_text(size = 12),
                          legend.title = element_blank(),
                          axis.title = element_text(size = 12,color = "black"),
                          text = element_text(color = "black" )
                        )
  )[[1]] +
    scale_color_manual(
      values = color_paired$color,
      labels = color_paired$group
    ) -> p
  ggsave(filename = paste(sur_name,signif(kmp, 2),"png",sep = "."), plot = p, path = result_path,device = "png",height = h,width = w)
  ggsave(filename = paste(sur_name,signif(kmp, 2),"pdf",sep = "."), plot = p, path = result_path,device = "pdf",height = h,width = w)
}

## fn TIL compare
fn_compare <- function(group, cancer_types, data, value, facet, ylab, title, data_type, result_path,C, h, w){
  
  print(paste("Compare",data_type,cancer_types))
  
  group %>%
    dplyr::inner_join(data, by = "barcode") %>%
    dplyr::rename("value" = value) -> plot_ready
  
  assign("group1",facet)
  color_paired <- tibble::tibble(group = c(1:5),
                                 color = c("pink1", "skyblue1", "darkseagreen1", "darkgreen", "dodgerblue4")) %>%
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
    fig_name <- paste("C",C,cancer_types, sep = "_")
    # result_path <- file.path(result_path, data_type)
    ggsave(filename = paste(fig_name,"tiff",sep = "."), path = result_path,device = "tiff",width = w,height = h)
    ggsave(filename = paste(fig_name,"pdf",sep = "."), path = result_path,device = "pdf",width = w,height = h)
  }
}

fn_mutation_burden_all <- function(data,group,value,color,xlab,m_a_name,result_path,w=4,h=3){
  data %>%
    ggpubr::ggboxplot(x = group, y = value,fill="white",alpha = 0,width = 0.1,
                      color = group #add = "jitter",#, palette = "npg"
    ) +
    geom_violin(aes(fill = group),alpha = 0.5) +
    geom_jitter(aes(color = group),alpha = 0.2,width = 0.1,size=0.1) +
    geom_boxplot(fill="white",alpha = 0,width = 0.1) +
    scale_x_discrete(#breaks = c(1:10),
      labels = unique(data$group)
      # expand = c(0.2,0.2,0.2)
    ) +
    # facet_wrap(~ cancer_types, strip.position = "bottom", scales = "free") +
    # scale_color_manual(
    #   values = color
    # )+
    # scale_fill_manual(
    #   values = color
    # ) +
    # ylim(4,12) +
    ylab(xlab) +
    xlab("Group") +
    theme(legend.position = "none",
          # axis.title.x = element_blank(),
          strip.background = element_rect(fill = "white",colour = "white"),
          text = element_text(size = 12, colour = "black"),
          strip.text = element_text(size = 12))+
    ggpubr::stat_compare_means(method = "anova") 
  # ggpubr::stat_compare_means(label.y = 14,paired = TRUE) +
  # ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label = "p.signif")
  ggsave(filename =paste(m_a_name,"png",sep="."), path = result_path,device = "png",width = w,height = h)
  ggsave(filename =paste(m_a_name,"pdf",sep="."), path = result_path,device = "pdf",width = w,height = h)
}
