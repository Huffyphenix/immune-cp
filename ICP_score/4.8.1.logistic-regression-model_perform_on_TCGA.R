##############################################
# apply logistic model get from clinical data on TCGA samples
##############################################
library(magrittr)
library(tidyverse)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(MASS)
library(pROC)

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
TCGA_path <- file.path("S:/study/生存分析/免疫检查点project/liucj_tcga_process_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"ICP_score/4.logistics-model-on-TCGA/by_GSVAScore")

# load data ---------------------------------------------------------------
## TCGA expression data from stad and skcm
exp_data <- readr::read_rds(file.path(TCGA_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types %in% c("SKCM","STAD"))
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header = T)
gene_list$symbol <- as.character(gene_list$symbol)

ICP_family <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/checkpoint/ICP_gene_family.txt"))%>%
  dplyr::inner_join(gene_list, by = "symbol")
ICP_ligand_receptor_pair <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/checkpoint/ICP_gene_ligand_receptor_pairs.txt")) %>%
  dplyr::mutate(pairs = paste("pair",pairs)) %>%
  dplyr::inner_join(gene_list, by = "symbol")
# gene id transfer --------------------------------------------------------
ensembl_gene <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","Homo_sapiens.gene_info.geneid.symbol.ensembl"))

ICPs_exp_in_TI.cancer.site_groups <- 
  readr::read_tsv(file.path(immune_res_path, "ICP_score/2.1.GSVA-ICPs_exp_site_5_feature","ICPs_exp_in_TI.cancer.site_groups.tsv")) %>%
  dplyr::inner_join(gene_list,by = "symbol") %>%
  dplyr::select(symbol, Tissues, `log2FC(I/T)`, TCGA_Cancer, ICP_exp_group, GeneID, type, functionWithImmune)

## clinical info of TCGA samples
survival_data <- readr::read_rds(file.path(TCGA_path,"TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>%
  dplyr::rename("cancer_types" = "type") %>%
  dplyr::rename("survival" = "data")

## clinical expression and response data
clinical_data_for_logistic <- readr::read_rds(file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical","genelist_exp_for_logistic.rds.gz"))

clinical_survival_info <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","RNAseq-sample_info_complete.tsv")) %>%
  dplyr::filter(Library_strategy == "RNA-Seq") %>%
  dplyr::filter(Biopsy_Time == "pre-treatment") %>%
  dplyr::select(Run, Survival_time, Survival_status) %>%
  dplyr::filter(!is.na(Survival_time)) %>%
  dplyr::mutate(Survival_status = ifelse(Survival_status=="Dead",1,0)) %>%
  dplyr::rename("time"="Survival_time","status"="Survival_status")

clinical_TIL <- readr::read_tsv(file.path(immune_res_path,"ICP_score","xCell_TIL_of_clinical_data.tsv")) %>%
  tidyr::gather(-Run, key="Cell_type", value="TIL")

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

# get GSVA score of all possible features of ICP ----------------------------------------------
# calculation of GSVA score -----------------------------------------------
fn_get_gene_list_feature_bysite  <- function(tissue){
  ICPs_exp_in_TI.cancer.site_groups %>%
    dplyr::filter(Tissues %in% tissue) -> .data
  
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

fn_GSVA <- function(tissue, exp){
  print(paste("GSVA",tissue))
  
  #### gene list feature by gene exp site #####
  genelist <- fn_get_gene_list_feature_bysite(tissue)
  #### gene list feature by gene function in immune system #####
  for(fun_class in unique(gene_list$functionWithImmune)){
    genelist[[fun_class]] <- gene_list %>%
      dplyr::filter(functionWithImmune == fun_class) %>%
      .$GeneID
  }
  #### gene list feature by ligand-receptor pairs #####
  for(p in unique(ICP_ligand_receptor_pair$pairs)){
    genelist[[p]] <- ICP_ligand_receptor_pair %>%
      dplyr::filter(pairs == p) %>%
      .$GeneID
  }
  #### gene list feature by gene family #####
  for(f in unique(ICP_family$family)){
    genelist[[f]] <- ICP_family %>%
      dplyr::filter(f == f) %>%
      .$GeneID
  }
  #### do GSVA ####
  if (length(genelist)>0) {
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
      tidyr::gather(-feature, key = "Run",value = "GSVA_score") %>%
      tidyr::spread(key = "feature", value = "GSVA_score") -> gsva.score
  } else{
    gsva.score <- tibble::tibble()
  }
  
  
  # fn_compare_TN(gsva.score, cancer_types, result_path = file.path(res_path,"TN_compare"))
  
  # fn_survival.calculate(gsva.score, cancer_types)
  
  # fn_heatmap(res.gsva)
  gsva.score
}
fn_get_TCGA_sample_class <- function(.x){
  metastic_sample <- colnames(.x)[substr(colnames(.x),14,15) == "06"]
  primary_sample <- colnames(.x)[substr(colnames(.x),14,15) == "01"]
  
  metastic_data <- .x[,c("entrez_id",metastic_sample)]
  primary_data <- .x[,c("entrez_id",primary_sample)]
  
  
  tibble::tibble(Cancer_types = c("metastic","primary"), data_spread = list(metastic_data,primary_data), 
                 sample_n = c(length(metastic_sample),length(primary_sample)))
}

exp_data %>%
  dplyr::mutate(data_spread = purrr::map(expr,.f = fn_get_TCGA_sample_class)) %>%
  dplyr::select(-expr) %>%
  tidyr::unnest() %>%
  dplyr::filter(sample_n > 10) %>%
  dplyr::mutate(tissue = ifelse(cancer_types=="STAD", "stomach", "skin")) %>%
  # head(1) %>%
  dplyr::mutate(GSVA = purrr::map2(tissue, data_spread, fn_GSVA)) %>%
  dplyr::select(-data_spread) -> GSVA.score

save.image(file.path(res_path,".Rdata"))
load(file.path(res_path,".Rdata"))

## features filtered by logistic regression of clinical GSVA score
final_feature.5 <- readr::read_tsv(file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical","by_gsva_score/class_metastic_type-and-gsvascore_with_class_metastic_type/repeat_times_feature","6_final_feature.tsv")) %>%
  tidyr::nest(features,.key="features")

## GSVA data from clinical
clinical_data_for_logistic <- readr::read_rds(file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical","by_gsva_score/class_metastic_type-and-gsvascore_with_class_metastic_type","genelist_exp_for_logistic.rds.gz"))

# combine data for analysis------------------------------------------------
GSVA.score %>%
  dplyr::mutate(Cancer_type=c("metastatic melanoma","advanced melanoma","metastatic gastric cancer")) -> GSVA.score
final_feature.5 %>%
  dplyr::mutate(cancer_types = ifelse(Cancer.y == "gastric cancer","STAD","SKCM")) %>%
  dplyr::inner_join(clinical_data_for_logistic, by = c("Cancer.y","Cancer_type")) %>%
  dplyr::inner_join(GSVA.score,by=c("Cancer_type","cancer_types")) %>% 
  dplyr::inner_join(survival_data,by="cancer_types") %>%
  dplyr::select(-sample_n,-tissue,-sample_group,-Cancer_types) %>%
  purrr::pwalk(.f=fn_run, analysis = "cluster")


# functions ---------------------------------------------------------------

fn_run <- function(Cancer.y,Cancer_type,features,response,data_spread,cancer_types,GSVA,survival,analysis){
  print(paste(Cancer.y,Cancer_type))
  color <- data.frame(group = c("1","2","no","yes","5","6"),
                      color = c("pink1","skyblue1","darkseagreen1","darkgreen","dodgerblue4","tan2"))
  C=2
  survival %>%
    dplyr::rename("barcode" = "bcr_patient_barcode") %>%
    dplyr::select(barcode,PFS,PFS.time) %>%
    dplyr::rename("status" = "PFS", "time" = "PFS.time") -> survival
  # ## get the logistic model
  # model <- fn_logistic_model(data_spread,response,features)
  # model.summary <- summary(model)
  # features.sig.coef <- model.summary$coefficients %>%
  #   as.data.frame() %>%
  #   dplyr::as.tbl() %>%
  #   dplyr::mutate(symbol = rownames(model.summary$coefficients)) %>%
  #   dplyr::filter(abs(Estimate)>0.2) %>%
  #   dplyr::select(symbol,Estimate)
  
  ## use features to class samples
  # GSVA %>%
  #   dplyr::select(-entrez_id) %>%
  #   dplyr::mutate(symbol = gsub("-",".",symbol)) %>%
  #   dplyr::inner_join(features.sig.coef,by="symbol") %>%
  #   tidyr::gather(-symbol,-Estimate,key="barcode",value="exp") %>%
  #   dplyr::filter(substr(barcode,14,15)=="01") %>%
  #   dplyr::group_by(symbol) %>%
  #   dplyr::mutate(group = ifelse(exp > quantile(exp,0.5),"High", "Low")) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::group_by(barcode)
  #   dplyr::mutate(Estimate)
  #   tidyr::spread(key="symbol",value="exp") %>%
  #   dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  #   tidyr::gather(-barcode,key="symbol",value="exp") %>%
  #   tidyr::spread(key="barcode",value="exp") -> filter_exp.ready
  filter_exp.ready <- GSVA
  colnames(filter_exp.ready) <- gsub(" ",".",colnames(filter_exp.ready))
  
  if(analysis == "lm_model"){
    res_path.1 <- file.path(res_path,"by_model")
    dir.create(res_path.1)
    dir.create(file.path(res_path.1,"survival"))
    dir.create(file.path(res_path.1,"TIL"))
    dir.create(file.path(res_path.1,"MB"))
    # get the logistic model
    model <- fn_logistic_model(data_spread,response,features)
    
    # get best cutoff
    # colnames(data_spread) <- gsub(" ",".",colnames(data_spread))
    probabilities <- model %>% predict(type = "response")
    predicted.classes <- ifelse(probabilities >= 0.5, "yes", "no") %>% as.character() %>% as.factor()
    observed.classes <- response$Response
    
    res.roc <- roc(observed.classes, probabilities)
    best_thresholds <- data_frame(
      thresholds = res.roc$thresholds,
      sensitivity = res.roc$sensitivities,
      specificity = res.roc$specificities
    ) %>% 
      dplyr::as.tbl() %>%
      dplyr::mutate(sum = sensitivity + specificity) %>%
      top_n(n = 1, wt = sum) %>%
      .$thresholds
    
    # clinical survival 
    colnames(data_spread) <- gsub(" ",".",colnames(data_spread))
    model %>% predict(data_spread,type = "response") %>%
      as.data.frame() %>%
      dplyr::mutate(Run = data_spread$Run) %>%
      dplyr::inner_join(clinical_survival_info,by="Run") %>%
      dplyr::rename("probabilities" = ".") %>%
      dplyr::mutate(group = ifelse(probabilities> best_thresholds, "yes", "no")) %>%
      dplyr::arrange(group) %>%
      dplyr::mutate(group = as.character(group)) -> clinical_sample_info
    if(nrow(clinical_sample_info)>0){
      fn_survival(data=clinical_sample_info,title=paste(Cancer_type,sep="_"),color=color,
                  "group",sur_name=paste("clinical",cancer_types,"from_features",Cancer_type,Cancer_type,sep="_"),
                  xlab = "Time (days)",result_path = file.path(res_path.1, "clinical/survival"),h=3,w=4,lx=0.8,ly=0.6)
    }
    ## TIL diff between response or not predicted by clinical GSVA score
    model %>% predict(data_spread,type = "response") %>%
      as.data.frame() %>%
      dplyr::mutate(Run = data_spread$Run) %>%
      dplyr::rename("probabilities" = ".") %>%
      dplyr::mutate(group = ifelse(probabilities> best_thresholds, "yes", "no")) %>%
      dplyr::inner_join(clinical_TIL,by="Run") -> clinical_sample_info.TIL
    if(nrow(clinical_sample_info.TIL)>0){
      fn_compare(plot_ready = clinical_sample_info.TIL,cancer_types=Cancer_type, facet = "Cell_type", value="TIL",
                 ylab = "TIL", C = C,
                 title  = cancer_types, data_type = "TIMER",result_path = file.path(res_path.1, "clinical/TIL"), h = 6, w = 8)
    }
    # use model to predict TCGA data
    probabilities <- model %>% predict(filter_exp.ready[,-1],type = "response")
    probabilities %>%
      as.data.frame() %>%
      dplyr::mutate(barcode = filter_exp.ready$Run) %>%
      dplyr::mutate(barcode = substr(barcode,1,12)) %>%
      dplyr::inner_join(survival,by="barcode") %>%
      dplyr::rename("probabilities" = ".") %>%
      dplyr::mutate(group = ifelse(probabilities> best_thresholds, "yes", "no")) %>%
      dplyr::arrange(group) %>%
      dplyr::mutate(group = as.character(group)) -> sample_info
  } else {
    res_path.1 <- file.path(res_path,"by_cluster")
    dir.create(res_path.1)
    dir.create(file.path(res_path.1,"survival"))
    dir.create(file.path(res_path.1,"TIL"))
    dir.create(file.path(res_path.1,"MB"))
    ## cluster analysis
    filter_exp.ready.m <- as.matrix(as.data.frame(filter_exp.ready[,-1]))
    rownames(filter_exp.ready.m) <- filter_exp.ready$Run
    # results = ConsensusClusterPlus(filter_exp.ready.m,maxK = 5,reps = 100,pItem = 0.8,
    #                                pFeature = 1, 
    #                                title = file.path(res_path, cancer_types),
    #                                distance = "binary", plot = "pdf",
    #                                clusterAlg = "hc",seed = 1262118388.71279)
    
    # Compute the dissimilarity matrix
    dist <- dist(filter_exp.ready.m, method = "euclidean")
    hc <- hclust(d = dist, method = "complete")
    C = 2
    group <- cutree(hc, k = C)
    # group <- results[[i]]$consensusClass
    # data.frame(Run = names(group),group = group) %>%
    #   readr::write_tsv(file.path(res_path, pubmed_ID, paste(C,"clusters_group_info.tsv",sep = "_"))) 
    
    data.frame(barcode = names(group),group = group) %>%
      dplyr::as.tbl() %>%
      plyr::mutate(barcode = substr(barcode,1,12)) -> group
    
    
    survival %>%
      dplyr::inner_join(group, by = "barcode") %>%
      # dplyr::select(Response, group) %>%
      dplyr::mutate(group = as.character(group)) %>%
      as.data.frame() -> sample_info
    
    
  }
  
  ## draw pic 
  
  # fn_survival(data=sample_info,title=paste(cancer_types,Cancer_type,sep="_"),color=color,
  #             "group",sur_name=paste(C,cancer_types,"from_features",Cancer_type,Cancer_type,sep="_"),
  #             xlab = "Time (days)",result_path = file.path(res_path.1,"survival"),h=3,w=4,lx=0.8,ly=0.6)

  fn_compare(plot_ready = sample_info %>% dplyr::inner_join(TIMER_immunity_onlyTumor,by="barcode"),
             cancer_types = paste(cancer_types,Cancer_type,Cancer_type),
             value = "TIL", facet = "Cell_type", ylab = "TIL", C = C,
             title  = cancer_types, data_type = "TIMER",result_path = file.path(res_path.1, "TIL"), h = 6, w = 8)

  # fn_mutation_burden_all(data = sample_info %>%
  #                          dplyr::mutate(group = as.integer(group)) %>%
  #                          dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  #                          dplyr::mutate(sm_count = log2(sm_count)),
  #                        group = "group",value = "sm_count", xlab = "log2(Mutation burden)",
  #                        m_a_name = paste(C,cancer_types,Cancer_type,sep="_"),
  #                        result_path = file.path(res_path.1, "MB"), h = 6, w = 8)
}

## construct logistic model
fn_logistic_model <- function(data_spread,response,features){
  data_spread %>%
    dplyr::inner_join(response, by = "Run") %>%
    dplyr::select(-Run) %>%
    dplyr::mutate(Response=as.factor(Response)) -> data.ready
  data.ready[is.na(data.ready)] <- 0
  colnames(data.ready) <- gsub(" ",".",colnames(data.ready))
  for(i in 1:nrow(features)){
    if(i == 1){
      formula <- paste("Response", features$features[i], sep = "~")
    } else{
      formula <- paste(formula, features$features[i], sep = "+")
    }
  }
  ## do logistic regression
  model <- glm( as.formula(formula), data = data.ready, family = binomial)
  
  step.model <- stepAIC(model, direction = "backward",trace = F)
  model <- step.model
}

## function to draw survival
fn_survival <- function(data,title,color,group,sur_name,xlab,result_path,h,w,lx=0.8,ly=0.6){
  if(!dir.exists(result_path)){
    dir.create(result_path,recursive = T)
  }
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
fn_compare <- function(plot_ready, cancer_types, value, facet, ylab, title, data_type, result_path,C, h, w){
  if(!dir.exists(result_path)){
    dir.create(result_path,recursive = T)
  }
  print(paste("Compare",data_type,cancer_types))
  
  # group %>%
  #   dplyr::inner_join(data, by = "barcode") %>%
  #   dplyr::rename("value" = value) -> plot_ready
  
  # assign("group1",facet)
  color_paired <- tibble::tibble(group = c("1","2","no","yes","5","6"),
                                 color = c("pink1","skyblue1","pink1","skyblue1","dodgerblue4","tan2")) %>%
    dplyr::inner_join(plot_ready, by = "group") %>%
    dplyr::select(group, color) %>%
    unique()
  anno_text <- plot_ready %>%
    dplyr::group_by(group,Cell_type) %>%
    dplyr::mutate(n = n(), y = min(TIL) - max(TIL)*0.05) %>%
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
      ggpubr::ggboxplot(x = "group", y = value, fill = "white",alpha = 0,width = 0.1,
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
