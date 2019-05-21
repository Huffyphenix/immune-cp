###########################################################
# analysis of the result of combined cluster analysis for scaled exp, methylation clusters, and CNV raw data cluster results
###########################################################

library(SNFtool)
library(ComplexHeatmap)
library(magrittr)

# get best K for cluster results  -----------------------------------------
basic_path <- "/home/huff/project"
source(file.path(basic_path,"github/immune-cp/subtype_analysis/funtions_to_draw_pic.R"))

result_path <- file.path(basic_path,"immune_checkpoint/result_20171025/subtype_analysis/coca_scaled-methyalltag/Get_best_clutser_20")
data_result_path <- file.path(basic_path,"immune_checkpoint/genelist_data")

results <- readr::read_rds(path = file.path(data_result_path, ".recode.combined.ConsensusClusterplus-exp-cnv-methyalltags.rds.gz"))

load(file = file.path(data_result_path, ".rda_genelist_data_survival_cancer.info.rda"))

# get best K for cluster results  -----------------------------------------
source("/project/huff/huff/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

mutation_burden_class <- readr::read_rds(file.path("/project/huff/huff/data/TCGA","classfication_of_26_cancers_by_mutation_burden192.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::rename("cancer_types"="Cancer_Types")
time_status %>%
  dplyr::full_join(mutation_burden_class,by="barcode") %>%
  dplyr::full_join(genelist_methy_mutaion_class.cancer_info,by="barcode") %>%
  dplyr::mutate(cancer_types=ifelse(is.na(cancer_types.x),cancer_types.y,cancer_types.x)) %>%
  dplyr::select(-cancer_types.x,-cancer_types.y) %>%
  dplyr::rename("time"="PFS.time","status"="PFS")-> clinical_data

TIL_data <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/data/immunity","immuneEstimation.txt")) %>%
  dplyr::mutate(barcode = substr(barcode,1,12))
color_20 <- c("#CDC0B0", "#838B8B", "#000000", "#0000FF", "#00008B", "#8A2BE2", "#A52A2A", "#FF4040", "#98F5FF", "#53868B", "#EEAD0E", "#458B00", "#EEA2AD", "#E066FF", "#EE3A8C", "#00FF00", "#FFFF00", "#5CACEE", "#8B6914", "#FF7F24")

C <- 4
W <- results[[C]][['consensusMatrix']]
W %>%
  readr::write_rds(file.path(result_path,paste(C,"consensusMatrix.rds.gz")),compress = "gz")

for (i in 2:20){
  print(paste(i,1,sep = "."))
  data_type <- "COCA"
  C <- i
  W <- results[[C]][['consensusMatrix']]
  group <- results[[C]][['consensusClass']]
  W <- matrix(W,nrow = length(group),dimnames = list(names(group),names(group)))
  data.frame(sample = names(group),group=group) %>%
    readr::write_tsv(file.path(result_path,paste(C,"clusters_group_info.tsv",sep="_")))
  # Figure 1
  # displayClusters(W, group)
  print(paste(i,2,sep = "."))
  
  # Figure 2. survival
  
  group_statistic <- group %>% table()
  less_than_10 <- names(group_statistic[group_statistic<10])
  all_clusters <- names(group_statistic)
  more_than_10 <- setdiff(all_clusters,less_than_10)
  
  data.frame(sample = names(group),group=group) %>%
    dplyr::as.tbl() %>%
    # dplyr::filter(group %in% more_than_10) %>%
    dplyr::rename("barcode"="sample") %>%
    dplyr::left_join(clinical_data,by="barcode") -> group_survival_data
  
  color_list = tibble::tibble(group = more_than_10,color = color_20[1:length(more_than_10)]) %>%
    dplyr::mutate(group = as.integer(group))
  title <- paste(data_type,"Survival for",C,"Clusters",sep=" ")
  sur_name <-  paste(data_type,"Survival_for",C,"Clusters.png",sep="_")
  if (length(more_than_10)>=2) {
    group_survival_data %>%
      dplyr::filter(group %in% more_than_10) %>%
      fn_survival(title,color_list,group="group",sur_name,xlab = "Time (days)",result_path=result_path,h=3,w=4)
  }else{
    print("Too small groups")
  }
  
  print(paste(i,3,sep = "."))
  
  # Figure 3. mutation burden in each cancers
  
  print(paste(i,4,sep = "."))
  if (length(more_than_10)>=2) {
    comp_mat <- combn(more_than_10,2)
    comp_list <- list()
    m_name <-  paste(data_type,"mutation_for",C,"Clusters",sep="_")
    for (col in 1:ncol(comp_mat)) {
      comp_list[[col]] <- comp_mat[,col]  }
    
    group_survival_data %>%
      dplyr::filter(group %in% more_than_10) %>%
      dplyr::mutate(sm_count = ifelse(is.na(sm_count),0,sm_count)) %>%
      dplyr::mutate(sm_count = log2(sm_count)) %>%
      fn_mutation_burden(group = "group",facet="~ cancer_types",value = "sm_count",color = color_list,
                         xlab = "log2(Mutation Burden)",comp_list = comp_list,m_name = m_name,
                         result_path = result_path,w=12)
    
    # Figure 4. mutation burden in all samples
    print(paste(i,5,sep = "."))
    m_a_name <-  paste(data_type,"mutation_all_for",C,"Clusters",sep="_")
    group_survival_data %>%
      dplyr::filter(group %in% more_than_10) %>%
      dplyr::mutate(sm_count = ifelse(is.na(sm_count),0,sm_count)) %>%
      dplyr::mutate(sm_count = log2(sm_count)) %>%
      fn_mutation_burden_all(group = "group",value = "sm_count",color = color_list,
                             xlab = "log2(Mutation Burden)",comp_list = comp_list,m_a_name = m_a_name,
                             result_path = result_path)
    
    # Figure 5. TIL in all samples
    print(paste(i,6,sep = "."))
    T_name <-  paste(data_type,"TIL_for",C,"Clusters",sep="_")
    group_survival_data %>%
      dplyr::filter(group %in% more_than_10) %>%
      dplyr::inner_join(TIL_data,by = "barcode") %>%
      fn_mutation_burden(group = "group",facet="~ cancer_types",value = "CD8_Tcell",color = color_list,
                             xlab = "CD8_Tcell",comp_list = comp_list,m_name = T_name,
                             result_path = result_path,w=12)
    
    # Figure 6. TIL in all samples
    print(paste(i,7,sep = "."))
    T_a_name <-  paste(data_type,"TIL_all_for",C,"Clusters",sep="_")
    group_survival_data %>%
      dplyr::filter(group %in% more_than_10) %>%
      dplyr::inner_join(TIL_data,by = "barcode") %>%
      fn_mutation_burden_all(group = "group",value = "CD8_Tcell",color = color_list,
                             xlab = "CD8_Tcell",comp_list = comp_list,m_a_name = T_a_name,
                             result_path = result_path)
    
  } else{
    print("Too small groups")
  }
  # Figure 5. cancer distrbution in each clusters
  print(paste(i,5,sep = "."))
  group_survival_data %>%
    dplyr::group_by(cancer_types,group) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::select(group,cancer_types,n) %>%
    unique()  %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cancer_types) %>%
    dplyr::mutate(n_a = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(`Sample Composition (%)` = 100*n/n_a) -> cluster_cancers_statistic
  
  cluster_cancers_statistic %>%
    ggplot(aes(x=group,y=cancer_types,fill = `Sample Composition (%)`)) +
    geom_tile(color = "grey") +
    geom_text(aes(label = n)) +
    scale_fill_gradient2(
      limit = c(0, 100),
      breaks = seq(0, 100, 25),
      label = c("0", "25","50","75","100"),
      high = "red",
      na.value = "white"
    ) +
    labs(x = "Cluster", y = "Cancer Types") +
    theme(
      panel.border = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black", size = 0.5), 
      panel.background = element_rect(fill = "white"),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.text = element_text(size = 12, colour = "black"),
      axis.text = element_text(size = 12, colour = "black"),
      legend.title = element_text(angle=90),
      axis.title = element_text(size = 12,color = "black")
    ) +
    guides(fill = guide_colorbar(title.position = "left"))
  ggsave(filename = paste("Sample_composition_for",C,"Clusters.png",sep="_"), path = result_path,device = "png",height = 6,width = 4)
  print(paste(i,6,sep = "."))
  #   return(1)
}
