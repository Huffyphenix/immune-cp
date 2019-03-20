##########################################
# detatiled analysis of scaled COCA analysis result
##########################################

library(SNFtool)
library(ComplexHeatmap)
library(magrittr)

# data path ---------------------------------------------------------------

data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"
result_path <- file.path(paste("/project/huff/huff/immune_checkpoint/result_20171025/subtype_analysis/coca_scaled/cluster",C,sep="_"))

load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))
source("/project/huff/huff/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

# load COCA result
results <- readr::read_rds(path = file.path(data_result_path, ".nonscale_data.recode.combined.ConsensusClusterplus.rds.gz"))

# best K we got from 4.6.get_best_clusterL_COCA_scaled.R
C = 6

W <- results[[C]][['consensusMatrix']]
group <- results[[C]][['consensusClass']]
W <- matrix(W,nrow = length(group),dimnames = list(names(group),names(group)))

clinical_tcga <- readr::read_rds(file.path("/project/huff/huff/TCGA_survival/data","Pancan.Merge.clinical.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::select(-cancer_types) %>%
  unique() %>%
  dplyr::mutate(OS=as.numeric(OS),Status=as.numeric(Status),Age=as.numeric(Age)) %>%
  dplyr::group_by(barcode) %>%
  dplyr::mutate(OS= max(OS)) %>%
  dplyr::mutate(Status =  max(Status)) %>%
  dplyr::ungroup()

survival_info <- time_status %>%
  dplyr::full_join(clinical_tcga,by="barcode") %>%
  dplyr::select(-Age,-Stage) %>%
  unique() %>%
  dplyr::filter(barcode %in% names(group)) %>% # colnames(W) change to names(group) when doing single data set(expr, cnv and methy data)
  dplyr::mutate(color=ifelse(PFS==1,"red","blue")) %>%
  dplyr::mutate(os_color=ifelse(Status==1,"red","blue")) %>%
  dplyr::mutate(PFS=ifelse(PFS==1,"Dead","Alive")) %>%
  dplyr::mutate(Status = ifelse(Status==1,"Dead","Alive")) 

mutation_info <- mutation_burden_class %>%
  dplyr::right_join(survival_info,by="barcode") %>%
  dplyr::select(barcode,mutation_status) %>%
  dplyr::mutate(color = ifelse(mutation_status=="NA","white","grey")) %>%
  dplyr::mutate(color = ifelse(mutation_status=="high_muation_burden","pink",color)) %>%
  dplyr::mutate(color = ifelse(is.na(color),"white",color)) %>%
  dplyr::mutate(mutation_status = ifelse(is.na(mutation_status),"NA","High")) %>%
  dplyr::mutate(mutation_status = ifelse(color=="grey","Low",mutation_status))

purity_info <- tumor_purity %>%
  dplyr::right_join(survival_info,by="barcode") %>%
  dplyr::select(barcode,purity) %>%
  dplyr::mutate(purity = ifelse(is.na(purity),0,purity))
cancer_color <- readr::read_tsv(file.path("/data/shiny-data/GSCALite","02_pcc.tsv"))
cancer_info <- gene_list_expr.cancer_info %>%
  rbind(gene_list_cnv_gistic.cancer_info) %>%
  rbind(genelist_methy_mutaion_class.cancer_info) %>%
  unique() %>%
  dplyr::select(cancer_types,barcode) %>%
  unique() %>%
  dplyr::filter(barcode %in% names(group)) %>%
  dplyr::inner_join(cancer_color,by = "cancer_types")
cluster_color <- data.frame(Cluster=c(1:C),color = colors()[seq(1,656,32)[1:C+1]])
cluster_info <- group %>%
  as.data.frame() %>%
  dplyr::rename("Cluster" = ".") %>%
  dplyr::mutate(barcode = names(group)) %>%
  dplyr::inner_join(cluster_color,by="Cluster")
M_label=cbind(names(group),group,survival_info$PFS,survival_info$Status,mutation_info$mutation_status,cancer_info$cancer_types,purity_info$purity)
colnames(M_label)=c("barcode","Cluster","Relapse","OS","mutation_burden_class","cancer_types","TumorPurity")

cancer_info$cancer_types %>% unique() %>% length() -> cancer_n
M_label_colors=cbind("Cluster"=cluster_info$color,
                     "Relapse"=survival_info$color,
                     "OS"=survival_info$os_color,
                     "mutation_burden_class"=mutation_info$color,
                     "cancer_types"=cancer_info$color,
                     "tumor_purity" = circlize::colorRamp2(c(0,
                                                             min(purity_info$purity[purity_info$purity>0]),
                                                             max(purity_info$purity[purity_info$purity>0])),
                                                           c("white","grey100","grey0"))(purity_info$purity)
)


# complex heatmap prepare -------------------------------------------------

normalize <- function(X) X/rowSums(X)
ind <- sort(as.vector(group), index.return = TRUE)
ind <- ind$ix

diag(W) <- median(as.vector(W))
W <- normalize(W)
W <- W + t(W)

M_label=data.frame('barcode' = names(group),"Cluster" = group,"Relapse"=survival_info$PFS,"OS"=survival_info$Status,"mutation_burden_class"=mutation_info$mutation_status,"cancer_types"=cancer_info$cancer_types,"TumorPurity"=as.numeric(purity_info$purity))
# rownames(M_label)=colnames(W) # To add if the spectralClustering function
c(W[ind,ind] %>% colnames()) %>%
  as.data.frame() %>%
  dplyr::rename("barcode"=".") %>%
  dplyr::left_join(M_label,by="barcode") %>%
  # dplyr::mutate(Group = ifelse(Cluster %in% c(1,2,3), "Group1","Group2")) %>% # for 6 clusters
  # dplyr::mutate(Group = ifelse(Cluster == 6, "Group3", Group)) %>%
  dplyr::arrange(Cluster,cancer_types) %>%
  dplyr::select(barcode,Cluster,cancer_types,Relapse,OS,mutation_burden_class,TumorPurity)-> M_label
# rownames(M_label)=M_label$barcode# To add if the spectralClustering function

M_label_colors[,"cancer_types"] %>% unique() -> cancer_anno
names(cancer_anno)= c(M_label$cancer_types %>% as.character() %>% unique())

M_label_colors[,"mutation_burden_class"] %>% unique() -> mutaion_anno
names(mutaion_anno)= c(M_label$mutation_burden_class %>% as.character() %>% unique())

M_label_colors[,"Relapse"] %>% unique() -> survival_anno
names(survival_anno)= c(M_label$Relapse %>% as.character() %>%unique())

M_label_colors[,"OS"] %>% unique() -> os_survival_anno
names(os_survival_anno)= c(M_label$OS %>% as.character() %>%unique())

M_label_colors[,"Cluster"] %>% unique() -> cluster_anno
names(cluster_anno)= c(M_label$Cluster %>% unique())

# c("#CD2626", "#CD9B1D", "#00B2EE") -> group_anno
# names(group_anno) <- paste("Group",1:3,sep="")

col_anno <- HeatmapAnnotation(df=M_label[,-1],
                              col = list("Cluster"=cluster_anno,
                                         "Relapse"=survival_anno,
                                         "OS"=os_survival_anno,
                                         "mutation_burden_class"=mutaion_anno,
                                         "cancer_types"=cancer_anno,
                                         "TumorPurity"=circlize::colorRamp2(c(0,
                                                                              min(M_label$TumorPurity[M_label$TumorPurity>0]),
                                                                              max(M_label$TumorPurity[M_label$TumorPurity>0])),
                                                                            c("white","grey100","grey0"))),
                              gap = unit(0.5, "mm"),
                              width = unit(0.5, "cm"))
draw(col_anno,1:8255)

library(circlize)
he = Heatmap(W[ind,ind],
             col = RColorBrewer::brewer.pal(9,"YlOrRd"),
             show_row_names = FALSE, 
             show_column_names = FALSE,
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             show_row_dend = FALSE, # whether show row clusters.
             top_annotation = col_anno,show_heatmap_legend = F
             # heatmap_legend_param = list(title = c("Scaled Exp."))
)
pdf(file.path(result_path,paste("clusters",C,"_group_heatmap.pdf",sep = "_")),width = 6,height = 6)
png(filename = file.path(result_path,paste("clusters",C,"_group_heatmap.png",sep = "_")),
    width = 600, height = 600, units = "px", pointsize = 12,
    bg = "white")
he
dev.off()


# sample distribution in cancers ------------------------------------------
data.frame(barcode = names(group),group=group) %>%
  dplyr::as.tbl() %>%
  dplyr::left_join(cancer_info,by="barcode") -> group_cancer_data
group_cancer_data %>%
  dplyr::group_by(cancer_types,group) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::select(group,cancer_types,n) %>%
  unique()  %>%
  dplyr::ungroup() %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(n_a = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(`Sample Composition (%)` = 100*n/n_a) %>%
  dplyr::mutate(tcga = "tcga") %>%
  dplyr::group_by(tcga) %>%
  dplyr::mutate(all_n = sum(n_a)) %>%
  dplyr::mutate(all_ratio = 100*n_a/all_n) -> cluster_cancers_statistic

cluster_cancers_statistic %>%
  dplyr::mutate(group=as.character(group)) %>%
  ggplot(aes(x=group,y=cancer_types,fill = `Sample Composition (%)`)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = n)) +
  scale_x_discrete(limit = as.character(c(1:6))) +
  scale_fill_gradient2(
    limit = c(0, 100),
    breaks = seq(0, 100, 25),
    label = c("0", "25","50","75","100"),
    high = "red",
    na.value = "white"
  ) +
  labs(x = "Cluster", y = "Cancer Types") +
  theme(
    # panel.border = element_blank(), 
    # panel.grid.major = element_blank(), 
    # panel.grid.minor = element_blank(), 
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
ggsave(filename =paste("Sample_composition_for",C,"Clusters-heatmap.png",sep="_"), path = result_path,device = "png",height = 6,width = 4)
ggsave(filename =paste("Sample_composition_for",C,"Clusters-heatmap.pdf",sep="_"), path = result_path,device = "pdf",height = 6,width = 4)

cluster_cancers_statistic %>%
  dplyr::mutate(group=as.character(group)) %>%
  dplyr::mutate(width = n_a/max(n_a)) %>%
  dplyr::mutate(width = ifelse(width == min(width),0.05333333,width)) %>%
  ggplot(aes(x=cancer_types,y=`Sample Composition (%)`)) +
  geom_bar(aes(fill=group,width=width),stat="identity",color= "grey")+
  scale_fill_manual(
    name = "Clusters",
    breaks = as.character(c(1:20)),
    values = c("#000000", "#0000FF","#8A2BE2", "#A52A2A", "#98F5FF", "#53868B", "#EEAD0E", "#458B00", "#EEA2AD", "#E066FF", "#EE3A8C", "#00FF00", "#FFFF00", "#5CACEE", "#8B6914", "#FF7F24")
  ) +
  scale_x_discrete(position = "top") +
  theme(
    axis.text.x.top = element_text(angle = 90, vjust=0.5, hjust=0),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(size = 10, colour = "black"),
    axis.text = element_text(size = 10,color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 8, colour = "black"),
    # axis.line = element_blank(),
    # axis.title = element_blank(),
    # axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    axis.line.x = element_blank()
  )
ggsave(filename =paste("Sample_composition_for",C,"Clusters-stacked.png",sep="_"), path = result_path,device = "png",height = 4,width = 12)
ggsave(filename =paste("Sample_composition_for",C,"Clusters-stacked.pdf",sep="_"), path = result_path,device = "pdf",height = 4,width = 12)

# survival ----------------------------------------------------------------
# survival function
fn_survival <- function(data,title,color,sur_name,result_path,h,w,lx=0.8,ly=0.6){
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
  tibble::tibble(group = group, color = color) %>%
    dplyr::inner_join(data,by="group") %>%
    dplyr::mutate(group = paste("C",group,sep="")) %>%
    dplyr::select(group,color) %>%
    unique() -> color_paired
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        data = data,
                        surv.median.line = "hv",
                        title = paste(title,", p =", signif(kmp, 2)), # change it when doing diff data
                        xlab = "Survival in days",
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
                          legend.text = element_text(size = 8, colour = "black"),
                          axis.text = element_text(size = 12, colour = "black"),
                          legend.title = element_blank(),
                          axis.title = element_text(size = 12,color = "black"),
                          legend.key.size = unit(0.2, "cm")
                        )
  ) +
    scale_color_manual(
      values = color_paired$color,
      labels = color_paired$group
    )
  ggsave(filename = paste(sur_name,signif(kmp,2),"png",sep="."), path = result_path,device = "png",height = h,width = w)
  ggsave(filename = paste(sur_name,signif(kmp,2),"pdf",sep="."), path = result_path,device = "pdf",height = h,width = w)
}

res_path <- file.path(result_path)
M_label %>%
  dplyr::as.tbl() -> C6_survival

color_list <- c("#000000", "#0000FF","#8A2BE2", "#A52A2A", "#98F5FF", "#53868B")
group <- c(1,2,3,4,5,6)


# cluster with groups ------ 5 years PFS
C6_survival %>%
  dplyr::rename("group" = "Cluster") %>%
  dplyr::inner_join(time_status,by="barcode") %>%
  dplyr::rename("time"="PFS.time","status"="PFS") -> C6_survival.PFS

sur_name <- paste("PFS_Survival-5years_for",C,"Groups",sep="_")
C6_survival.PFS %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("6 Clutsers, PFS",color_list,sur_name,res_path,3,4,0.9,0.9)

sur_name <- paste("PFS_Survival-all_years_for",C,"Groups",sep="_")
C6_survival.PFS %>%
  fn_survival("6 Clutsers, PFS",color_list,sur_name,res_path,3,4,0.9,0.9)

# cluster with groups for each cancers
C6_survival.PFS %>%
  dplyr::group_by(cancer_types,group) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,group,n) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::arrange(cancer_types) %>%
  dplyr::select(cancer_types) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::filter(Freq >=2) %>%
  .$. %>% as.character() -> cancers_do_survival

for(cancer in cancers_do_survival){
  sur_name <- paste("Survival-5years","PFS",cancer,3,"Groups",sep="_")
  
  C6_survival.PFS %>%
    dplyr::filter(time<=1825) %>%
    dplyr::filter(cancer_types == cancer) %>%
    fn_survival(paste(cancer," PFS",sep=","),color_list,sur_name,file.path(res_path,"PFS"),3,4)
}

for(cancer in cancers_do_survival){
  sur_name <- paste("Survival-all_years","PFS",cancer,3,"Groups",sep="_")
  
  C6_survival.PFS %>%
    dplyr::filter(cancer_types == cancer) %>%
    fn_survival(paste(cancer," PFS",sep=","),color_list,sur_name,file.path(res_path,"PFS"),3,4)
}

# cluster with groups ------ 5 years OS
C6_survival %>%
  dplyr::rename("group" = "Cluster") %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::rename("time"="OS.y","status"="Status") %>%
  dplyr::select(barcode,group,cancer_types,time,status) %>%
  unique() -> C6_survival.OS

sur_name <- paste("OS_Survival-all_years_for",C,"Groups",sep="_")
C6_survival.OS %>%
  fn_survival("6 Clutsers, OS",color_list,sur_name,res_path,3,4,0.9,0.9)

sur_name <- paste("OS_Survival-5years_for",C,"Groups",sep="_")
C6_survival.OS %>%
  dplyr::filter(time<=1825) %>%
  fn_survival("6 Clutsers, OS",color_list,sur_name,res_path,3,4,0.9,0.9)

# cluster with groups for each cancers
C6_survival.OS %>%
  dplyr::group_by(cancer_types,group) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,group,n) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::arrange(cancer_types) %>%
  dplyr::select(cancer_types) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::filter(Freq >=2) %>%
  .$. %>% as.character() -> cancers_do_survival

for(cancer in cancers_do_survival){
  sur_name <- paste("Survival-5years","OS",cancer,3,"Groups",sep="_")
  
  C6_survival.OS %>%
    dplyr::filter(time<=1825) %>%
    dplyr::filter(cancer_types == cancer) %>%
    fn_survival(paste(cancer," OS",sep=","),color_list,sur_name,file.path(res_path,"OS"),3,4)
}

for(cancer in cancers_do_survival){
  sur_name <- paste("Survival-all_years","OS",cancer,3,"Groups",sep="_")
  
  C6_survival.OS %>%
    dplyr::filter(cancer_types == cancer) %>%
    fn_survival(paste(cancer," OS",sep=","),color_list,sur_name,file.path(res_path,"OS"),3,4)
}

# mutation burden ---------------------------------------------------------

C6_survival.PFS %>%
  # dplyr::rename("cluster" = "group","barcode"="sample") %>%
  dplyr::left_join(mutation_burden_class,by="barcode") %>%
  dplyr::mutate(sm_count = ifelse(is.na(sm_count),0,sm_count)) %>%
  dplyr::rename("cancer_types"="cancer_types.x") -> group_cluster_mutation
comp_mat <- combn(length(unique(group_cluster_mutation$group)),2)
comp_list <- list()

m_name <-  paste("Combined_mutation_for",C,"Groups",sep="_")
for (col in 1:ncol(comp_mat)) {
  comp_list[[col]] <- comp_mat[,col]
}
# color_list <- c("#CD2626", "#CD9B1D", "#00B2EE")
group_cluster_mutation %>%
  dplyr::mutate(sm_count = log2(sm_count)) %>%
  fn_mutation_burden("group",facet="~ cancer_types","sm_count",color_list,"log2 (MutationBurden)",comp_list,m_name,res_path,6,6)

# Figure 4
m_a_name <-  paste("Combined_all_mutation_for",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::mutate(sm_count = log2(sm_count)) %>%
  fn_mutation_burden_all("group","sm_count",color_list,"log2 (MutationBurden)",comp_list,m_a_name,res_path)


# tumor purity ------------------------------------------------------------
m_name <-  paste("Combined_TumorPurity_for",C,"Groups",sep="_")
group_cluster_mutation %>%
  fn_mutation_burden("group",facet="~ cancer_types","TumorPurity",color_list,"Tumor Purity",comp_list,m_name,res_path,8,8)

m_a_name <-  paste("Combined_all_TumorPurity_for",C,"Groups",sep="_")
group_cluster_mutation %>%
  fn_mutation_burden_all("group","TumorPurity",color_list,"Tumor Purity",comp_list,m_a_name,res_path)


# immune infiltration  -----------------------------------------------------
# miao's data (TCAP)
TCGA_infiltration_data <- readr::read_rds(file.path("/project/huff/huff/data/TCGA/immune_infiltration/miao_TCAP_prediction_for_all_samples/All_TCGA_sample_TIL.rds.gz")) %>%
  dplyr::select(barcode,InfiltrationScore)

m_name <-  paste("Combined_ImmuneInfiltration_for",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::left_join(TCGA_infiltration_data,by="barcode") %>%
  dplyr::filter(!is.na(InfiltrationScore)) %>%
  dplyr::mutate(InfiltrationScore = as.numeric(as.character(InfiltrationScore))) %>%
  fn_mutation_burden("group",facet="~ cancer_types","InfiltrationScore",color_list,"Immune Infiltration Score",comp_list,m_name,res_path,h=12)

m_a_name <-  paste("Combined_all_ImmuneInfiltration_for",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(InfiltrationScore = as.numeric(as.character(InfiltrationScore))) %>%
  fn_mutation_burden_all("group","InfiltrationScore",color_list,"Immune Infiltration Score",comp_list,m_a_name,res_path)

all_TCGA_infiltration_data <- readr::read_rds(file.path("/project/huff/huff/data/TCGA/immune_infiltration/miao_TCAP_prediction_for_all_samples/All_TCGA_sample_TIL.rds.gz")) 

m_name <-  paste("Combined_TIL_for","CD4_naive",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(CD4_naive = as.numeric(as.character(CD4_naive))) %>%
  fn_mutation_burden_all("group","CD4_naive",color_list,"CD4 naive",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","CD8_naive",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(CD8_naive = as.numeric(as.character(CD8_naive))) %>%
  fn_mutation_burden_all("group","CD8_naive",color_list,"CD8 naive",comp_list,m_name,res_path)

  m_name <-  paste("Combined_TIL_for","Cytotoxic",c,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Cytotoxic = as.numeric(as.character(Cytotoxic))) %>%
  fn_mutation_burden_all("group","Cytotoxic",color_list,"Cytotoxic",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","Exhausted",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Exhausted = as.numeric(as.character(Exhausted))) %>%
  fn_mutation_burden_all("group","Exhausted",color_list,"Exhausted",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","Tr1",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Tr1 = as.numeric(as.character(Tr1))) %>%
  fn_mutation_burden_all("group","Tr1",color_list,"Tr1",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","nTreg",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(nTreg = as.numeric(as.character(nTreg))) %>%
  fn_mutation_burden_all("group","nTreg",color_list,"nTreg",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","iTreg",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(iTreg = as.numeric(as.character(iTreg))) %>%
  fn_mutation_burden_all("group","iTreg",color_list,"iTreg",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","Th1",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Th1 = as.numeric(as.character(Th1))) %>%
  fn_mutation_burden_all("group","Th1",color_list,"Th1",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","Th2",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Th2 = as.numeric(as.character(Th2))) %>%
  fn_mutation_burden_all("group","Th2",color_list,"Th2",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","Th17",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Th17 = as.numeric(as.character(Th17))) %>%
  fn_mutation_burden_all("group","Th17",color_list,"Th17",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","Tfh",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Tfh = as.numeric(as.character(Tfh))) %>%
  fn_mutation_burden_all("group","Tfh",color_list,"Tfh",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","Central_memory",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Central_memory = as.numeric(as.character(Central_memory))) %>%
  fn_mutation_burden_all("group","Central_memory",color_list,"Central_memory",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","Effector_memory",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Effector_memory = as.numeric(as.character(Effector_memory))) %>%
  fn_mutation_burden_all("group","Effector_memory",color_list,"Effector_memory",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","NKT",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(NKT = as.numeric(as.character(NKT))) %>%
  fn_mutation_burden_all("group","NKT",color_list,"NKT",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","MAIT",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(MAIT = as.numeric(as.character(MAIT))) %>%
  fn_mutation_burden_all("group","MAIT",color_list,"MAIT",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","DC",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(DC = as.numeric(as.character(DC))) %>%
  fn_mutation_burden_all("group","DC",color_list,"DC",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","Bcell",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Bcell = as.numeric(as.character(Bcell))) %>%
  fn_mutation_burden_all("group","Bcell",color_list,"Bcell",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","Monocyte",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Monocyte = as.numeric(as.character(Monocyte))) %>%
  fn_mutation_burden_all("group","Monocyte",color_list,"Monocyte",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","Macrophage",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Macrophage = as.numeric(as.character(Macrophage))) %>%
  fn_mutation_burden_all("group","Macrophage",color_list,"Macrophage",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","NK",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(NK = as.numeric(as.character(NK))) %>%
  fn_mutation_burden_all("group","NK",color_list,"NK",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","Neutrophil",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Neutrophil = as.numeric(as.character(Neutrophil))) %>%
  fn_mutation_burden_all("group","Neutrophil",color_list,"Neutrophil",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","Gamma_delta",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(Gamma_delta = as.numeric(as.character(Gamma_delta))) %>%
  fn_mutation_burden_all("group","Gamma_delta",color_list,"Gamma_delta",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","CD4_T",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(CD4_T = as.numeric(as.character(CD4_T))) %>%
  fn_mutation_burden_all("group","CD4_T",color_list,"CD4_T",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","CD4_T_cancers",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(CD4_T = as.numeric(as.character(CD4_T))) %>%
  fn_mutation_burden("group",facet = "~ cancer_types","CD4_T",color_list,"CD4_T",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","CD8_T",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(CD8_T = as.numeric(as.character(CD8_T))) %>%
  fn_mutation_burden_all("group","CD8_T",color_list,"CD8_T",comp_list,m_name,res_path)

m_name <-  paste("Combined_TIL_for","all",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(all_TCGA_infiltration_data,by="barcode") %>%
  dplyr::select(-Relapse,-OS,-mutation_burden_class,-TumorPurity ,-status,-time,-cancer_types.y,-sm_count,-mutation_status,-cancers) %>%
  tidyr::gather(-group,-cancer_types,-barcode,key="TIL_class",value="score") %>%
  dplyr::mutate(score = as.numeric(score)) %>%
  fn_mutation_burden("group","score",facet="~ TIL_class",color_list,"score",comp_list,m_name,res_path)

# TIMER data
TIMER_immunity <- readr::read_tsv(file.path("/project/huff/huff/immune_checkpoint/data/immunity","immuneEstimation.txt")) %>%
  dplyr::mutate(TIL = B_cell+CD4_Tcell+CD8_Tcell+Neutrophil+Macrophage+Dendritic) %>%
  dplyr::mutate(group =  substr(barcode,14,14)) %>%
  dplyr::filter(group == "0") %>%
  dplyr::select(-group) %>%
  dplyr::mutate(barcode = substr(barcode,1,12))

m_name <-  paste("Combined_TIMER_for",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::left_join(TIMER_immunity,by="barcode") %>%
  dplyr::filter(!is.na(TIL)) %>%
  fn_mutation_burden("group",facet="~ cancer_types","TIL",color_list,"TIMER Score",comp_list,m_name,res_path,h=6)

m_a_name <-  paste("Combined_TIMER-B_cell-for",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  dplyr::filter(!is.na(TIL)) %>%
  fn_mutation_burden_all("group","B_cell",color_list,"TIMER B_cell Score",comp_list,m_a_name,res_path)

m_a_name <-  paste("Combined_TIMER-CD4_Tcell-for",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  dplyr::filter(!is.na(TIL)) %>%
  fn_mutation_burden_all("group","CD4_Tcell",color_list,"TIMER CD4_Tcell Score",comp_list,m_a_name,res_path)

m_a_name <-  paste("Combined_TIMER-CD8_Tcell-for",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  dplyr::filter(!is.na(TIL)) %>%
  fn_mutation_burden_all("group","CD8_Tcell",color_list,"TIMER CD8_Tcell Score",comp_list,m_a_name,res_path)

m_a_name <-  paste("Combined_TIMER-Neutrophil-for",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  dplyr::filter(!is.na(TIL)) %>%
  fn_mutation_burden_all("group","Neutrophil",color_list,"TIMER Neutrophil Score",comp_list,m_a_name,res_path)

m_a_name <-  paste("Combined_TIMER-Macrophage-for",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  dplyr::filter(!is.na(TIL)) %>%
  fn_mutation_burden_all("group","Macrophage",color_list,"TIMER Macrophage Score",comp_list,m_a_name,res_path)

m_a_name <-  paste("Combined_TIMER-Dendritic-for",C,"Groups",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  dplyr::filter(!is.na(TIL)) %>%
  fn_mutation_burden_all("group","Dendritic",color_list,"TIMER Dendritic Score",comp_list,m_a_name,res_path)

# 10 cancer related pathways score ----------------------------------------

TCGA_rppa <- readr::read_rds(file.path("/data/TCGA/TCGA_data","pancan32_rppa_score.rds.gz")) %>%
  tidyr::unnest() 

group_cluster_mutation %>%
  dplyr::inner_join(TCGA_rppa,by="barcode") %>%
  dplyr::select(barcode,group,cancer_types.x,group,pathway,score) -> pathway_score_for_groups

fn_ttest_pathway <- function(p_c,.data,i=1){
  print(p_c)
  # print(pathway)
  table(.data$group) %>%
    as.data.frame() %>%
    dplyr::as.tbl() %>%
    dplyr::filter(Freq > 3) %>%
    .$Var1 -> valid_g
  length(valid_g) -> n
  if(n<=1){
    return(tibble::tibble())
  } else{
    if(n==2){
      t.test(score ~ group,data = .data %>% dplyr::filter(group %in% valid_g)) %>%
        broom::tidy() %>%
        dplyr::select(p.value)
    }else{
      oneway.test(score ~ group, data = .data %>% dplyr::filter(group %in% valid_g)) %>%
        broom::tidy() %>%
        dplyr::select(p.value)
    }
  }
}
pathway_score_for_groups %>%
  dplyr::as.tbl() %>%
  tidyr::nest(-c("pathway","cancer_types.x")) -> pathway_score_for_groups.cancer.sep
pathway_score_for_groups.cancer.sep %>% 
  dplyr::mutate(p_c = paste(cancer_types.x,pathway,sep="_")) %>%
  dplyr::group_by(p_c) %>%
  dplyr::mutate(ttest = purrr::map2(p_c,data,fn_ttest_pathway)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() -> pathway_score_for_groups.diff.pva
pathway_score_for_groups.diff.pva %>%
  tidyr::unnest() %>%
  dplyr::filter(p.value <=0.05) -> pathway_score_for_groups.diff.pva.sig

pathway_score_for_groups %>%
  dplyr::as.tbl() %>%
  tidyr::nest(-c("pathway")) -> pathway_score_for_groups.pathway.sep
pathway_score_for_groups.pathway.sep %>% 
  dplyr::group_by(pathway) %>%
  dplyr::mutate(ttest = purrr::map2(pathway,data,fn_ttest_pathway)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() -> pathway_score_for_all.diff.pva
pathway_score_for_all.diff.pva %>%
  tidyr::unnest() %>%
  dplyr::filter(p.value <=0.05) -> pathway_score_for_all.diff.pva.sig


pathway_score_for_groups %>%
  dplyr::filter(pathway %in% pathway_score_for_all.diff.pva.sig$pathway) %>%
  ggpubr::ggboxplot(x = "group", y = "score",
                    color = "white" #add = "jitter",#, palette = "npg"
  ) +
  geom_jitter(aes(color = group),size=0.2,alpha=0.2) +
  geom_violin(fill="white",alpha = 0) +
  geom_boxplot(fill="white",alpha = 0,width=0.1) +
  # geom_point(aes(x=as.numeric(group)+b,y=sm_count,color=group),alpha = 0.5) +
  scale_x_discrete(#breaks = c(1:3),
    labels = c(1:6)
    # expand = c(0.2,0.2,0.2)
  ) +
  # facet_grid(cancer_types.x~ pathway, scales = "free") +
  facet_wrap(~ pathway, strip.position = "bottom", scales = "free") +
  scale_color_manual(
    values = color_list
  )+
  # ylim(4,12) +
  ylab("Pathway Score") +
  xlab("Group") +
  theme(strip.background = element_rect(fill = "white",colour = "white"),
        legend.position = "none",
        text = element_text(size = 12,colour = "black"),
        strip.text = element_text(size = 10,colour = "black")) +
  ggpubr::stat_compare_means(label.y = 20,paired = TRUE)
  # ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label = "p.signif")

# p_name <- paste("Combined_PathwayScore_for",C,"Groups",sep="_")
# ggsave(filename = paste(p_name,"png",sep="."), path = res_path,device = "png",width = 8,height = 6)
# ggsave(filename = paste(p_name,"pdf",sep="."), path = res_path,device = "pdf",width = 8,height = 6)

p_a_name <- paste("Combined_all_PathwayScore_for",C,"Groups.pdf",sep="_")
ggsave(filename = paste(p_a_name,"png",sep="."), path = res_path,device = "png",width = 8,height = 6)
ggsave(filename = paste(p_a_name,"pdf",sep="."), path = res_path,device = "pdf",width = 8,height = 6)


# virus related pathways score ----------------------------------------
TCGA_virus <- readr::read_rds(file.path("/project/huff/huff/data/TCGA/virus_level","TCGA_virus_level.rds.gz")) %>%
  dplyr::select(-Provenance)

group_cluster_mutation %>%
  dplyr::select(barcode,cluster,cancer_types,group) %>%
  dplyr::left_join(TCGA_virus,by="barcode") %>%
  tidyr::gather(-c("barcode","cluster","cancer_types","group"),key=virus,value=virus_levels) %>% 
  tidyr::drop_na(virus_levels) %>%
  # dplyr::mutate(virus_levels=ifelse(is.na(virus_levels),0,virus_levels)) %>%
  dplyr::mutate(virus_levels=as.numeric(virus_levels)) -> virus_for_groups

virus_for_groups %>%
  dplyr::group_by(group,virus) %>%
  dplyr::mutate(n=n(),y=max(virus_levels)) %>%
  dplyr::mutate(groupx = substr(as.character(group),6,6)) %>%
  dplyr::select(groupx,virus,n,y) %>%
  unique() -> virus.text

virus_for_groups %>%
  dplyr::filter(virus %in% c("EBV","HPV")) %>%
  dplyr::mutate(virus_levels = as.numeric(as.character(virus_levels))) %>%
  ggpubr::ggboxplot(x = "group", y = "virus_levels",
                    color = "white" #add = "jitter",#, palette = "npg"
  ) +
  geom_jitter(aes(color = group)) +
  # geom_text(aes(x=group,y=y,label = n),data = virus.text) +
  # geom_point(aes(x=as.numeric(group)+b,y=sm_count,color=group),alpha = 0.5) +
  scale_x_discrete(#breaks = c(1:3),
    labels = c(1:3)
    # expand = c(0.2,0.2,0.2)
  ) +
  # facet_grid(cancer_types~ virus, scales = "free") +
  facet_wrap(~ virus, strip.position = "bottom", scales = "free") +
  scale_color_manual(
    values = color_list
  )+
  # ylim(4,12) +
  ylab("Virus Level") +
  xlab("Group") +
  theme(strip.background = element_rect(fill = "white",colour = "white"),
        legend.position = "none",
        text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  # ggpubr::stat_compare_means(label.y = 14,paired = TRUE) +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label = "p.signif")

v_name <- paste("Combined_all_virus_levels_for",3,"Groups",sep="_")
ggsave(filename =paste(v_name,"pdf",sep = "."), path = res_path,device = "pdf",width = 4,height = 3)
ggsave(filename =paste(v_name,"png",sep = "."), path = res_path,device = "png",width = 4,height = 3)

virus_for_groups %>% dplyr::filter(virus_levels>0) %>% .$cancer_types %>% table() %>% as.data.frame() %>% dplyr::as.tbl() %>%
  dplyr::filter(Freq > 5) %>% .$. %>% as.character()-> index

virus_for_groups %>%
  dplyr::filter(virus %in% c("EBV","HPV")) %>%
  dplyr::filter(cancer_types %in% index) %>%
  ggpubr::ggboxplot(x = "group", y = "virus_levels",
                    color = "group" #add = "jitter",#, palette = "npg"
  ) +
  # geom_text(aes(x=group,y=y,label = n),data = virus.text) +
  # geom_point(aes(x=as.numeric(group)+b,y=sm_count,color=group),alpha = 0.5) +
  scale_x_discrete(#breaks = c(1:3),
    labels = c(1:3)
    # expand = c(0.2,0.2,0.2)
  ) +
  facet_grid(cancer_types~ virus, scales = "free") +
  # facet_wrap(~ virus, strip.position = "bottom", scales = "free") +
  scale_color_manual(
    values = color_list
  )+
  # ylim(4,12) +
  ylab("Virus Level") +
  xlab("Group") +
  theme(strip.background = element_rect(fill = "white",colour = "white"),
        legend.position = "none",
        text = element_text(size = 5),
        strip.text = element_text(size = 8)) +
  # ggpubr::stat_compare_means(label.y = 14,paired = TRUE) +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label = "p.signif")
v_a_name <- paste("Combined_virus_levels_for",3,"Groups.pdf",sep="_")
ggsave(filename =v_a_name, path = res_path,device = "pdf",width = 6,height = 12)

save.image(file = file.path(res_path, ".rda_combined_clusters_into_group_analysis.rda"))
load(file = file.path(res_path, ".rda_combined_clusters_into_group_analysis.rda"))
