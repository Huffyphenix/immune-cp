library(SNFtool)
library(ComplexHeatmap)

# data path ---------------------------------------------------------------

result_path <- file.path("/project/huff/huff/immune_checkpoint/result_20171025/subtype_analysis/combine/Get_best_clutser_20")
load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))

# cluster K ---------------------------------------------------------------

C <- 6
print(C)
group <- readr::read_tsv(file.path(result_path,paste(C,"clusters_group_info.tsv",sep="_")))

results <- readr::read_rds(file.path("/project/huff/huff/immune_checkpoint/genelist_data",".rds_PanCan28_combine-expr-cnv-methy_snf_20.rds.gz"))
W <- results$distanceMatrix

displayClusters(W, group$group)

survival_info <- time_status %>%
  dplyr::filter(barcode %in% colnames(W)) %>% # colnames(W) change to names(group) when doing single data set(expr, cnv and methy data)
  dplyr::mutate(color=ifelse(PFS==1,"red","blue")) %>%
  dplyr::mutate(PFS=ifelse(PFS==1,"Dead","Alive")) 

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
cancer_info <- cnv_merge_snv_data %>%
  dplyr::select(cancer_types,barcode) %>%
  unique() %>%
  dplyr::filter(barcode %in% colnames(W)) %>%
  dplyr::inner_join(cancer_color,by = "cancer_types")
M_label=cbind(group,survival_info$PFS,mutation_info$mutation_status,cancer_info$cancer_types,purity_info$purity)
colnames(M_label)=c("barcode","spectralClustering","PFS.status","mutation_burden_class","cancer_types","TumorPurity")
rownames(M_label)=colnames(W)
cancer_info$cancer_types %>% unique() %>% length() -> cancer_n
M_label_colors=cbind("spectralClustering"=getColorsForGroups(M_label[,"spectralClustering"],
                                                             colors=rainbow(C)),
                     "PFS.status"=survival_info$color,
                     "mutation_burden_class"=mutation_info$color,
                     "cancer_types"=cancer_info$color,
                     "tumor_purity" = circlize::colorRamp2(c(0,
                                                             min(purity_info$purity[purity_info$purity>0]),
                                                             max(purity_info$purity[purity_info$purity>0])),
                                                           c("white","grey100","grey0"))(purity_info$purity)
)


# complex heatmap prepare -------------------------------------------------

normalize <- function(X) X/rowSums(X)
ind <- sort(as.vector(group$group), index.return = TRUE)
ind <- ind$ix

diag(W) <- median(as.vector(W))
W <- normalize(W)
W <- W + t(W)

M_label=data.frame(group,"PFS.status"=survival_info$PFS,"mutation_burden_class"=mutation_info$mutation_status,"cancer_types"=cancer_info$cancer_types,"TumorPurity"=as.numeric(purity_info$purity))
rownames(M_label)=colnames(W) # To add if the spectralClustering function

M_label_colors[,"cancer_types"] %>% unique() -> cancer_anno
names(cancer_anno)= c(M_label$cancer_types %>% as.character() %>% unique())

M_label_colors[,"mutation_burden_class"] %>% unique() -> mutaion_anno
names(mutaion_anno)= c(M_label$mutation_burden_class %>% as.character() %>% unique())

M_label_colors[,"PFS.status"] %>% unique() -> survival_anno
names(survival_anno)= c(M_label$PFS.status %>% as.character() %>%unique())

M_label_colors[,"spectralClustering"] %>% unique() -> cluster_anno
names(cluster_anno)= c(M_label$spectralClustering %>% unique())

col_anno <- HeatmapAnnotation(df=M_label,
                              col = list("spectralClustering"=cluster_anno,
                                         "PFS.status"=survival_anno,
                                         "mutation_burden_class"=mutaion_anno,
                                         "cancer_types"=cancer_anno,
                                         "TumorPurity"=circlize::colorRamp2(c(0,
                                                                              min(M_label$TumorPurity[M_label$TumorPurity>0]),
                                                                              max(M_label$TumorPurity[M_label$TumorPurity>0])),
                                                                            c("white","grey100","grey0"))),
                              
                              width = unit(0.5, "cm"))
draw(col_anno,1:20)

library(circlize)
he = Heatmap(W[ind, ind],
             col = RColorBrewer::brewer.pal(9,"YlOrRd"),
             show_row_names = FALSE, 
             show_column_names = FALSE,
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             show_row_dend = FALSE, # whether show row clusters.
             top_annotation = col_anno,show_heatmap_legend = F
             # heatmap_legend_param = list(title = c("Scaled Exp."))
)
tiff(file.path("/project/huff/huff/immune_checkpoint/result_20171025/subtype_analysis/combine",paste(C,"clusters_group_heatmap.tiff",sep = "_")),width = 6,height = 6)
he
dev.off()


# survival ----------------------------------------------------------------
res_path <- "/project/huff/huff/immune_checkpoint/result_20171025/subtype_analysis/combine"
M_label %>%
  dplyr::rename("cluster" = "group","barcode"="sample") %>%
  dplyr::mutate(group = ifelse(cluster %in% c(1,3,5,6), "Group1","Group2")) %>%
  dplyr::mutate(group = ifelse(cluster == 2, "Group3", group)) %>%
  dplyr::select(barcode,cluster,group) -> C6_into_Group3_by_survival

C6_into_Group3_by_survival %>%
  dplyr::inner_join(time_status,by="barcode") %>%
  dplyr::rename("time"="PFS.time","status"="PFS") -> C6_into_Group3_by_survival.PFS

# survival with cluster
color_list = rainbow(6)
sur_name <- paste("Combined_Survival_for",6,"Cluster.png",sep="_")
C6_into_Group3_by_survival.PFS %>%
  dplyr::select(-group) %>%
  dplyr::rename("group"="cluster") %>%
  fn_survival("6 Clsuter PFS",color_list,sur_name,res_path)

# cluster with groups
color_list = rainbow(6)[2:4]
sur_name <- paste("Combined_Survival_for",3,"Groups.png",sep="_")
fn_survival(C6_into_Group3_by_survival.PFS,"3 Group PFS",color_list,sur_name,res_path)


# mutation burden ---------------------------------------------------------

M_label %>%
  dplyr::rename("cluster" = "group","barcode"="sample") %>%
  dplyr::left_join(mutation_burden_class,by="barcode") %>%
  dplyr::mutate(sm_count = ifelse(is.na(sm_count),0,sm_count)) %>%
  dplyr::rename("cancer_types"="cancer_types.x") %>%
  dplyr::mutate(group = ifelse(cluster %in% c(1,3,5,6), "Group1","Group2")) %>%
  dplyr::mutate(group = ifelse(cluster == 2, "Group3", group))-> group_cluster_mutation
comp_mat <- combn(length(unique(group_cluster_mutation$group)),2)
comp_list <- list()
m_name <-  paste("Combined_mutation_for",3,"Groups.png",sep="_")
for (col in 1:ncol(comp_mat)) {
  comp_list[[col]] <- comp_mat[,col]
}
color_list = rainbow(6)[2:4]
group_cluster_mutation %>%
  dplyr::mutate(sm_count = log2(sm_count)) %>%
  fn_mutation_burden("group","sm_count",color_list,"log2MutationBurden",comp_list,m_name,res_path)

# Figure 4
m_a_name <-  paste("Combined_all_mutation_for",3,"Groups.png",sep="_")
group_cluster_mutation %>%
  dplyr::mutate(sm_count = log2(sm_count)) %>%
  fn_mutation_burden_all("group","sm_count",color_list,"log2MutationBurden",comp_list,m_a_name,res_path)


# tumor purity ------------------------------------------------------------
m_name <-  paste("Combined_TumorPurity_for",3,"Groups.pdf",sep="_")
group_cluster_mutation %>%
  fn_mutation_burden("group","TumorPurity",color_list,"Tumor Purity",comp_list,m_name,res_path)

m_a_name <-  paste("Combined_all_TumorPurity_for",3,"Groups.pdf",sep="_")
group_cluster_mutation %>%
  fn_mutation_burden_all("group","TumorPurity",color_list,"Tumor Purity",comp_list,m_a_name,res_path)


# immune infiltration  -----------------------------------------------------
# wait for miao's data (TCAP)
TCGA_infiltration_data <- readr::read_rds(file.path("/project/huff/huff/data/TCGA/immune_infiltration/miao_TCAP_prediction","pancan33_immune_infiltration_by_TCAP.rds.gz")) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::select(barcode,InfiltrationScore)

m_name <-  paste("Combined_ImmuneInfiltration_for",3,"Groups.pdf",sep="_")
group_cluster_mutation %>%
  dplyr::left_join(TCGA_infiltration_data,by="barcode") %>%
  dplyr::filter(!is.na(InfiltrationScore)) %>%
  dplyr::mutate(InfiltrationScore = as.numeric(InfiltrationScore)) %>%
  fn_mutation_burden("group","InfiltrationScore",color_list,"Immune Infiltration Score",comp_list,m_name,res_path,h=12)

m_a_name <-  paste("Combined_all_ImmuneInfiltration_for",3,"Groups.pdf",sep="_")
group_cluster_mutation %>%
  dplyr::inner_join(TCGA_infiltration_data,by="barcode") %>%
  dplyr::mutate(InfiltrationScore = as.numeric(InfiltrationScore)) %>%
  fn_mutation_burden_all("group","InfiltrationScore",color_list,"Immune Infiltration Score",comp_list,m_a_name,res_path)


# 10 cancer related pathways score ----------------------------------------

TCGA_rppa <- readr::read_rds(file.path("/data/TCGA/TCGA_data","pancan32_rppa_score.rds.gz")) %>%
  tidyr::unnest() 

group_cluster_mutation %>%
  dplyr::inner_join(TCGA_rppa,by="barcode") %>%
  dplyr::select(barcode,cluster,cancer_types.x,group,pathway,score) -> pathway_score_for_groups

pathway_score_for_groups %>%
  ggpubr::ggboxplot(x = "group", y = "score",
                    color = "group" #add = "jitter",#, palette = "npg"
  ) +
  # geom_point(aes(x=as.numeric(group)+b,y=sm_count,color=group),alpha = 0.5) +
  scale_x_discrete(#breaks = c(1:3),
                   labels = c(1:3)
                   # expand = c(0.2,0.2,0.2)
  ) +
  facet_grid(cancer_types.x~ pathway, scales = "free") +
  # facet_wrap(~ pathway, strip.position = "bottom", scales = "free") +
  scale_color_manual(
    values = color_list
  )+
  # ylim(4,12) +
  ylab("Pathway Score") +
  xlab("Group") +
  theme(strip.background = element_rect(fill = "white",colour = "white"),
        legend.position = "none",
        text = element_text(size = 5),
        strip.text = element_text(size = 8)) +
  # ggpubr::stat_compare_means(label.y = 14,paired = TRUE) +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label = "p.signif")

p_name <- paste("Combined_PathwayScore_for",3,"Groups.pdf",sep="_")
ggsave(filename =p_name, path = res_path,device = "pdf",width = 6,height = 6)

p_a_name <- paste("Combined_all_PathwayScore_for",3,"Groups.pdf",sep="_")
ggsave(filename =p_a_name, path = res_path,device = "pdf",width = 6,height = 15)

