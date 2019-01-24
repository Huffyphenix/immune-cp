#load library
library(methods)
library(magrittr)
library(CancerSubtypes)
library(ConsensusClusterPlus)
library(SNFtool)
# load data ---------------------------------------------------------------
data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"

## get from 1.data_prepare.R
load(file = file.path(data_result_path, ".rda_genelist_data_survival_cancer.info.rda"))

# prepare clinical data ---------------------------------------------------
mutation_burden_class <- readr::read_rds(file.path("/project/huff/huff/data/TCGA","classfication_of_26_cancers_by_mutation_burden192.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::rename("cancer_types"="Cancer_Types")
time_status %>%
  dplyr::full_join(mutation_burden_class,by="barcode") %>%
  dplyr::full_join(genelist_methy_mutaion_class.cancer_info,by="barcode") %>%
  dplyr::mutate(cancer_types=ifelse(is.na(cancer_types.x),cancer_types.y,cancer_types.x)) %>%
  dplyr::select(-cancer_types.x,-cancer_types.y) %>%
  dplyr::rename("time"="PFS.time","status"="PFS")-> clinical_data

# load all methy data ------
immune_path <- "/project/huff/huff/immune_checkpoint"
tcga_path <- file.path(immune_path,"/data/TCGA_data")

mthy <- readr::read_rds(file.path(tcga_path,"pancan33_meth.rds.gz"))

result_path <- "/project/huff/huff/immune_checkpoint/result_20171025/subtype_analysis/methy/comfirm_use_random_genelist"

# gene pool ----
mthy$methy[[1]] %>% .$symbol -> gene_pool

# samples pool ----
# include 10 cancers: BLCA,BRCA,HNSC,KIRC,LGG,LUAD,PRAD,THCA,UCEC,SKCM ----
load(file = file.path(data_result_path, ".rda_genelist_data_survival_cancer.info.rda"))
cancers_needed <- c("BLCA","BRCA","HNSC","KIRC","LGG","LUAD","PRAD","THCA","UCEC","SKCM")
genelist_methy_mutaion_class.cancer_info %>%
  dplyr::filter( cancer_types %in% cancers_needed) %>%
  dplyr::inner_join(time_status,by="barcode") %>%
  .$barcode -> sample_pool

set.seed(464134343)
seed.pool <- sample(10000:100000, size = 10)

for (i in seed.pool) {
  # get a random gene list with 67 population ----
  set.seed(seed.pool[i])
  genelist <- gene_pool[sample(1:length(gene_pool), size = 67)]
  
  # get a random samples with 1000 population ----
  set.seed(seed.pool[i])
  samplelist <- sample_pool[sample(1:length(sample_pool), size = 500)]
  
  # get methy data matrix ----
  mthy %>%
    dplyr::filter(cancer_types %in% cancers_needed) %>%    # some cancers in TCGA merged serveral cancers into one
    tidyr::unnest() %>%
    dplyr::filter(symbol %in% genelist) %>%
    dplyr::select(-gene) %>%
    tidyr::gather(-cancer_types,-symbol,key=barcode,value=methy) %>%
    dplyr::mutate(barcode = substr(barcode,1,12)) %>%
    dplyr::filter(barcode %in% samplelist) %>%
    # dplyr::inner_join(mutation_burden_class,by="barcode") %>%
    dplyr::filter(! is.na(methy)) %>% 
    dplyr::arrange(barcode) -> genelist_methy_mutaion_class
  
  genelist_methy_mutaion_class %>%
    dplyr::select(symbol,barcode,methy) %>%
    dplyr::group_by(symbol,barcode) %>%
    dplyr::mutate(methy = mean(methy)) %>%
    unique() %>%
    dplyr::ungroup() %>%
    tidyr::spread(key = barcode,value = methy) -> gene_list_methy_spread
  
  gene_list_methy_matrix <- as.matrix(gene_list_methy_spread[,-1])
  rownames(gene_list_methy_matrix) <- gene_list_methy_spread$symbol
  
  # data matrix manage -----
  d <- gene_list_methy_matrix
  index=which(is.na(d)) 
  d=data.imputation(d,fun="median")
  
  dc = sweep(d,1, apply(d,1,median,na.rm=T)) ## median center genes
  
  
  # same as above but with pre-computed distance matrix, useful for large datasets (>1,000's of items)
  dt = as.dist(1-cor(dc,method="pearson"))
  

# do cluster --------------------------------------------------------------
  results = ConsensusClusterPlus(dt,maxK=10,reps=100,pItem=0.8,pFeature=1,title="methy_CC",distance="pearson",clusterAlg="hc",seed=1262118388.71279)
  result_name <- paste(i,"th.comfirm.CC10.result.rds.gz",sep = "")
  results %>%
    readr::write_rds(file.path(result_path,result_name),compress = "gz")
  
# cluster result analysis -------------------------------------------------
  C=10
  W <- results[[C]][['consensusMatrix']]
  group <- results[[C]][['consensusClass']]
  W <- matrix(W,nrow = length(group),dimnames = list(names(group),names(group)))
  data.frame(sample = names(group),group=group) %>%
    readr::write_tsv(file.path(result_path,paste(i,"th_comfirm_for",C,"clusters_group_info.tsv",sep=".")))
  
  # statistic of the cancer composition -----
  data.frame(barcode = names(group),group=group) %>%
    dplyr::as.tbl() %>%
    dplyr::left_join(clinical_data,by="barcode") %>%
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
  
  # draw ggplot pic, stacked bar plot -----
  cluster_cancers_statistic %>%
    dplyr::mutate(group=as.character(group)) %>%
    ggplot(aes(x=cancer_types,y=`Sample Composition (%)`,width=n_a)) +
    geom_bar(aes(fill=group),stat="identity",color= "grey")+
    facet_wrap(~cancer_types,nrow = 1) +
    scale_fill_manual(
      name = "Clusters",
      breaks = as.character(c(1:10)),
      values = c("#CDC0B0", "#838B8B", "#000000", "#0000FF", "#00008B", "#8A2BE2", "#A52A2A", "#FF4040", "#98F5FF", "#EEAD0E")
    ) +
    theme(
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
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line.x = element_blank()
    )
  ggsave(filename =paste(i,"th_comfirm_for","Sample_composition_for",C,"Clusters-stacked.png",sep="."), path = result_path,device = "png",height = 4,width = 6)
}


