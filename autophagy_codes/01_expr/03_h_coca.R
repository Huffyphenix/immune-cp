#load library
library(methods)
library(magrittr)
library(CancerSubtypes)


# Integrative clustering analysis
# Use mRNA, methylation and copy number
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
pancan_color <- readr::read_tsv(file.path(tcga_path, 'PanCanAtlasTumors_color_coded_by_organ_system_20170302.tsv')) 
pcc <- pancan_color %>% dplyr::pull(color)

tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
names(pcc) <- pancan_color %>% dplyr::pull(cancer_types)

# mRNA clustering
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")

gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))

atg_gene <- gene_list %>% dplyr::filter(type == "atg") %>%  dplyr::pull(symbol)
lys_gene <- gene_list %>% dplyr::filter(type == "lys") %>% dplyr::pull(symbol)

expr <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_a_gene_list_expr.rds.gz"))

fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
}
fn_filter_tumor <- function(filter_expr){
  # filter_expr <- te$filter_methy[[1]]
  # print(cancer_types)
  filter_expr %>% 
    # dplyr::select(-entrez_id) %>%
    # dplyr::select(-gene) %>%
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::distinct(symbol, barcode, .keep_all = T) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type != 11) %>% 
    dplyr::select(-type) %>% 
    tidyr::spread(key = barcode, value = expr) 
}
fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
}

# expr %>%
#   dplyr::mutate(filter_expr = purrr::map(.x = filter_expr, .f = fn_filter_tumor)) %>%
#   tibble::deframe() %>%
#   dplyr::bind_cols() %>%
#   dplyr::select(-dplyr::matches('symbol\\d+')) -> expr_mat
# 
# expr_mat %>% readr::write_rds(path = file.path(expr_path_a, ".rds_03_h_coca_expr_mat.rds.gz"), compress = 'gz')
expr_mat <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_h_coca_expr_mat.rds.gz"))

# copy number data
cnv_path <- "/home/cliu18/liucj/projects/6.autophagy/03_cnv"
cnv <- readr::read_rds(path = file.path(cnv_path, ".rds_02_cnv_a_gene_list_raw.rds.gz"))

# cnv %>%
#   dplyr::mutate(filter_cnv = purrr::map(.x = filter_cnv, .f = fn_filter_tumor)) %>%
#   tibble::deframe() %>%
#   dplyr::bind_cols() %>%
#   dplyr::select(-dplyr::matches('symbol\\d+')) -> cnv_mat
# cnv_mat %>% readr::write_rds(path = file.path(expr_path_a, ".rds_03_h_coca_cnv_mat.rds.gz"), compress = 'gz')
cnv_mat <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_h_coca_cnv_mat.rds.gz"))

# methylation data
methy_path <- "/home/cliu18/liucj/projects/6.autophagy/06_methylation"
methy <- readr::read_rds(path = file.path(methy_path, ".rds_02_methy_a_gene_list_methy.rds.gz"))

# methy %>% 
#   dplyr::mutate(filter_methy = purrr::map(.x = filter_methy, .f = fn_filter_tumor)) %>%
#   tibble::deframe() %>% 
#   dplyr::bind_cols() %>% 
#   dplyr::select(-dplyr::matches('symbol\\d+')) -> methy_mat
# 
# methy_mat %>% readr::write_rds(path = file.path(expr_path_a, ".rds_03_h_coca_methy_mat.rds.gz"), compress = 'gz')
methy_mat <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_h_coca_methy_mat.rds.gz"))

# unify samples
colnames(expr_mat) %<>% fun_barcode
expr_mat[!colnames(expr_mat) %>% duplicated()]-> expr_mat_unique
colnames(cnv_mat) %<>% fun_barcode 
cnv_mat[!colnames(cnv_mat) %>% duplicated()]-> cnv_mat_unique
colnames(methy_mat) %<>% fun_barcode 
methy_mat[!colnames(methy_mat) %>% duplicated()] -> methy_mat_unique

common_names <- intersect(colnames(expr_mat_unique), intersect(colnames(cnv_mat_unique), colnames(methy_mat_unique)))

expr_matrix <- 
  expr_mat_unique[common_names] %>% 
  dplyr::filter(!is.na(`TCGA-P7-A5NY`)) %>% 
  .[-1] %>% as.matrix() %>% 
  data.normalization(type = "feature_zscore", log2 = T)

cnv_matrix <- 
  cnv_mat_unique[common_names] %>% 
  dplyr::filter(!is.na(`TCGA-P7-A5NY`)) %>%
  .[-1] %>% as.matrix()

methy_matrix <- 
  methy_mat_unique[common_names] %>% 
  dplyr::filter(!is.na(`TCGA-P7-A5NY`)) %>%
  .[-1] %>% as.matrix() %>% 
  data.imputation(fun = "median")

coca <- list(mrna = expr_matrix, cnv = cnv_matrix, methy = methy_matrix)

# cc for every one
# expr_cc <- ExecuteCC(clusterNum = 3, d = expr_matrix, maxK = 10, clusterAlg = "hc", distance = "pearson", title = "expr_matrix", plot = F)
# expr_cc %>% readr::write_rds(file.path(expr_path_a, ".rds_03_h_coca_expr_cc.rds.gz"), compress = 'gz')
# 
# cnv_cc <- ExecuteCC(clusterNum = 3, d = cnv_matrix, maxK = 10, clusterAlg = "hc", distance = "pearson", title = "cnv_matrix", plot = F)
# cnv_cc %>% readr::write_rds(file.path(expr_path_a, ".rds_03_h_coca_cnv_cc.rds.gz"), compress = 'gz')
# 
# methy_cc <- ExecuteCC(clusterNum = 3, d = methy_matrix, maxK = 10, clusterAlg = "hc", distance = "pearson", title = "methy_matrix", plot = F)
# methy_cc %>% readr::write_rds(file.path(expr_path_a, ".rds_03_h_coca_methy_cc.rds.gz"), compress = 'gz')



# coca_cc <- ExecuteCC(clusterNum=3,d=coca,maxK=10,clusterAlg="hc",distance="pearson")
# coca_snf <- ExecuteSNF(coca, clusterNum=3, K=20, alpha=0.5, t=20, plot = F)
# coca_snf %>% readr::write_rds(path = file.path(expr_path_a, ".rds_03_h_coca_coca_snf.rds.gz"), compress = "gz")
coca_snf <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_h_coca_coca_snf.rds.gz"))

d <- coca_snf$distanceMatrix %>% as.matrix()

fn_best_clust <- function(k, d){
  # k = 3
  group <- spectralClustering(d, K = k)
  
  diag(d) <- 0
  diag(d) <- max(d)
  
  tibble::tibble(
    name = paste("V", c(1:ncol(d)), sep = ""), 
    group = group) %>% 
    dplyr::arrange(group)  -> .rank_sample
  
  d %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    tibble::add_column(
      sample1 = paste("V", c(1:ncol(d)), sep = ""), 
      .before = 1) %>% 
    tidyr::gather(key = sample2, value = sim, -sample1) -> .plot_ready
  
  .plot_ready %>% 
    dplyr::mutate(sim = ifelse(sim > quantile(sim, 0.75), 1, 0)) %>%
    ggplot(aes(x= sample1, y = sample2, fill = sim)) +
    geom_tile() +
    scale_x_discrete(limits = .rank_sample$name) +
    scale_y_discrete(limits = .rank_sample$name) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      
      legend.position = "none"
    )  +
    geom_vline(xintercept = table(group) %>% as.vector() %>% purrr::accumulate(`+`), color = "red") +
    geom_hline(yintercept = table(group) %>% as.vector() %>% purrr::accumulate(`+`), color = "red") -> p
  
  ggsave(filename = paste(k, "coca_clust.tif", sep = "_"), plot = p, device = "tiff", width = 8, height = 8, path = "/home/cliu18/liucj/github/RstudioWithGit/autophagy_codes/01_expr/coca_clust_trial")
  
  return(1)
}

cluster <- multidplyr::create_cluster(6)
tibble::tibble(k = 2:7) %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("SNFtool") %>%
  multidplyr::cluster_library("ggplot2") %>%
  multidplyr::cluster_assign_value("fn_best_clust", fn_best_clust)  %>%
  multidplyr::cluster_assign_value("d", d)  %>%
  dplyr::mutate(a = purrr::walk(.x = k, .f = fn_best_clust, d = d)) 
parallel::stopCluster(cluster)

SNFtool::spectralClustering(d, K = ) -> group
tibble::tibble(
  sample = common_names[-1], 
  name = paste("V", c(1:ncol(coca_snf$distanceMatrix)), sep = ""), 
  group = coca_snf$group) %>% 
  dplyr::arrange(group)  -> rank_sample

class(coca_snf$distanceMatrix) <- "matrix"

coca_snf$distanceMatrix %>% 
  as.data.frame() %>% 
  tibble::as_tibble() %>% 
  tibble::add_column(
    sample1 = paste("V", c(1:ncol(coca_snf$distanceMatrix)), sep = ""), 
    .before = 1) %>% 
  tidyr::gather(key = sample2, value = sim, -sample1) -> plot_ready

plot_ready %>% 
  ggplot(aes(x= sample1, y = sample2, fill = sim)) +
  geom_tile() +
  scale_x_discrete(limits = rank_sample$name) +
  scale_y_discrete(limits = rank_sample$name) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    
    legend.position = "none"
  ) -> p

ggsave(filename = "coca_snf.tif", plot = p, device = "tiff", path = expr_path_a, width = 8, height = 8)
#

#
save.image(file = file.path(expr_path_a, ".rda_03_h_coca.rda"))
load(file = file.path(expr_path_a, ".rda_03_h_coca.rda"))
