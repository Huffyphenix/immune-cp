library(methods)
library(magrittr)
library(SNFtool)
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
pancan_color <- readr::read_tsv(file.path(tcga_path, 'PanCanAtlasTumors_color_coded_by_organ_system_20170302.tsv'))
pcc <- pancan_color %>% dplyr::pull(color)

tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
names(pcc) <- pancan_color %>% dplyr::pull(cancer_types)

# mRNA clustering
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")
load(file = file.path(expr_path_a, ".rda_03_h_coca.rda"))
coca_snf <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_h_coca_coca_snf.rds.gz"))

class(coca_snf$distanceMatrix) <- "matrix"

d <- coca_snf$distanceMatrix

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
