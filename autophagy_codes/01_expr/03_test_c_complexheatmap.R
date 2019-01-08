library(magrittr)
library(ggplot2)

expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
heatmap_path <- file.path(expr_path, "03_test_c_complexheatmap")
expr_path <- file.path(expr_path, "03_a_gene_expr")

gene_list <- readr::read_rds(file.path(expr_path, "rds_03_at_ly_comb_gene_list.rds.gz"))
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))

fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
} #get short barcode from long barcode
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
}

fn_extract_pair <- function(cancer_types, filter_expr){
  print(cancer_types)
  filter_expr %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type %in% c("01", "11")) %>%
    dplyr::mutate(sample = fun_barcode(barcode)) %>% 
    dplyr::mutate(type = plyr::revalue(
      x = type,
      replace = c("01" = "Tumor", "11" = "Normal"),
      warn_missing = F
    )) %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(n() >= 2, length(unique(type)) == 2) %>%
    dplyr::ungroup() -> d
  
  if(nrow(d) == 0 ){return(NULL)}
  
  d %>% 
    dplyr::mutate(expr = log2(expr + 1)) %>% 
    dplyr::select(-barcode) %>% 
    tidyr::unite(col = sample, type, sample, sep = "_") %>% 
    tidyr::spread(key = sample, value = expr) %>% 
    dplyr::mutate_all(.funs = dplyr::funs(ifelse(is.na(.), 0.001, .))) -> samples
  
  if (samples %>% colnames() %>% unique() %>% length() -1 < 20) {
    return(NULL)
  }
  
  tibble::tibble(cancer_types = cancer_types, filter_expr = list(samples))
  
  }


gene_list_expr %>% 
  purrr::pmap(.f = fn_extract_pair) %>% 
  dplyr::bind_rows() -> gene_list_expr_pair
  


gene_list_expr_pair %>% 
  head(1) %>% 
  tidyr::unnest() -> d
  
mat <- as.matrix(d[, grep("TCGA", colnames(d))])
rownames(mat) <- d$symbol
base_mean = rowMeans(mat)
mat_scaled = t(apply(mat, 1, scale))
type = stringr::str_extract(string = colnames(mat), pattern = "\\w+_") %>% 
  stringr::str_replace("_","")

library(ComplexHeatmap)
library(circlize)
ha = HeatmapAnnotation(df = data.frame(type = type))
Heatmap(
  head(mat_scaled, 100),
  name = "expression",
  km = 2,
  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
  top_annotation = ha,
  top_annotation_height = unit(4, "mm"),
  show_row_names = T,
  show_column_names = T
) 
save.image(file = file.path(heatmap_path, ".rda_03_test_c_complexheatmap.rda"))
load(file = file.path(heatmap_path, ".rda_03_test_c_complexheatmap.rda"))

