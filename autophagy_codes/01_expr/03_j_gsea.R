library(methods)
library(magrittr)

tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
#output path
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")
# script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
script_path <- "/home/cliu18/liucj/github/RstudioWithGit/autophagy_codes/01_expr"
#
gene_list_expr <- readr::read_rds(file.path(expr_path_a, ".rds_03_a_gene_list_expr.rds.gz"))
gene_list_fc_pvalue_simplified <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_a_atg_lys_fc_pvalue_simplified.rds.gz"))

gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))

gene_sets <- list(
  lc3 = dplyr::filter(gene_list, detail == "LC3"),
  ulk_complex = dplyr::filter(gene_list, detail == "ULK complex"),
  pi3k = dplyr::filter(gene_list, detail == "PI3K III complex"),
  atg_lys = dplyr::filter(gene_list, status %in% c("a", "l")),
  atg = dplyr::filter(gene_list, type == "atg"),
  lys = dplyr::filter(gene_list, type == "lys"),
  pathway = dplyr::filter(gene_list, type == "pathway"),
  atg_core = dplyr::filter(gene_list, pathway == "autophagesome formation-core"),
  atg_form = dplyr::filter(gene_list, pathway == "autophagesome formation"),
  lys_com = dplyr::filter(gene_list, pathway == "lysosome composition"),
  lys_deg = dplyr::filter(gene_list, pathway == "lysosome degradation")
) %>% purrr::map("symbol")

gsea_path <- file.path(expr_path_a, "gsea")
if (!file.exists(gsea_path)) { dir.create(gsea_path) }

# create gmt ----------------------------------------------
fn_gmt <- function(.x, .path = gsea_path){
  .file_name <- file.path(.path,"atg_lys_pathway.gmt")
  .x %>% purrr::map_int(length) %>% unlist() %>% max() -> .max_length
  
  .x %>% 
    purrr::map(.f = function(x, .m){.y <- c(x, rep(NA, .m - length(x))); names(.y) <- 1:.m; .y}, .m = .max_length) %>% 
    purrr::map(.f = tibble::enframe) %>% 
    tibble::enframe() %>% 
    tidyr::unnest() %>% 
    tidyr::spread(key = name1, value = value) %>% 
    dplyr::select(1, as.character(1:.max_length)) %>% 
    tibble::add_column(des = "cj_atg_lys", .before = 2) -> .x_write
  
  .x_write %>% readr::write_delim(path = .file_name, delim = "\t", col_names = F, na = "")
}
gene_sets %>% fn_gmt(.path = gsea_path)

# create gct cls file------------------------------------
fn_gct_cls <- function(.x, .y, .path = gsea_path){
  # .x <- te$cancer_types
  # .y <- te$filter_expr[[1]]
  .y %>% tidyr::drop_na() -> .y
  
  tibble::tibble(barcode = names(.y)[-c(1,2)], type = stringr::str_sub(names(.y)[-c(1,2)], start = 14, 15)) %>% 
    dplyr::filter(type %in% c("01", "11")) %>% 
    dplyr::arrange(type) %>% 
    dplyr::pull(barcode) -> .names
  
  .y %>% 
    dplyr::select(1,2, .names) %>% 
    dplyr::rename(NAME = symbol, DESCRIPTION = entrez_id) -> .y_gct
  
  .gct_filename <- file.path(.path, paste(.x,'mRNA_expression.gct', sep = "_"))
  
  # write gct files
  readr::write_lines(x = "#1.2", path = .gct_filename)
  readr::write_tsv(x = data.frame(gene = nrow(.y), sample = sum(.tn)), path = .gct_filename, append = T, col_names = F)
  readr::write_tsv(x = .y_gct, path = .gct_filename, append = T, col_names = T)
  
  .cls_filename <- file.path(.path, paste(.x, "mRNA_expression.cls", sep = "_"))
  readr::write_delim(x = data.frame(sample = length(.names), 2, 1), path = .cls_filename, col_names = F,  delim  = " ")
  readr::write_delim(x = data.frame("#", "Tumor", "Normal"), path = .cls_filename, col_names = F, append = T, delim = " ")
  .names %>% 
    stringr::str_sub(start = 14, 15) %>% 
    dplyr::recode("01" = "T", "11" = "N") %>% 
     as.matrix() %>% t() %>% as.data.frame() -> .label
  readr::write_delim(x = .label, path = .cls_filename, delim = " ", col_names = F, append = T)
}
purrr::walk2(.x = gene_list_expr$cancer_types, .y = gene_list_expr$filter_expr, .f = fn_gct_cls, .path = gsea_path)

# run gsea--------------------------------------------------------------------

fn_atleast2_normal <- function(cancer_types, filter_expr){
  names(filter_expr)[-c(1,2)] ->.barcode
  
  sum(.barcode %>% stringr::str_sub(14,15) %>% stringr::str_detect("11")) -> .n_normal

  if(.n_normal >= 10){
    return(cancer_types)
  } else{
    return(NULL)
  }
}
cancers <- gene_list_expr %>% purrr::pmap(.f = fn_atleast2_normal) %>% unlist()

fn_gsea <- function(.ds, .cls, .db, .output, .doc){
  # .ds is the gct format file
  # .cls is the phenotype format
  # .db is the gmt file format
  
  GSEA( # Input/Output Files :-------------------------------------------
        input.ds =  .ds,               # Input gene expression Affy dataset file in RES or GCT format
        input.cls = .cls,               # Input class vector (phenotype) file in CLS format
        gs.db = .db,           # Gene set database in GMT format
        output.directory = .output,            # Directory where to store output and results (default: "")
        #  Program parameters :----------------------------------------------------------------------------------------------------------------------------
        doc.string            = .doc,     # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
        non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
        reshuffling.type      = "gene.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels"
        nperm                 = 1000,            # Number of random permutations (default: 1000)
        weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
        nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
        fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
        fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
        topgs                 = 50,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
        adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
        gs.size.threshold.min = 10,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
        gs.size.threshold.max = 1000,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
        reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
        preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
        random.seed           = 111,             # Random number generator seed. (default: 123456)
        perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
        fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
        replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
        save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
        OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
        use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
  )
}
fn_run_gsea <- function(.x, .path = gsea_path, script_path = script_path){
  .gct <- file.path(.path, paste(.x, "mRNA_expression.gct", sep = "_"))
  .cls <- file.path(.path, paste(.x, "mRNA_expression.cls", sep = "_"))
  .gmt <- file.path(.path, "atg_lys_pathway.gmt")
  .output_dir <- file.path(.path, "gsea_result/")
  .doc <- paste(.x, "GSEA.analysis", sep = ".")
  if(!file.exists(.output_dir)){ dir.create(.output_dir)}
  
  source(file.path(script_path, "GSEA.1.0.r"))
  fn_gsea(.ds = .gct, .cls = .cls, .db = .gmt, .output = .output_dir, .doc = .doc)
}

tibble::tibble(cancer_types = cancers) %>% 
  head(1) %>% 
  dplyr::mutate(res = purrr::walk(.x = cancer_types, .f = fn_run_gsea, .path = gsea_path, script_path = script_path))
# source GSEA

cluster <- multidplyr::create_cluster(length(cancers))
tibble::tibble(cancer_types = cancers) %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_run_gsea", fn_run_gsea)  %>%
  multidplyr::cluster_assign_value("fn_gsea", fn_gsea)  %>%
  multidplyr::cluster_assign_value("gsea_path", gsea_path)  %>%
  multidplyr::cluster_assign_value("script_path", script_path)  %>%
  dplyr::mutate(res = purrr::walk(.x = cancer_types, .f = fn_run_gsea, .path = gsea_path, script_path = script_path)) %>% 
  dplyr::collect()
parallel::stopCluster(cluster)

gsea_result_path <- file.path(gsea_path, "gsea_result")

fn_load_es <- function(.x, .path = gsea_result_path){
  .filename_T <- paste(.x, "GSEA.analysis.SUMMARY.RESULTS.REPORT.T.txt", sep = ".")
  .t <- readr::read_tsv(file = file.path(.path, .filename_T)) %>% 
    dplyr::select(GS, ES, NES, p_val = `NOM p-val`, q_val = `FDR q-val`)
  
  .filename_N <- paste(.x, "GSEA.analysis.SUMMARY.RESULTS.REPORT.N.txt", sep = ".")
  .n <- readr::read_tsv(file = file.path(.path, .filename_N)) %>% 
    dplyr::select(GS, ES, NES, p_val = `NOM p-val`, q_val = `FDR q-val`)
  
  dplyr::bind_rows(.t, .n) %>% 
    tibble::add_column(cancer_types = .x, .before = 1)
}

cancers %>% 
  purrr::map(.f = fn_load_es, .path = gsea_result_path) %>% 
  dplyr::bind_rows() -> cancers_es

readr::write_rds(cancers_es, path = file.path(expr_path_a, ".rds_02_i_gsea_cancer_es.rds.gz"), compress = "gz")
library(ggplot2)
cancers_es %>% 
  dplyr::filter(q_val < 0.05) %>%
  ggplot(aes(x = cancer_types, y = GS)) +
  geom_point(aes(size = -log10(q_val), color = NES)) +
  scale_color_gradient2(low = 'blue', high = 'red', mid = 'white')

save.image(file = file.path(expr_path_a, ".rda_03_j_gsea.rda"))
load(file = file.path(expr_path_a, ".rda_03_j_gsea.rda"))

