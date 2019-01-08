
# library -----------------------------------------------------------------
library(magrittr)
library(ggplot2)
library(GSVA)
# library(msigdf)
library(biomaRt)

load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c1_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c3_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c4_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c6_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c5_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c7_v5p2.rdata"))
# mouse data
# load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_H_v5p2.rdata"))
# load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c1_v5p2.rdata"))
# load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p2.rdata"))
# load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c3_v5p2.rdata"))
# load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c4_v5p2.rdata"))
# load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c6_v5p2.rdata"))
# load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c5_v5p2.rdata"))
# load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c7_v5p2.rdata"))
msigdf <- list()
for (i in ls(pattern = "^Hs\\.|^Mm\\.")) {
  message(i)
  msigdf[[i]] <- 
    eval(parse(text = i)) %>%
    head(100) %>% 
    plyr::ldply(function(x) dplyr::data_frame(entrez = x), .id = "geneset") %>%
    dplyr::filter(entrez != "-") %>%
    dplyr::mutate(entrez = as.integer(entrez), geneset = as.character(geneset)) %>%
    dplyr::tbl_df()
}; rm(i)

msigdf <- 
  msigdf %>%
  dplyr::bind_rows(.id = "x") %>%
  tidyr::separate(x, into = c("org", "collection")) %>%
  dplyr::mutate(org = dplyr::recode(org, Hs = "human", Mm = "mouse")) %>%
  dplyr::mutate(collection = gsub("H", "hallmark", collection))
msigdf.human <- msigdf %>% dplyr::filter(org == "human") %>% dplyr::select(-org)

# path --------------------------------------------------------------------
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")

# load data ---------------------------------------------------------------
expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz"))
gene_list_gsva <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_i_gsva_gene_list_gsva.rds.gz"))

# MsigDB entrez to symbol -------------------------------------------------
# name c1-c7
msigdf.human$collection %>% unique() -> collection
name <- c("C1_POSITIONAL","C2_CURATED","C3_MOTIF","C4_COMPUTATIONAL","C5_GENE_ONTOLOGY","C6_ONCOGENIC_SIGNATURES","C7_IMMUNOLOGIC_SIGNATURES","HALLMARK")
setNames(name, collection) -> name

# entrez to symbol
GENES = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")
GENES.ATTRIBUTES <- listAttributes(GENES)
GENES.FILTERS <- listFilters(GENES)
msigdf.human$entrez %>% unique() -> entrez
getBM(
  attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"), 
  filters = "entrezgene", values = entrez, mart = GENES) %>% 
  tidyr::drop_na() %>% 
  tibble::as_tibble() %>%
  dplyr::mutate(entrez = as.integer(entrezgene)) %>% 
  dplyr::select(symbol = hgnc_symbol, entrez) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(symbol != "") -> entrez_symbol 
msigdf.human %>% 
  dplyr::mutate(name = plyr::revalue(collection, name)) %>% 
  dplyr::inner_join(entrez_symbol, by = "entrez") -> human
readr::write_rds(human, path = file.path(tcga_path, "MSigDB_human_collections.rds.gz"), compress = "gz")

C2_CURATED <- dplyr::filter(human, collection == "c2")
C5_GENE_ONTOLOGY <- dplyr::filter(human, collection == "c5")
C6_ONCOGENIC_SIGNATURES <- dplyr::filter(human, collection == "c6")
C7_IMMUNOLOGIC_SIGNATURES <- dplyr::filter(human, collection == "c7")
HALLMARK <- dplyr::filter(human, collection == "hallmark")

# MSigDB set gsva ---------------------------------------------------------
fn_geneset <- function(.g){
  .g %>% 
    dplyr::select(geneset, entrez) %>% 
    tidyr::nest(entrez) %>% 
    tibble::deframe() %>% purrr::map("entrez")
}
human %>% dplyr::filter(collection %in% c("c2", "c5","c6", "c7", "hallmark")) -> human_sets
human_sets %>% fn_geneset() -> gene_sets
human_sets %>% dplyr::select(collection = name, geneset) %>% dplyr::distinct() -> cg
fn_gsva <- function(.x, .y, gene_sets, cg){
  # .x <- .te$cancer_types
  # .y <- .te$expr[[1]]
  print(.x)
  
  .y %>% 
    tidyr::drop_na() %>% 
    dplyr::select(-symbol) %>% 
    dplyr::rename(entrez = entrez_id) %>% 
    dplyr::distinct(entrez, .keep_all = T) -> .d
  
  .d_mat <- as.matrix(.d[,-1])
  rownames(.d_mat) <- .d$entrez
  
  .es_dif <- gsva(.d_mat, gene_sets, method = "gsva", mx.diff = TRUE, verbose = FALSE, parallel.sz = 1)
  
  .es_dif$es.obs %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    tibble::add_column(geneset = rownames(.es_dif$es.obs), .before = 1) -> .d_es
  
  cg %>% dplyr::inner_join(.d_es, by = "geneset")
}

# multidplyr
cluster <- multidplyr::create_cluster(nrow(expr))
expr %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("GSVA") %>%
  multidplyr::cluster_assign_value("fn_gsva", fn_gsva)  %>%
  multidplyr::cluster_assign_value("gene_sets", gene_sets) %>% 
  multidplyr::cluster_assign_value("cg", cg) %>% 
  dplyr::mutate(gsva = purrr::map2(.x = cancer_types, .y = expr, .f = fn_gsva, gene_sets = gene_sets, cg = cg)) %>%
  dplyr::select(-expr) %>% 
  dplyr::collect() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) -> MSigDB_GSVA_score
parallel::stopCluster(cluster)

readr::write_rds(x = MSigDB_GSVA_score, path = file.path(tcga_path, "pancan33_MSigDB_GSVA_score.rds.gz"), compress = "gz")

# MSigDB set gsva correlate with autophagy --------------------------------
gene_list_gsva <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_i_gsva_gene_list_gsva.rds.gz"))

gene_list_gsva %>%
  dplyr::inner_join(MSigDB_GSVA_score, by = "cancer_types") %>% 
  dplyr::rename(atg = gsva.x, msig = gsva.y) -> atg_msig

fn_cor <- function(.x, .y){
  # .x <- .te$atg[[1]]
  # .y <- .te$msig[[1]]
  
  .x %>% tidyr::gather(key = "barcode", value = "atg_gsva", -c(set)) -> .dx
  .y %>% tidyr::gather(key = "barcode", value = "msig_gsva", -c(collection, geneset)) -> .dy
  .dy %>% 
    dplyr::inner_join(.dx, by = "barcode") %>% 
    dplyr::filter(stringr::str_sub(barcode, 14, 15) != "11") -> .d
  
  .d %>%
    # dplyr::filter(collection == "HALLMARK") %>% 
    dplyr::group_by(collection, geneset, set) %>% 
    dplyr::do(
      broom::tidy(
        tryCatch(
          cor.test(~msig_gsva + atg_gsva, data = ., method = "spearman"),
          error = function(e){1}
        )
      )
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(collection, geneset, set, coef = estimate, pval = p.value) %>% 
    dplyr::mutate(fdr = p.adjust(pval, method = "fdr")) -> .d_corr
  
  .d_corr
}

cluster <- multidplyr::create_cluster(nrow(expr))
atg_msig %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("GSVA") %>%
  multidplyr::cluster_assign_value("fn_cor", fn_cor)  %>%
  dplyr::mutate(corr = purrr::map2(.x = atg, .y = msig, .f = fn_cor)) %>% 
  dplyr::collect() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) -> atg_msig_corr
parallel::stopCluster(cluster)
readr::write_rds(x = atg_msig_corr, path = file.path(expr_path_a, ".rda_04_gsva_atg_msig_corr.rds.gz"), compress = "gz")

atg_msig_corr <- readr::read_rds(path = file.path(expr_path_a, ".rda_04_gsva_atg_msig_corr.rds.gz"))

# Filter significante and plot --------------------------------------------
atg_msig_corr %>% 
  tidyr::unnest(corr) %>% 
  dplyr::filter(fdr < 0.05, abs(coef) > 0.3) -> atg_msig_corr_un_sig

atg_msig_corr_un_sig %>% 
  dplyr::filter(collection %in% c("HALLMARK")) %>%
  # dplyr::filter(stringr::str_detect(geneset, stringr::coll("MYC", ignore_case = T))) %>%
  dplyr::filter(set %in% c("atg", "lys")) %>%
  dplyr::select(1,2, 3, 4, 5,7) %>% 
  dplyr::mutate(direction = ifelse(coef > 0, "p", "n")) %>%
  dplyr::group_by(geneset, set, direction) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(key = direction, value = n, fill = 0) %>% 
  dplyr::mutate(t = p - n) %>% 
  dplyr::arrange(-t, -p, -n) %>% 
  dplyr::filter(n + p > 5) -> plot_ready

plot_ready %>% 
  dplyr::group_by(geneset) %>% 
  dplyr::summarise(m = mean(t)) %>% 
  dplyr::arrange(-m) %>% 
  dplyr::pull(geneset) -> geneset_rank

CPCOLS <- c("#000080", "#F5F5F5", "#CD0000")
plot_ready %>% 
  dplyr::mutate(s = glue::glue("neg: {n}, pos: {p}")) %>% 
  ggplot(aes(x = set, y = geneset)) +
  geom_tile(aes(fill = t)) +
  # geom_text(aes(label = s)) +
  scale_y_discrete(limits = geneset_rank) +
  scale_x_discrete(label = c("ATG", "LYS")) +
  scale_fill_gradient2(
    name = "# of Cancer types",
    low = CPCOLS[1],
    mid = CPCOLS[2],
    high = CPCOLS[3]
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  theme_bw()

# save the result ---------------------------------------------------------

save.image(file.path(expr_path_a, ".rda_04_gsva_other_pathway.rda"))
