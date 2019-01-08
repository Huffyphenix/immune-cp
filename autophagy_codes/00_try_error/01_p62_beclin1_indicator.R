
# library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)


# path --------------------------------------------------------------------

tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path <- file.path(expr_path, "03_a_gene_expr")
rppa_path <- "/extraspace/liucj/projects/6.autophagy/00_try_error/01_rppa_try"

# load data ---------------------------------------------------------------
rppa_expr <- readr::read_rds(file.path(tcga_path, "pancan32_rppa_expr.rds.gz"))
rppa_name <- readr::read_rds(file.path(tcga_path, "rppa_name_symbol.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path, "rds_03_a_atg_lys_gene_list.rds.gz"))
rppa_name %>% dplyr::semi_join(gene_list, by = "symbol") -> atg_rppa

atg_rppa %>% 
  dplyr::inner_join(gene_list, by = "symbol") %>% 
  dplyr::select(1,2, process) %>% 
  dplyr::arrange(process) -> sym_func
knitr::kable(sym_func)


# atg involved rppa ---------------------------------------------------------
fn_select_marker <- function(.x, sym_func){
  # rppa_expr$protein_expr[[1]] -> .x
  sym_func %>% dplyr::inner_join(.x, by = c("protein", "symbol"))
}

rppa_expr %>% 
  dplyr::mutate(protein_expr = purrr::map(.x = protein_expr, .f = fn_select_marker, sym_func)) -> atg_rppa_expr

# p62 vs. beclin1 ratio test---------------------------------------------------------------------

fn_class <- function(.x, .y, .gene){
  # atg_rppa_expr$protein_expr[[3]] -> .x
  # .gene <- c("SQSTM1", "BECN1")
  print(.y)
 
  .x %>% 
    dplyr::filter(symbol %in% .gene) %>%
    dplyr::select(-protein, -process) %>% 
    tidyr::gather(key = barcode, value = rppa, -symbol) %>% 
    tidyr::spread(key = symbol, value = rppa) %>% 
    tidyr::drop_na() %>% 
    dplyr::mutate(beclin1 = exp(BECN1), p62 = exp(SQSTM1)) %>% 
    dplyr::mutate(ratio = p62 / beclin1) -> .d
  if ( nrow(.d) < 10 ) return(NULL)
  
  dplyr::bind_rows(
    cor.test(~beclin1 + p62, data = .d, method = "spearman") %>% 
      broom::tidy() %>% 
      dplyr::select(coef = estimate, pval = p.value) %>% 
      tibble::add_column(compare = "beclin1 vs. p62", .before = 1),
    cor.test(~beclin1 + ratio, data = .d, method = "spearman") %>% 
      broom::tidy() %>% 
      dplyr::select(coef = estimate, pval = p.value) %>% 
      tibble::add_column(compare = "beclin1 vs. ratio", .before = 1),
    cor.test(~p62 + ratio, data = .d, method = "spearman") %>% 
      broom::tidy() %>% 
      dplyr::select(coef = estimate, pval = p.value) %>% 
      tibble::add_column(compare = "p62 vs. ratio", .before = 1)
  ) %>% 
    tibble::as_tibble()
  # k = 2
  # x <- seq(range(.d$BECN1)[1], range(.d$BECN1)[2], length.out = k)
  # y <- seq(range(.d$SQSTM1)[1], range(.d$SQSTM1)[2], length.out = k)
  # z <- seq(range(.d$ratio)[1], range(.d$ratio)[2], length.out = k)
  # 
  # fn_fz <- function(.zz, .xx, .yy, .d){
  #   .d %>% 
  #     dplyr::mutate(beclin1 = ifelse(BECN1 > .xx, "H", "L")) %>% 
  #     dplyr::mutate(p62 = ifelse(SQSTM1 > .yy, "L", "H")) %>% 
  #     dplyr::mutate(rat = ifelse(ratio > .zz, "L", "H"))  %>% 
  #     dplyr::mutate(bp = beclin1 == p62) %>% 
  #     dplyr::mutate(br = beclin1 == rat) %>% 
  #     dplyr::mutate(pr = p62 == rat) %>% 
  #     dplyr::mutate(bpr = dplyr::case_when(
  #       beclin1 == "H" & p62 == "H" & rat == "H" ~ TRUE,
  #       beclin1 == "L" & p62 == "L" & rat == "L" ~ TRUE,
  #       TRUE ~ FALSE
  #     )) %>% 
  #     dplyr::summarise(bpp = sum(bp) / n(), brp = sum(br) / n(), prp = sum(pr) / n(), bprp = sum(bpr) / n()) %>% 
  #     unlist() %>% 
  #     sum() -> total_percent
  #     tibble::tibble(.xx, .yy, .zz, total_percent)
  # }
  # 
  # fn_fy <- function(.yy, .xx, .zz){
  #   .zz %>% purrr::map(.f = fn_fz, .xx, .yy, .d)
  # }
  # 
  # fn_fx <- function(.xx, .yy, .zz){
  #   .yy %>% purrr::map(.f = fn_fy, .xx, .zz)
  # }
  # 
  # x %>% 
  #   purrr::map(.f = fn_fx, y, z) %>% 
  #   unlist(recursive = F) %>% 
  #   unlist(recursive = F) %>% 
  #   dplyr::bind_rows() %>% 
  #   dplyr::arrange(-total_percent) %>% 
  #   dplyr::slice(1) %>% 
  #   dplyr::select(-total_percent) %>% 
  #   unlist(use.names = F) -> params
  # 
  # .d %>% 
  #   dplyr::mutate(beclin1 = ifelse(BECN1 > params[1], "H", "L")) %>% 
  #   dplyr::mutate(p62 = ifelse(SQSTM1 > params[2], "L", "H")) %>% 
  #   dplyr::mutate(rat = ifelse(ratio > params[3], "L", "H"))  %>% 
  #   dplyr::mutate(bp = beclin1 == p62) %>% 
  #   dplyr::mutate(br = beclin1 == rat) %>% 
  #   dplyr::mutate(pr = p62 == rat) %>% 
  #   dplyr::mutate(bpr = dplyr::case_when(
  #     beclin1 == "H" & p62 == "H" & rat == "H" ~ TRUE,
  #     beclin1 == "L" & p62 == "L" & rat == "L" ~ TRUE,
  #     TRUE ~ FALSE
  #   )) %>% 
  #   dplyr::summarise(bpp = sum(bp) / n(), brp = sum(br) / n(), prp = sum(pr) / n(), bprp = sum(bpr) / n()) %>%
  #   unlist()

}

atg_rppa_expr %>% 
  dplyr::mutate(class = purrr::map2(.x = protein_expr, .y = cancer_types, .f = fn_class, c("SQSTM1", "BECN1"))) %>% 
  dplyr::filter(!purrr::map_lgl(.x = class, .f = is.null)) %>% 
  tidyr::unnest(class) -> corr_atg

corr_atg %>% 
  dplyr::filter(pval < 0.1, abs(coef) > 0.3) %>% 
  ggplot(aes(x = compare, y = cancer_types, size = -log10(pval), color = coef)) +
  geom_point() +
  scale_size(name = "P-value") +
  scale_color_gradient2(name = "Coef", low = "blue", mid = "white", high = "red") +
  theme_bw() +
  labs(x = "Correlation of three markers", y = "Cancer Types", 
       title = "Protein marker selection", 
       subtitle = "The ratio is p62 / beclin1 of protein levels. \nIn general, p62 should be negative correlates with beclin1.") -> p
ggsave(filename = "Correlation_of_markers.pdf", plot = p, device = "pdf", path = rppa_path, width = 5, height = 6)

# Comment 1 ---------------------------------------------------------------

# Is it reasonable to use beclin1 and ratio of p62 / beclin1 as autophagy flux marker.

# use becn1 and ratio to classification -----------------------------------
fn_intersection <- function(.x, .y){
  # .x <- atg_rppa_expr$protein_expr[[3]]
  .x %>% 
    dplyr::filter(symbol %in% .gene) %>%
    dplyr::select(-protein, -process) %>% 
    tidyr::gather(key = barcode, value = rppa, -symbol) %>% 
    tidyr::spread(key = symbol, value = rppa) %>% 
    tidyr::drop_na() %>% 
    dplyr::mutate(beclin1 = exp(BECN1), p62 = exp(SQSTM1)) %>% 
    dplyr::mutate(ratio =  p62 / beclin1) %>% 
    dplyr::mutate(bg = ifelse(beclin1 > median(beclin1), "H", "L")) %>% 
    dplyr::mutate(pg = ifelse(p62 < median(p62), "H", "L")) %>% 
    dplyr::mutate(rg = ifelse(ratio < median(ratio), "H", "L")) %>% 
    dplyr::mutate(bp = bg == pg) %>% 
    dplyr::mutate(br = bg == rg) %>%
    dplyr::mutate(pr = pg == rg) %>%
    dplyr::mutate(bpr = dplyr::case_when(
          bg == "H" & pg == "H" & rg == "H" ~ TRUE,
          bg == "L" & pg == "L" & rg == "L" ~ TRUE,
          TRUE ~ FALSE
        )) -> .d
  if ( nrow(.d) < 10 ) return(NULL)
  
  tibble::tibble(
    `beclin1` = nrow(.d),# - sum(.d$bp) - sum(.d$br) + sum(.d$bpr), # beclin1
    `p62` = nrow(.d),# - sum(.d$bp) - sum(.d$pr) +  sum(.d$bpr), # p62
    `ratio` = nrow(.d),# - sum(.d$pr) - sum(.d$br) + sum(.d$bpr), # ratio
    `beclin1 & p62` = sum(.d$bp), #- sum(.d$bpr), # beclin1 & p62
    `p62 & ratio` = sum(.d$pr), #- sum(.d$bpr), # p62 & ratio
    `beclin1 & ratio` = sum(.d$br), #- sum(.d$bpr), # beclin1 & ratio
    `three inter` = sum(.d$bpr) # three inter
  ) -> .inter
  

}

fn_draw_pie <- function(.inter){
  .inter %>% 
    unlist() %>% 
    setNames(c("area1", "area2", "area3", "n12", "n23", "n13", "n123")) %>% 
    as.list() -> .inter_list
  
  .inter_list$catergory <- c("beclin1", "p62", "ratio")
  .inter_list$lty <- "blank"
  .inter_list$fill <- c("skyblue", "pink1", "mediumorchid")

  purrr::lift(VennDiagram::draw.triple.venn)(.inter_list) -> p
}

atg_rppa_expr %>% 
  dplyr::mutate(inter = purrr::map(.x = protein_expr, .f = fn_intersection)) %>% 
  dplyr::filter(!purrr::map_lgl(.x = inter, .f = is.null)) %>% 
  dplyr::mutate(p = purrr::map(.x = inter, .f = fn_draw_pie)) -> atg_rppa_expr_plot

atg_rppa_expr_plot %>% 
  tidyr::unnest(inter) %>% 
  dplyr::mutate(
    `beclin1 & p62` = `beclin1 & p62` / beclin1,
    `p62 & ratio` =  `p62 & ratio` / beclin1, 
    `beclin1 & ratio` = `beclin1 & ratio` / beclin1, 
    `three inter` = `three inter` / beclin1) %>% 
  dplyr::select(cancer_types, `beclin1 & p62`, `p62 & ratio`, `beclin1 & ratio`, `three inter`) %>% 
  dplyr::arrange(`p62 & ratio`) %>% 
  knitr::kable()
















