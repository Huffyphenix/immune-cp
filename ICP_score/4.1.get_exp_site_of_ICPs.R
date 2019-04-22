###############################
# class the ICP genes into 5 features by expression mode in immune and tumor cell
# [-∞, -5]：genes express only in tumor cell
# [-5, -2]: genes express mainly in tumor cell
# [-2, 2]: genes express in both tumor and immune cell
# [2, 5]: genes express mainly in immune cell
# [5, +∞]: genes express only in immune cell

library(magrittr)

# data path ---------------------------------
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"ICP_score/2.GSVA-ICPs_exp_site_5_feature")

# laod data ---------------------------------
ICPs_exp_in_TI <- readr::read_tsv(file.path(immune_res_path,"ICP_exp_patthern","ICP_exp_pattern_in_immune_tumor_cell.detailed_tissues.tsv")) %>%
  dplyr::rename("Tissues" = "Characteristics[Tissue]", "log2FC(I/T)" = "log2FC") %>%
  dplyr::select(symbol, Tissues, `log2FC(I/T)`)

TCGA_tissue <- readr::read_tsv("/home/huff/project/data/TCGA/TCGA_cancer_tissue_classification.txt")

ICPs_exp_in_TI %>%
  dplyr::inner_join(TCGA_tissue, by = "Tissues") -> ICPs_exp_in_TI.cancer


# classify the genes in each cancer ---------------------------------------
fn_group <- function(.x) {
  if (.x <= (-5)) {
    .group <- "IRs_express_only_in_tumor_cell"
  } else if (.x > (-5) & .x <= (-2)) {
    .group <- "IRs_express_mainly_in_tumor_cell"
  } else if (.x > (-2) & .x < 2) {
    .group <- "IRs_express_in_both_tumor_and_immune_cell"
  } else if (.x >=2 & .x < 5) {
    .group <- "IRs_express_mainly_in_immune_cell"
  } else if (.x >= 5) {
    .group <- "IRs_express_only_in_immune_cell"
  }
  .group
}

ICPs_exp_in_TI.cancer %>%
  dplyr::mutate(ICP_exp_group = purrr::map(`log2FC(I/T)`, fn_group)) %>%
  tidyr::unnest() -> ICPs_exp_in_TI.cancer.site_groups

ICPs_exp_in_TI.cancer.site_groups %>%
  readr::write_tsv(file.path(res_path, "ICPs_exp_in_TI.cancer.site_groups.tsv"))
