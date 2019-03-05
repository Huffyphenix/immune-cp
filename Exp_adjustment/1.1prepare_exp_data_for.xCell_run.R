# Get immune fraction from immune + stromal data ---------------------------------

library(magrittr)

# data path ---------------------------------------------------------------
xcell_path <- "/project/huff/huff/data/TCGA/immune_infiltration/xCell"
purity_path <- "/project/huff/huff/data/TCGA_tumor_purity_from_ncomms9971"
immune_path <- "/project/huff/huff/immune_checkpoint"
tcga_path <- file.path(immune_path,"/data/TCGA_data")


# load data ---------------------------------------------------------------

TCGA_purity <- readr::read_tsv(file.path(purity_path,"ncomms9971-s2.txt")) %>%
  dplyr::select(`Cancer type`,`Sample ID`,`CPE`) %>%
  dplyr::mutate(`Sample ID`=substr(`Sample ID`,1,15))

all_TCGA_exp <- readr::read_rds(file.path(tcga_path,"pancan33_expr.rds.gz"))

### Importantly, xCell performs best with heterogenous dataset. Thus it is recommended to use all data combined in one run, and not break down to pieces (especially not cases and control in different runs).
### So there is no need to run samples not in the xCell data set.
TCGA_purity %>%
  dplyr::filter(!`Sample ID` %in% xCell_TCGA_RSEM.immune_stroma.ratio$`Sample ID`) -> TCGA_sample_with_no_xCell_data


fn_get_subset_samples <- function(.cancer,.x){
  TCGA_purity %>%
    dplyr::filter(`Cancer type` %in% .cancer) %>%
    dplyr::select(`Sample ID`) -> .samples
  .x %>%
    dplyr::filter(!symbol =="?") %>%
    dplyr::select(-entrez_id) %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate_all(mean)%>%  # rows with same gene symbols will be averaged 
    unique() %>%
    t()  %>%
    as.data.frame()-> .x.t
  .x.t %>%
    dplyr::mutate(barcode = substr(rownames(.x.t),1,15)) %>%
    dplyr::filter(barcode %in% .samples$`Sample ID`) %>%
    dplyr::select(barcode,everything()) %>%
    t() %>%
    as.data.frame()-> .x.fiter
  
  colnames(.x.fiter) <- t(.x.fiter["barcode",]) 
  .x.fiter <- .x.fiter[-1,]
  rownames(.x.fiter) <- t(.x.t["symbol",])
  print(.cancer)
  .x.fiter
}


all_TCGA_exp %>%
  dplyr::filter(cancer_types %in% unique(TCGA_purity$`Cancer type`)) %>%
  dplyr::mutate(expr_filter=purrr::map2(cancer_types,expr,.f=fn_get_subset_samples)) %>%
  dplyr::select(-expr) -> smaple_exp.matrix

