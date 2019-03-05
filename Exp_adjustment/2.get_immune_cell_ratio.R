## Data source : xCell: digitally portraying the tissue cellular heterogeneity landscape
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5688663/
## /project/huff/huff/data/TCGA/immune_infiltration/xCell/xCell_TCGA_RSEM.txt

# Get immune fraction from immune + stromal data ---------------------------------

library(magrittr)

# data path ---------------------------------------------------------------
xcell_path <- "/project/huff/huff/data/TCGA/immune_infiltration/xCell"
purity_path <- "/project/huff/huff/data/TCGA_tumor_purity_from_ncomms9971"
immune_path <- "/project/huff/huff/immune_checkpoint"
tcga_path <- file.path(immune_path,"/data/TCGA_data")

# load data ---------------------------------------------------------------

xCell_TCGA_RSEM <- readr::read_rds(file.path(xcell_path,"Pan21_xCell_results.rds.gz")) 
xCell_TCGA_RSEM %>%
  dplyr::mutate(xCell_gather = purrr::map(xCell,.f=function(.x){
    .x %>%
      as.data.frame() %>%
      dplyr::mutate(`Cell types` = rownames(.x)) %>%
      tidyr::gather(-`Cell types`,key="barcode",value="score")
  })) %>%
  dplyr::select(-xCell) %>%
  tidyr::unnest() -> xCell_TCGA_RSEM.gather

xCell_immune_type <- readr::read_tsv(file.path(xcell_path,"cell_type_classification.13059_2017_1349_MOESM1_ESM.txt")) %>%
  dplyr::select(`Full name`,`Cell types`,Subgroup) %>%
  dplyr::filter(Subgroup %in% c("Lymphoid","Myeloid"))
xCell_cell_type <- readr::read_tsv(file.path(xcell_path,"cell_type_classification.13059_2017_1349_MOESM1_ESM.txt")) %>%
  dplyr::select(`Full name`,`Cell types`,Subgroup) 

xCell_Fantom_cellname_adjust <- readr::read_tsv(file.path(xcell_path,"xCell_Fantom5_cellname_uniform.txt")) %>%
  dplyr::rename("Cell types" = "xCell_Cell types") %>%
  dplyr::inner_join(xCell_cell_type,by="Cell types") %>%
  dplyr::mutate(Type = ifelse(Subgroup == "Lymphoid" | Subgroup == "Myeloid","Immune Cell","Stromal Cell")) %>%  # fantom5 and xCell overlapped cells
  dplyr::mutate(Type = ifelse(is.na(Fantom5_cellname),"Other stromal",Type)) %>% # cell in xCell but not in fantom5 are calss into other stromal.
  dplyr::select(`Cell types`,Type) %>%
  unique()

TCGA_purity <- readr::read_tsv(file.path(purity_path,"ncomms9971-s2.txt")) %>%
  dplyr::select(`Cancer type`,`Sample ID`,`CPE`) %>%
  dplyr::mutate(`Sample ID`=substr(`Sample ID`,1,15))

# classification of tcga cancers
TCGA_tissue <- readr::read_tsv("/project/huff/huff/data/TCGA/TCGA_cancer_tissue_classification.txt")

# data processing ---------------------------------------------------------

xCell_TCGA_RSEM.gather %>%
  dplyr::inner_join(xCell_Fantom_cellname_adjust,by="Cell types") %>%
  dplyr::rename("Group"= "Type") %>%
  dplyr::group_by(barcode,Group) %>%
  dplyr::mutate(group_Score = sum(score)) %>%
  dplyr::select(barcode,group_Score,Group) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key=Group,value=group_Score) -> xCell_TCGA_RSEM.immune_stroma.score
 

# get ratio of each cell type ---------------------------------------------

xCell_TCGA_RSEM.immune_stroma.score %>%
  dplyr::mutate(`Sample ID`=gsub(pattern = "\\.",replacement = "-",x = barcode)) %>%
  dplyr::inner_join(TCGA_purity,by="Sample ID") %>%
  dplyr::mutate(non_tumor_ratio = 1-CPE) %>%
  dplyr::mutate(Immune_ratio=non_tumor_ratio*`Immune Cell`/(`Stromal Cell`+`Immune Cell`+`Other stromal`)) %>%
  dplyr::mutate(Stromal_ratio=non_tumor_ratio*`Stromal Cell`/(`Stromal Cell`+`Immune Cell`+`Other stromal`)) %>%
  dplyr::mutate(other_Stromal_ratio=non_tumor_ratio*`Other stromal`/(`Stromal Cell`+`Immune Cell`+`Other stromal`)) -> xCell_TCGA_RSEM.immune_stroma.ratio

xCell_TCGA_RSEM.immune_stroma.ratio %>%
  readr::write_rds(file.path(immune_path,"genelist_data","Pancan21.tcga.xCell.immune_stroma.ratio.rds.gz"))



