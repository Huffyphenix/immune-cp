## Data source : xCell: digitally portraying the tissue cellular heterogeneity landscape
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5688663/
## /project/huff/huff/data/TCGA/immune_infiltration/xCell/xCell_TCGA_RSEM.txt

# Get immune fraction from immune + stromal data ---------------------------------

library(magrittr)

# data path ---------------------------------------------------------------
xcell_path <- "/project/huff/huff/data/TCGA/immune_infiltration/xCell"
purity_path <- "/project/huff/huff/data/TCGA_tumor_purity_from_ncomms9971"

# load data ---------------------------------------------------------------

xCell_TCGA_RSEM <- readr::read_tsv(file.path(xcell_path,"xCell_TCGA_RSEM.txt")) %>%
  tidyr::gather(-X1,key="barcode",value="score")

xCell_immune_type <- readr::read_tsv(file.path(xcell_path,"cell_type_classification.13059_2017_1349_MOESM1_ESM.txt")) %>%
  dplyr::select(`Full name`,`Cell types`,Subgroup) %>%
  dplyr::filter(Subgroup %in% c("Lymphocytes","Myeloid"))

TCGA_purity <- readr::read_tsv(file.path(purity_path,"ncomms9971-s2.txt")) %>%
  dplyr::select(`Sample ID`,`CPE`) %>%
  dplyr::mutate(`Sample ID`=substr(`Sample ID`,1,15))


# data processing ---------------------------------------------------------

xCell_TCGA_RSEM %>%
  dplyr::mutate(Group = ifelse(X1 %in% xCell_immune_type$`Cell types`,"Immune cell","Other")) %>%
  dplyr::group_by(barcode,Group) %>%
  dplyr::mutate(group_Score = sum(score)) %>%
  dplyr::select(barcode,group_Score,Group) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key=Group,value=group_Score) -> xCell_TCGA_RSEM.immunecell_other
 
xCell_TCGA_RSEM.immunecell_other %>%
  dplyr::mutate(`Sample ID`=gsub(pattern = "\\.",replacement = "-",x = barcode)) %>%
  dplyr::inner_join(TCGA_purity,by="Sample ID") %>%
  dplyr::mutate(non_tumor = 1-CPE) %>%
  dplyr::mutate(Immune=non_tumor*`Immune cell`/(Other+`Immune cell`)) -> xCell_TCGA_RSEM.immunecell_ratio


# TCGA samples dont included in cCell data --------------------------------
library(xCell)
exprMatrix = read.table(expr_file,header=TRUE,row.names=1, as.is=TRUE)
xCellAnalysis(exprMatrix)

