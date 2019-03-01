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

xCell_TCGA_RSEM <- readr::read_tsv(file.path(xcell_path,"xCell_TCGA_RSEM.txt")) %>%
  tidyr::gather(-X1,key="barcode",value="score")

xCell_immune_type <- readr::read_tsv(file.path(xcell_path,"cell_type_classification.13059_2017_1349_MOESM1_ESM.txt")) %>%
  dplyr::select(`Full name`,`Cell types`,Subgroup) %>%
  dplyr::filter(Subgroup %in% c("Lymphoid","Myeloid"))
xCell_cell_type <- readr::read_tsv(file.path(xcell_path,"cell_type_classification.13059_2017_1349_MOESM1_ESM.txt")) %>%
  dplyr::select(`Full name`,`Cell types`,Subgroup) 

xCell_Fantom_cellname_adjust <- readr::read_tsv(file.path(xcell_path,"xCell_Fantom5_cellname_uniform.txt")) %>%
  dplyr::rename("Cell types" = "xCell_Cell types") %>%
  dplyr::inner_join(xCell_cell_type,by="Cell types") %>%
  dplyr::mutate(Type = ifelse(Subgroup == "Lymphoid" | Subgroup == "Myeloid","Immune Cell","Stromal Cell")) %>%
  dplyr::mutate(Type = ifelse(is.na(Fantom5_cellname),"Other stromal",Type)) %>%
  dplyr::select(`Cell types`,Type) %>%
  unique()

TCGA_purity <- readr::read_tsv(file.path(purity_path,"ncomms9971-s2.txt")) %>%
  dplyr::select(`Cancer type`,`Sample ID`,`CPE`) %>%
  dplyr::mutate(`Sample ID`=substr(`Sample ID`,1,15))

# classification of tcga cancers
TCGA_tissue <- readr::read_tsv("/project/huff/huff/data/TCGA/TCGA_cancer_tissue_classification.txt")

# data processing ---------------------------------------------------------

xCell_TCGA_RSEM %>%
  dplyr::rename("Cell types"= "X1") %>%
  dplyr::inner_join(xCell_Fantom_cellname_adjust,by="Cell types") %>%
  dplyr::rename("Group"= "Type") %>%
  dplyr::group_by(barcode,Group) %>%
  dplyr::mutate(group_Score = sum(score)) %>%
  dplyr::select(barcode,group_Score,Group) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key=Group,value=group_Score) -> xCell_TCGA_RSEM.immune_stroma.score
 
xCell_TCGA_RSEM.immune_stroma.score %>%
  dplyr::mutate(`Sample ID`=gsub(pattern = "\\.",replacement = "-",x = barcode)) %>%
  dplyr::inner_join(TCGA_purity,by="Sample ID") %>%
  dplyr::mutate(non_tumor_ratio = 1-CPE) %>%
  dplyr::mutate(Immune_ratio=non_tumor_ratio*`Immune Cell`/(`Stromal Cell`+`Immune Cell`)) %>%
  dplyr::mutate(Stromal_ratio=non_tumor_ratio*`Stromal Cell`/(`Stromal Cell`+`Immune Cell`)) -> xCell_TCGA_RSEM.immune_stroma.ratio

xCell_TCGA_RSEM.immune_stroma.ratio %>%
  readr::write_rds(file.path(immune_path,"genelist_data","xCell_TCGA_RSEM.immune_stroma.ratio.rds.gz"))

# TCGA samples dont included in cCell data --------------------------------
# library(xCell)
# exprMatrix = read.table(expr_file,header=TRUE,row.names=1, as.is=TRUE)
# xCellAnalysis(exprMatrix)

### Importantly, xCell performs best with heterogenous dataset. Thus it is recommended to use all data combined in one run, and not break down to pieces (especially not cases and control in different runs).
### So there is no need to run samples not in the xCell data set.
# TCGA_purity %>%
#   dplyr::filter(!`Sample ID` %in% xCell_TCGA_RSEM.immune_stroma.ratio$`Sample ID`) -> TCGA_sample_with_no_xCell_data
# 
# all_TCGA_exp <- readr::read_rds(file.path(tcga_path,"pancan33_expr.rds.gz"))
# 
# fn_get_subset_samples <- function(.x){
#   .x %>%
#     dplyr::select(-entrez_id) %>%
#     dplyr::filter(!symbol =="?") %>%
#     tidyr::gather(-symbol,key="barcode",value="exp_tmp") %>%
#     dplyr::mutate(barcode = substr(barcode,1,15)) %>%
#     dplyr::group_by(symbol,barcode) %>%
#     dplyr::mutate(exp = mean(exp_tmp)) %>%
#     dplyr::select(symbol,barcode,exp) %>%
#     unique() %>%
#     dplyr::filter(barcode %in% TCGA_sample_with_no_xCell_data$`Sample ID`)
# }
# all_TCGA_exp %>%
#   head(1) %>%
#   dplyr::mutate(purrr::map(expr,fn_get_subset_samples)) %>%
#   dplyr::select(-expr,-cancer_types) %>%
#   tidyr::unnest() %>%
#   tidyr::spread(key=barcode,value=exp) -> noxCell_smaple_exp.matrix

