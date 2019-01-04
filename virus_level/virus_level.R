
# get virus of tcga samples -----------------------------------------------
library(magrittr)

# data path
data_path <- "/project/huff/huff/data/TCGA/virus_level"

# load data
sample_virus_level <- readr::read_tsv(file.path(data_path,"virus_levels_of_samples.txt")) %>%
  dplyr::rename("barcode" = "X1")

# data process
sample_virus_level %>%
  dplyr::filter(Provenance == "TCGA") %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  readr::write_rds(file.path(data_path,"TCGA_virus_level.rds.gz"),compress = "gz")
