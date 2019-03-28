library(magrittr)
library(methods)
library(magrittr)
library(CancerSubtypes)
library(ConsensusClusterPlus)
library(SNFtool)
### methylation data path
basic_path <- "/home/huff/project"
ICP_methy_path <- file.path(basic_path,"TCGA_survival/data/ICP_methy_pancan")

data_result_path <- file.path(basic_path,"immune_checkpoint/genelist_data")
### load methylation data
ICP_methy <- readr::read_rds(file.path(ICP_methy_path,"ICP_methy_450_pancan.rds.gz"))

### methylation data process
fn_data_process <- function(.x){
  .x %>%
    dplyr::select(-symbol) %>%
    tidyr::gather(-tag_symbol,key="barcode",value="methy")
}

ICP_methy %>%
  head(2) %>%
  dplyr::mutate(gather = purrr::map(methy,fn_data_process)) %>%
  dplyr::select(-methy) -> ICP_methy.gather

ICP_methy.gather %>%
  tidyr::unnest() %>%
  dplyr::filter(! is.na(methy)) %>%
  dplyr::filter(substr(barcode,14,14)!=1) -> ICP_methy.gather.unnest

# scaled by cancer types
ICP_methy.gather.unnest %>%
  dplyr::mutate(methy = as.numeric(methy)) %>%
  dplyr::group_by(cancer_types,tag_symbol) %>%
  dplyr::mutate(methy_scaled = scale(methy)) %>%
  dplyr::select(-methy) ->  ICP_methy.gather.unnest.scaled

ICP_methy.gather.unnest.scaled %>%
  dplyr::select(tag_symbol,barcode,methy) %>%
  dplyr::group_by(tag_symbol,barcode) %>%
  dplyr::mutate(methy = mean(methy)) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key = barcode,value = methy) -> PanCan_ICP_methy.spread

PanCan_ICP_methy.spread[1:10,1:2]

PanCan_ICP_methy_matrix <- as.matrix(PanCan_ICP_methy.spread[,-1])
rownames(PanCan_ICP_methy_matrix) <- PanCan_ICP_methy.spread$tag_symbol
PanCan_ICP_methy_matrix %>%
  readr::write_rds(file.path(data_result_path,"PanCan_ICP_methy_matrix.alltag-scaled.rds.gz"),compress = "gz")

# memary controling
rm(PanCan_ICP_methy.spread)
rm(ICP_methy.gather.unnest.scaled)
rm(ICP_methy.gather.unnest)
rm(ICP_methy)

### do cluster analysis

# reduce the dataset to the top 5,000 most variable genes, measured by median absolute deviation(mad). 绝对中位差
# d <- readr::read_rds(file.path(data_result_path,"PanCan_ICP_methy_matrix.alltag.rds.gz"))
d <- PanCan_ICP_methy_matrix
rm(PanCan_ICP_methy_matrix)

index=which(is.na(d))
d=data.imputation(d,fun="median")

dc = sweep(d,1, apply(d,1,median,na.rm=T)) ## median center genes

dt = as.dist(1-cor(dc,method="pearson"))
results = ConsensusClusterPlus(dt,maxK=20,reps=100,pItem=0.8,pFeature=1,title="methy_CC",distance="pearson",clusterAlg="hc",seed=1262118388.71279)

results %>%
  readr::write_rds(file.path(data_result_path,"genelist_methy_alltag-scaled_CC_20.rds.gz"),compress = "gz")

