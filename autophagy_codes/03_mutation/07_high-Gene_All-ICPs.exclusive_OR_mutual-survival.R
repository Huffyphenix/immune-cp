############################################
# get mutual exclusive mutation profile between all ICPs.
############################################

library(magrittr)
library(ggplot2)
out_path <- "/home/huff/project/immune_checkpoint/result_20171025"
snv_path <- file.path(out_path, "m_2_snv")


# LOAD DATA 
ICP_highGene_snv <- readr::read_rds(file.path(snv_path,"ICP_highGene_snv.rds.gz")) 

survival_data <- readr::read_rds(file.path("/home/huff/project/data/TCGA-survival-time/cell.2018.survival","TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>%
  dplyr::rename("cancer_types" = "type")

ICP_highGene_snv %>%
  dplyr::inner_join(survival_data, by = "cancer_types") %>%
  dplyr::rename("survival" = "data") -> ICP_highGene_snv.survival

# GET EXCLUSIVE MUTATION PROFILE BETWEEN ICPS IN EACH CANCERS -----------------------
# SURVIVAL DIFFERENCE BETWEEN EACH MUTATION GROUPS
library(showtext)
font_add("Arial","ARIAL.TTF") # sudo进入container，将/home/huff/fonts/中的字体拷贝到/usr/share/fonts/,then do:fc-cache

library(export)
source("/home/huff/project/github/immune-cp/autophagy_codes/03_mutation/fn_get_mutation_exclusive.R")
source("/home/huff/project/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

library(survival)
fn_survival_on <- function(V1, V2, .data1, .data2, cancer_types, .survival){
  # .data <- filter_cnv
  # V1 <- 'TP53'
  # V2 <- 'HLA-E'
  print("start fn_survival_on --------------")
  V1 <- as.character(V1)
  V2 <- as.character(V2)
  # print(V1)
  #  print(V2)
  .data2 %>%
    dplyr::filter(symbol %in% c(V1)) %>% 
    dplyr::select(-symbol) %>%
    tidyr::gather(key = barcode, value = mut) %>%
    dplyr::mutate(mut = ifelse(is.na(mut),0,1)) %>%
    dplyr::mutate(mut = as.integer(mut)) %>%
    dplyr::inner_join(.data1, by = "barcode") %>%
    tidyr::unite(group, c("mut", V2)) %>%
    dplyr::inner_join(.survival, by = "barcode") -> plot_ready
  
  # do survival diff analysis, get pvalue 
  fit <- survfit(survival::Surv(time, status) ~ group, data = plot_ready, na.action = na.exclude)
  diff <- survdiff(survival::Surv(time, status) ~ group, data = plot_ready, na.action = na.exclude)
  kmp <- 1 - pchisq(diff$chisq, df = length(levels(as.factor(plot_ready$group))) - 1)
  
  # plot if pvalue is significant
  if (kmp <= 0.05) {
    title <- paste(cancer_types, paste(V1, V2, sep = "-"))
    color_list <- tibble::tibble(group = c( "0_0", "0_1", "1_0", "1_1"), 
                   color = c("#00B2EE", "#CDAD00", "pink1","#CD2626"))
    sur_name <- paste(cancer_types, V1, V2, sep = "_")
    result_path <- file.path(snv_path,"mutual_exclusive_allICPs-highGene","survival_5y")
    fn_survival(plot_ready,title,color_list,"group",sur_name,xlab = "Time (years)",result_path,3,4,lx = 0.8,ly = 0.8)
    print("draw survival plot ----------")
  }
  
  name <- paste(c(cancer_types, V1, V2), collapse = "_")
  print(name)
  
  print("end fn_survival_on --------------")
  
  tibble::tibble(te = name, kmp = signif(kmp))
}

fn_survival_pre <- function(cancer_types, ICP_SNV, highGene_SNV, survival, cluster){
  ICP_SNV %>%
    tidyr::gather(-symbol,key = "barcode",value = "mut") %>%
    dplyr::mutate(mut = ifelse(is.na(mut),0,mut)) %>%
    dplyr::group_by(barcode) %>%
    dplyr::mutate(mut_overall = sum(mut)) %>%
    dplyr::mutate(mut_overall = ifelse(mut_overall > 0,1,mut_overall)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-symbol,-mut) %>%
    dplyr::rename("ICPs" = "mut_overall") %>%
    unique() -> ICP_overall_mut
  
  .survival <- survival %>%
    dplyr::select(bcr_patient_barcode, PFI.1, PFI.time.1) %>%
    dplyr::rename("barcode" = "bcr_patient_barcode","status" = "PFI.1","time" = "PFI.time.1") %>%
    dplyr::filter(time <= 1825) %>%
    dplyr::mutate(time = time/365)
  
  expand.grid(highGene_SNV$symbol,"ICPs") -> .gene_pairs
  
  .gene_pairs %>% 
    dplyr::mutate(rs = purrr::map2(Var1, Var2, .f = fn_survival_on, .data1 = ICP_overall_mut, .data2 = highGene_SNV,cancer_types = cancer_types, .survival = .survival)) %>% 
    dplyr::as_tibble() %>%
    dplyr::ungroup() %>%
    dplyr::select(rs) %>% 
    tidyr::unnest() -> .gene_pairs_pval
  
  .gene_pairs_pval 
}

# RUNNING
ICP_highGene_snv.survival %>% 
  purrr::pmap(.f = fn_survival_pre, cluster = cluster) %>% 
  dplyr::bind_rows() %>%
  tidyr::separate(col = te, into = c('cancer_types', 'g1', 'g2')) -> mutual_exclusive_survival

mutual_exclusive_survival %>%
  readr::write_tsv(file.path(snv_path,"mutual_exclusive_allICPs-highGene", "survival_5y", "allICPs_mutual_exclusive_highGeneSNV-survival.tsv"))

