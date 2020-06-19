#######################
# tcga cancer, RNAseq exp, sample size statistic
#######################
library(magrittr)

# processed path
basic_path <- "/home/huff/project"
expr_path <- file.path(basic_path,"immune_checkpoint/result_20171025/expr_rds")
out_path <- file.path(basic_path,"immune_checkpoint/result_20171025")

gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))

# calculation -------------------------------------------------------------

gene_list_expr %>%
  dplyr::mutate(data = purrr::map(filter_expr,.f=function(.x){
    .x %>%
      tidyr::gather(-symbol,-entrez_id,key="barcode",value="exp") %>%
      dplyr::select(barcode) %>%
      unique()
  })) %>%
  dplyr::select(-filter_expr) %>%
  tidyr::unnest() -> TCGA_rnaseq_sample

TCGA_rnaseq_sample %>% 
  dplyr::mutate(patient = substr(x = barcode, start = 1, stop = 12)) %>% 
  dplyr::mutate(code = substr(x = barcode, start = 14, stop = 15)) %>% 
  dplyr::group_by(cancer_types,code) %>%
  dplyr::mutate(n = dplyr::n()) %>% 
  dplyr::select(cancer_types,code,n) %>%
  unique() %>%
  tidyr::spread(key="code",value="n") %>%
  tidyr::gather(-cancer_types,key="code",value="n") %>%
  dplyr::mutate(n = ifelse(is.na(n),0,n))%>%
  tidyr::spread(key="code",value="n") %>%
  dplyr::ungroup() -> TCGA_rnaseq_sample_statis

TCGA_rnaseq_sample %>% 
  dplyr::group_by(cancer_types) %>% 
  tidyr::nest() %>% 
  # head() %>%
  dplyr::mutate(data = purrr::map(.x = data, .y=cancer_types,.f = function(.x,.y) {
    print(.y)
    .x %>% 
      unique() %>% 
      dplyr::mutate(patient = substr(x = barcode, start = 1, stop = 12)) %>% 
      dplyr::mutate(code = substr(x = barcode, start = 14, stop = 15)) %>% 
      dplyr::filter(code %in% c('01', '11')) %>% 
      dplyr::group_by(code,patient) %>%
      tidyr::nest() %>%
      dplyr::mutate(id = purrr::map(data,.f=function(.y){
        if(nrow(.y)>0){
          .y %>%
            dplyr::mutate(id = 1:nrow(.y))
        } else{
          .y %>%
            dplyr::mutate(id = NA)
        }
      })) %>%
      tidyr::unnest() %>%
      dplyr::select(patient,code,barcode,id) %>%
      tidyr::spread(key="code",value="barcode") -> .paired_sample
    if(!any(colnames(.paired_sample) =="11")){
      .paired_sample %>%
        dplyr::mutate(`11`=NA) -> .paired_sample
    }
    if(!any(colnames(.paired_sample) =="01")){
      .paired_sample %>%
        dplyr::mutate(`01`=NA) -> .paired_sample
    }
    .paired_sample %>%
      dplyr::filter(!is.na(`01`) & !is.na(`11`)) -> .paired_sample
    .paired_sample %>%
      dplyr::select(patient) %>%
      unique() %>%
      nrow() -> .n
    
    if (nrow(.paired_sample)==0) {
      tibble::tibble(patient="NA", `11`="NA", `01` ="NA",num_paired_TN = .n) -> .paired_sample
    } else{
      .paired_sample %>%
        dplyr::ungroup() %>%
        dplyr::select(-id) %>%
        dplyr::mutate(num_paired_TN = .n) -> .paired_sample
    }
    .paired_sample  %>%
      dplyr::rename("primary_solid_tumor" = "01",
                    "solid_normal_tissue" = "11")%>% 
      tidyr::nest(-num_paired_TN,.key="paired_samples_detials")
  })) %>% 
  tidyr::unnest(cols = data)  ->  TCGA_rnaseq_Paired_sample_statis  

TCGA_rnaseq_Paired_sample_statis %>%
  dplyr::inner_join(TCGA_rnaseq_sample_statis,by="cancer_types") %>%
  dplyr::rename("num_primary_solid_tumor" = "01",
                "num_recurrent_solid_tumor" = "02",
                "num_primary_blood_derived_cancer" = "03",
                "num_new_primary" = "05",
                "num_metastatic" = "06",
                "num_additional_metastatic" = "07",
                "num_solid_tissue_normal" = "11") -> TCGA_rnaseq_sample_overall



# output ------------------------------------------------------------------

TCGA_rnaseq_sample_overall %>%
  dplyr::select(-paired_samples_detials) %>%
  readr::write_tsv(file.path(out_path,"e_6_exp_profile","TCGA_rnaseq_sample_overall_statistic.tsv"))

TCGA_rnaseq_sample_overall %>%
  readr::write_rds(file.path(out_path,"e_6_exp_profile","TCGA_rnaseq_sample_overall_statistic_and_sample_details.rds"),compress = "gz")
