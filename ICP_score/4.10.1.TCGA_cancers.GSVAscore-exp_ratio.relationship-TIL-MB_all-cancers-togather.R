### mean expression of each gene features in TCGA
library(magrittr)
library(tidyverse)

# data path ---------------------------------------------------------------
# server 1
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path(basic_path,"immune_checkpoint/data/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")

res_path <- file.path(immune_res_path,"TCGA_GSVAScore/GSVA_add_exp_ratio")

exp_ratio_path <- file.path(immune_res_path,"ICP_score.new")
GSVA_path <- file.path(immune_res_path,"TCGA_GSVAScore")
tcga_path <- file.path(basic_path,"/data/TCGA")

# load image --------------------------------------------------------------
# load(file.path(res_path,"TCGA.GSVA_score.ICP_features-all_cancer_togather.rda"))

# load data ---------------------------------------------------------------
GSVA.score.alltumor <- 
  readr::read_rds(file.path(GSVA_path,"TCGA.All-cancer-Tumor_GSVA.score_ICPs_features.rds.gz")) %>%
  tidyr::unnest() %>%
  tidyr::gather(-cancer_types,-barcode,key="Features",value="value") %>%
  dplyr::mutate(Features=gsub(" ","_",Features)) %>%
  tidyr::spread(key="Features",value="value")
exp_ratio <- 
  readr::read_rds(file.path(exp_ratio_path,"TCGA_mean_fold_ratio_features_value.rds.gz"))  %>%
  tidyr::gather(-barcode,key="Features",value="value") %>%
  # dplyr::mutate(Features=gsub("\\(","",Features)) %>%
  # dplyr::mutate(Features=gsub("\\)","",Features)) %>%
  # dplyr::mutate(Features=gsub("/","_",Features)) %>%
  dplyr::mutate(Features=gsub(" ","_",Features)) %>%
  dplyr::mutate(Features=gsub("-",".",Features)) %>%
  tidyr::spread(key="Features",value="value")


# combine data ------------------------------------------------------------
GSVA.score.alltumor %>%
  dplyr::inner_join(exp_ratio, by="barcode") %>%
  tidyr::nest(-cancer_types, .key="GSVA") -> GSVA_exp_ratio_combine

# laod ICP feature info ---------------------------------------------------

gene_list <- readr::read_tsv(file.path(gene_list_path, "ICPs_all_info_class.tsv")) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Only_exp_on_Immune","Mainly_exp_on_Immune"),"Mainly_exp_on_Immune",Exp_site)) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Only_exp_on_Tumor","Mainly_exp_on_Tumor" ),"Mainly_exp_on_Tumor",Exp_site.1)) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("N"),"Not_sure",Exp_site.1))

# plot theme --------------------------------------------------------------

my_theme <-   theme(
  panel.background = element_rect(fill = "white",colour = "black"),
  panel.grid.major=element_line(colour=NA),
  axis.text.y = element_text(size = 10,colour = "black"),
  axis.text.x = element_text(size = 10,colour = "black"),
  # legend.position = "none",
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 12),
  legend.background = element_blank(),
  legend.key = element_rect(fill = "white", colour = "black"),
  plot.title = element_text(size = 20),
  axis.text = element_text(colour = "black"),
  strip.background = element_rect(fill = "white",colour = "black"),
  strip.text = element_text(size = 10),
  text = element_text(color = "black")
)
# 1. GSVA score correlation with TIL -----------------------------------------
# 1.1. load TIL ----
immunity_path_2 <- "/project/huff/huff/immune_checkpoint/data/immunity"
TIMER_immunity <- readr::read_tsv(file.path(immunity_path_2,"immuneEstimation.txt")) %>%
  dplyr::mutate(TIL = B_cell+CD4_Tcell+CD8_Tcell+Neutrophil+Macrophage+Dendritic)

# 1.2.function --------
fn_correlation_and_DE <- function(GSVA,TIL,key){
  res <- tibble::tibble()
  for(cell in key){
    TIL %>%
      dplyr::rename("IL" = cell) -> .TIL
    GSVA %>%
      dplyr::mutate(barcode = substr(barcode,1,15)) %>%
      tidyr::gather(-barcode,key="Features",value="GSVA_score")%>%
      dplyr::inner_join(.TIL,by="barcode") %>%
      dplyr::group_by(Features) %>%
      dplyr::mutate(group = ifelse(GSVA_score >= quantile(GSVA_score,0.5,na.rm = TRUE),"GSVA_high","GSVA_low")) %>%
      tidyr::nest(-Features) -> .tmp
    
    # correlation
    .tmp %>%
      dplyr::mutate(cor = purrr::map(data,.f=function(.x){
        cor.test(.x$IL,.x$GSVA_score,method = "spearman") %>%
          broom::tidy()
      })) %>%
      dplyr::select(-data) -> cor.res
    
    # DE
    .tmp %>%
      dplyr::mutate(DE = purrr::map2(data,Features,.f=function(.x,.y){
        print(.y)
        .x %>%
          dplyr::filter(group == "GSVA_high") %>%
          dplyr::mutate(mean = mean(cell)) %>%
          .$IL %>%
          mean() -> mean_high
        .x %>%
          dplyr::filter(group == "GSVA_low") %>%
          .$IL %>%
          mean() -> mean_low
        broom::tidy(wilcox.test(IL  ~ group, data = .x, alternative = "two.sided")) %>%
          dplyr::mutate(mean_high = mean_high, mean_low = mean_low) %>%
          dplyr::mutate(`log2FC(High/Low)` = log2((mean_high)/(mean_low))) 
      })) %>%
      dplyr::select(-data)  -> DE.res
    
    cor.res %>%
      dplyr::inner_join(DE.res,by="Features") %>%
      dplyr::mutate(cell_type = cell)-> .tmp
    rbind(res,.tmp)->res
  }
  res
}

# 1.3.calculation ------
GSVA_exp_ratio_combine %>% 
  # head(1) %>%
  dplyr::mutate(res = purrr::map(GSVA,.f=fn_correlation_and_DE,
                                 TIL = TIMER_immunity,
                                 key=c("B_cell","CD4_Tcell","CD8_Tcell","Neutrophil","Macrophage","Dendritic","TIL"))) %>%
  dplyr::select(-GSVA) %>%
  tidyr::unnest() -> GSVA.TIL.res
GSVA.TIL.res %>%
  tidyr::unnest() %>%
  readr::write_tsv(file.path(res_path,"1.TIL_with_GSVA_score","1.spm-cor.DE.TIL-GSVAscore.allcancers.tsv"))

# 1.4.ploting ---------

# 1.4.3.mean TIL DE in cancer types  --------
GSVA.TIL.res %>%
  tidyr::unnest() %>%
  dplyr::rename("cor"="estimate") %>%
  dplyr::mutate(FC = 2^`log2FC(High/Low)`) %>%
  dplyr::select(Features,FC,cor,cell_type) %>%
  unique() -> GSVA.TIL.cor.DE.plot

gene_list %>%
  dplyr::mutate(All_gene = "All_gene") %>%
  dplyr::select(symbol,functionWithImmune,family,Recepter_pairs,Exp_site.1,All_gene) %>%
  tidyr::gather(-symbol,key="Groups",value="Features") %>%
  tidyr::nest(-Features) %>%
  dplyr::filter(!Features %in% "Other") %>%
  dplyr::filter(!is.na(Features)) %>%
  dplyr::mutate(Features = gsub(" ","_",Features)) %>%
  dplyr::mutate(Feature_symbol = purrr::map(data,.f=function(.x){
    tibble::tibble(
      Feature_symbol=paste(.x$symbol,collapse="/"),
      n=length(.x$symbol))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(label = ifelse(n<5,paste("GSVA.",Feature_symbol,sep=""),paste("GSVA.",Features,sep=""))) %>%
  dplyr::right_join(GSVA.TIL.cor.DE.plot,by="Features") %>%
  dplyr::mutate(label=ifelse(is.na(label),Features,paste("GSVA.",Features,sep=""))) -> label.data

library(ggrepel)

label.data %>%
  dplyr::mutate(label = ifelse(cor >=0.4 & FC >=1.5, label, NA)) %>%
  dplyr::mutate(label_or_not = ifelse(cor >=0.3 & FC >=1.5, "yes", "no")) %>%
  ggplot(aes(x=cor,y=FC)) +
  geom_jitter(aes(color=label_or_not)) +
  facet_wrap(".~cell_type",scales = "free") +
  # geom_text_repel(aes(x=cor,y=FC,label=label),color="red",size=3) +
  scale_color_manual(
    values = c("black","red")
  ) +
  labs(x="Spearman correlation between \nTIL and GSVA score in tumors",
       y="TIL fold change between high\nand low GSVA score groups in tumors") +
  my_theme +
  theme(legend.position = "none")+
  geom_hline(yintercept = 1.5,linetype = 2) +
  geom_vline(xintercept = 0.3,linetype = 2)


ggsave(file.path(res_path,"1.TIL_with_GSVA_score","3.DE.Cor.TIL-GSVAscore.allcancers.pdf"),device = "pdf",height = 9, width=12)
ggsave(file.path(res_path,"1.TIL_with_GSVA_score","3.DE.Cor.TIL-GSVAscore.allcancers.png"),device = "png",height = 9, width=12)


# 2. GSVA score with mutation burden ----------------------------
# 2.1. load mutation data ------
TCGA_mutation_burden <- readr::read_rds(file.path("/project/huff/huff/data/TCGA/snv","pancan33_maf_syn7824274.sample.mut_burden.rds.gz"))

# 2.2. calculation --------
GSVA_exp_ratio_combine %>% 
  # head(1) %>%
  dplyr::mutate(res = purrr::map(GSVA,.f=fn_correlation_and_DE,
                                 TIL = TCGA_mutation_burden %>% 
                                   dplyr::select(sample,mutation_durden) %>%
                                   dplyr::rename("barcode"="sample"),
                                 key="mutation_durden")) %>%
  dplyr::select(-GSVA) %>%
  tidyr::unnest() -> GSVA.MB.res

GSVA.MB.res %>%
  tidyr::unnest() %>%
  readr::write_tsv(file.path(res_path,"2.MB_with_GSVA_score","1.spm-cor.DE.MB-GSVAscore.allcancers.tsv"))

# 2.3.ploting ---------

# readr::read_tsv(file.path(tcga_path,"02_pcc.tsv")) %>%
#   dplyr::select(cancer_types,color) %>%
#   unique() -> cancer.color.MB

# 2.3.3.mean MB DE and Cor. in cancer types  --------
GSVA.MB.res %>%
  tidyr::unnest() %>%
  dplyr::group_by(Features) %>%
  dplyr::mutate(cor = estimate,FC = 2^`log2FC(High/Low)`) %>%
  dplyr::select(Features,FC,cor) %>%
  dplyr::ungroup() %>%
  unique() -> GSVA.MB.cor.DE.plot

gene_list %>%
  dplyr::mutate(All_gene = "All_gene") %>%
  dplyr::select(symbol,functionWithImmune,family,Recepter_pairs,Exp_site.1,All_gene) %>%
  tidyr::gather(-symbol,key="Groups",value="Features") %>%
  tidyr::nest(-Features) %>%
  dplyr::filter(! Features %in% "Other") %>%
  dplyr::filter(!is.na(Features)) %>%
  dplyr::mutate(Features = gsub(" ","_",Features)) %>%
  dplyr::mutate(Feature_symbol = purrr::map(data,.f=function(.x){
    tibble::tibble(
      Feature_symbol=paste(.x$symbol,collapse="/"),
      n=length(.x$symbol))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(label = ifelse(n<5,Feature_symbol,Features)) %>%
  dplyr::right_join(GSVA.MB.cor.DE.plot,by="Features") %>%
  dplyr::mutate(label=ifelse(is.na(label),Features,paste("GSVA.",Features,sep=""))) -> label.data.meanMB

library(ggrepel)

label.data.meanMB %>%
  dplyr::mutate(label = ifelse(cor >=0.3 & FC >=1.5, label, NA)) %>%
  dplyr::mutate(label_or_not = ifelse(cor >=0.3 & FC >=1.5, "yes", "no")) %>%
  ggplot(aes(x=cor,y=FC)) +
  geom_jitter(aes(color=label_or_not)) +
  geom_text_repel(aes(x=cor,y=FC,label=label),color="red",size=3) +
  scale_color_manual(
    values = c("black","red")
  ) +
  labs(x="Spearman correlation between\n
       mutation burden and GSVA score in tumors",
       y="Mutation burden fold change between high\n
       and low  GSVA score groups) in tumors") +
  my_theme +
  theme(legend.position = "none")+
  geom_hline(yintercept = 1.5,linetype = 2) +
  geom_vline(xintercept = 0.3,linetype = 2)

ggsave(file.path(res_path,"2.MB_with_GSVA_score","3.DE.Cor.MB-GSVAscore.allcancers.pdf"),device = "pdf",height = 4, width=5)
ggsave(file.path(res_path,"2.MB_with_GSVA_score","3.DE.Cor.MB-GSVAscore.allcancers.png"),device = "png",height = 4, width=5)



# 4. GSVA score with CTL --------------------------------------------------
# 4.1.load CTL data -----
CTL_path <- "/home/huff/project/data/TCGA/CTL_level_estimated"
CTL_data <- readr::read_rds(file.path(CTL_path,"CTL_estimated_from_CD8A_CD8B_GZMA_GZMB_PRF1.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::rename("barcode" = "sample","CTL"="CTL(average_exp)") %>%
  dplyr::select(cancer_types,barcode,CTL)

# 4.2.calculate correlation ----
GSVA_exp_ratio_combine %>% 
  # head(1) %>%
  dplyr::mutate(res = purrr::map(GSVA,.f=fn_correlation_and_DE,
                                 TIL = CTL_data %>% 
                                   dplyr::select(barcode,CTL) %>%
                                   dplyr::mutate(barcode = substr(barcode,1,15)) ,
                                 key="CTL")) %>%
  dplyr::select(-GSVA) %>%
  tidyr::unnest() -> GSVA.CTL.res

GSVA.CTL.res %>%
  tidyr::unnest() %>%
  readr::write_tsv(file.path(res_path,"4.CTL_with_GSVA_score","1.spm-cor.DE.CTL-GSVAscore.allcancers.tsv"))

# 4.3.all CTL DE and cor in all cancer  --------
GSVA.CTL.res %>%
  tidyr::unnest() %>%
  dplyr::rename("Correlation"="estimate") %>%
  dplyr::mutate(FC=2^`log2FC(High/Low)`)-> GSVA.CTL.DE.plot.all

gene_list %>%
  dplyr::mutate(All_gene = "All_gene") %>%
  dplyr::select(symbol,functionWithImmune,family,Recepter_pairs,Exp_site.1,All_gene) %>%
  tidyr::gather(-symbol,key="Groups",value="Features") %>%
  tidyr::nest(-Features) %>%
  dplyr::filter(! Features %in% "Other") %>%
  dplyr::filter(!is.na(Features)) %>%
  dplyr::mutate(Features = gsub(" ","_",Features)) %>%
  dplyr::mutate(Feature_symbol = purrr::map(data,.f=function(.x){
    tibble::tibble(
      Feature_symbol=paste(.x$symbol,collapse="/"),
      n=length(.x$symbol))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::right_join(GSVA.CTL.DE.plot.all,by="Features") %>%
  dplyr::mutate(n = ifelse(is.na(n),10,n)) %>%
  dplyr::mutate(label = ifelse(n<5,Feature_symbol,Features)) %>%
  dplyr::mutate(label=ifelse(is.na(Feature_symbol),Features,paste("GSVA.",Features,sep=""))) -> label.data.CTL

label.data.CTL %>%
  dplyr::mutate(label = ifelse(Correlation  >=0.5 & FC >=4, label, NA)) %>%
  dplyr::mutate(label_or_not = ifelse(Correlation  >=0.3 & FC >=1.5, "yes", "no")) %>%
  ggplot(aes(x=Correlation,y=FC)) +
  geom_jitter(aes(color=label_or_not)) +
  geom_text_repel(aes(x=Correlation,y=FC,label=label),color="red",size=3) +
  scale_color_manual(
    values = c("black","red")
  ) +
  labs(x="Spearman correlation\nbetween CTL and GSVA score in tumors",
       y="CTL fold change between high\nand low  GSVA score groups) in tumors") +
  my_theme +
  theme(legend.position = "none")+
  geom_hline(yintercept = 1.5,linetype = 2) +
  geom_vline(xintercept = 0.3,linetype = 2)

ggsave(file.path(res_path,"4.CTL_with_GSVA_score","4.DE.Cor.TIL-GSVAscore.allcancers.pdf"),device = "pdf",height = 4, width=5)
ggsave(file.path(res_path,"4.CTL_with_GSVA_score","4.DE.Cor.TIL-GSVAscore.allcancers.png"),device = "png",height = 4, width=5)

# 3.survival analysis -------------------------------------------------------
# 3.1.load survival data ----
clinical_tcga <- readr::read_rds(file.path(basic_path,"TCGA_survival/data","Pancan.Merge.clinical.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::select(-cancer_types) %>%
  unique() %>%
  dplyr::mutate(OS=as.numeric(OS),Status=as.numeric(Status),Age=as.numeric(Age)) %>%
  dplyr::group_by(barcode) %>%
  dplyr::mutate(OS= max(OS)) %>%
  dplyr::mutate(Status =  max(Status)) %>%
  dplyr::ungroup()

clinical <- readr::read_rds(file.path(basic_path,"data/TCGA-survival-time/cell.2018.survival","TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::select(type,bcr_patient_barcode,PFS,PFS.time) %>%
  tidyr::drop_na() %>%
  dplyr::rename("barcode" = "bcr_patient_barcode") %>%
  unique()

clinical %>%
  dplyr::inner_join(clinical_tcga,by="barcode") %>%
  dplyr::rename("cancer_types"= "type") -> clinical_all_data

# 3.2.function to do cox survival analysis  ----
# 3.2.1.univariable cox analysis ---------
fn_survival_test <- function(data,feature){
  print(feature)
  .cox <- survival::coxph(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)
  summary(.cox) -> .z
  
  # KM pvalue
  kmp <- 1 - pchisq(survival::survdiff(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)$chisq,df = length(levels(as.factor(data$group))) - 1)
  
  # mearge results
  tibble::tibble(
    group = rownames(.z$conf.int),
    n = .z$n,
    coef = .z$coefficients[,1],
    hr = .z$conf.int[1],
    hr_l = .z$conf.int[3],
    hr_h = .z$conf.int[4],
    coxp = .z$waldtest[3],
    kmp = kmp) %>%
    dplyr::mutate(status = ifelse(hr > 1, "High_risk", "Low_risk"))
  # multi-variable cox analysis
}

# 3.2.2.multi-variable cox analysis ---------
fn_survival_test.multiCox <- function(data,uni_sig_feature){
  # covariates <- uni_sig_feature$Features %>% unique()
  uni_sig_feature %>%
    # head() %>%
    dplyr::group_by(Features) %>%
    dplyr::mutate(res = purrr::map(Features,.f=function(.x){
      multi_formulas <- as.formula(paste('Surv(time, status)~', paste(c(.x,"Age", "Stage", "B_cell","CD4_Tcell","CD8_Tcell","Neutrophil","Macrophage","Dendritic","TIL","CTL"),collapse = "+")))
      print(multi_formulas)
      model <- coxph(multi_formulas, data = data)
      model.s <- summary(model)
      
      # Extract data
      tibble::tibble(
        Features = rownames(model.s$coefficients),
        n = model.s$n,
        coef = model.s$coefficients[,1],
        hr = model.s$conf.int[,1],
        hr_l = model.s$conf.int[,3],
        hr_h = model.s$conf.int[,4],
        coxp = model.s$coefficients[,5]) %>%
        dplyr::mutate(status = ifelse(hr > 1, "High_risk", "Low_risk")) %>%
        head(1)
    })) %>%
    tidyr::unnest()
}
# 3.2.3.drawing cox plot ------
fn_cox_plot_1 <- function(data,filename,title,facet, dir,w=4,h=4){
  data %>% 
    # dplyr::mutate(functionWithImmune=functionWithImmune) %>%
    # dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
    ggplot(aes(y = hr, x = Features, ymin=hr_l,ymax=hr_h)) +
    geom_pointrange(aes(color=cox_sig),size=0.5) +
    scale_color_manual(values=c("red","black")) +
    geom_hline(aes(yintercept = 1), linetype =2) +
    scale_size(name = "p-value") +
    scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6),
                       labels = c("1/64","1/32","1/16","1/8","1/4","1/2",1,2,3,4,5,6)) +
    # facet_grid(as.formula(facet),scales = "free", space = "free") +
    facet_wrap(as.formula(facet),scales = "free") +
    coord_flip() +
    # ggthemes::theme_gdocs() +
    my_theme +
    labs(y = "Hazard Ratio (High vs. low GSVA score)", x = "Features",title = title)+
    theme(
      legend.position = "none",
      axis.line.y = element_line(color="black"),
      axis.text = element_text(color = "black",size=8),
      axis.title = element_text(color = "black",size=10),
      text = element_text(color = "black"),
      title = element_text(size = 12)
    )  -> p;p
  ggsave(file.path(dir,paste(filename,"png",sep=".")),device = "png",width = w,height = h)
  ggsave(file.path(dir,paste(filename,"pdf",sep=".")),device = "pdf",width = w,height = h)
}

fn_cox_plot_2 <- function(data,filename,title,facet, dir,w=4,h=4){
  data %>% 
    # dplyr::mutate(functionWithImmune=functionWithImmune) %>%
    # dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
    ggplot(aes(y = hr, x = Features, ymin=hr_l,ymax=hr_h)) +
    geom_pointrange(aes(color=cox_sig),size=0.5) +
    scale_color_manual(values=c("red","black")) +
    geom_hline(aes(yintercept = 1), linetype =2) +
    scale_size(name = "p-value") +
    scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6),
                       labels = c("1/64","1/32","1/16","1/8","1/4","1/2",1,2,3,4,5,6)) +
    facet_grid(as.formula(facet),scales = "free", space = "free") +
    # facet_wrap(as.formula(facet),scales = "free") +
    coord_flip() +
    # ggthemes::theme_gdocs() +
    my_theme +
    theme(
      legend.position = "none",
      axis.line.y = element_line(color="black"),
      axis.text = element_text(color = "black",size=8),
      axis.title = element_text(color = "black",size=10),
      text = element_text(color = "black")
    ) +
    labs(y = "Hazard Ratio (High vs. low value)", x = "Features",title = title) -> p;p
  ggsave(file.path(dir,paste(filename,"png",sep=".")),device = "png",width = w,height = h)
  ggsave(file.path(dir,paste(filename,"pdf",sep=".")),device = "pdf",width = w,height = h)
}
# 3.3. univariate survival analysis ------
# 3.3.1.PFS ----
stage_class <- tibble::tibble(Stage = c("i/ii nos","is","Not_applicable",
                                        "stage i","stage ia","stage ib",
                                        "stage ii","stage iia","stage iib","stage iic",
                                        "stage iii","stage iiia","stage iiib","stage iiic",
                                        "stage iv","stage iva","stage ivb","stage ivc",
                                        "stage x",NA),
                              stage_group = c("Not_applicable","Not_applicable","Not_applicable",
                                        "stage i","stage i","stage i",
                                        "stage ii","stage ii","stage ii","stage ii",
                                        "stage iii","stage iii","stage iii","stage iii",
                                        "stage iv","stage iv","stage iv","stage iv",
                                        "Not_applicable","Not_applicable"))
clinical_all_data %>%
  dplyr::inner_join(stage_class,by="Stage") %>%
  dplyr::select(-Stage) %>%
  dplyr::rename("Stage"="stage_group") %>%
  dplyr::filter(Stage!="Not_applicable")-> clinical_all_data

GSVA_exp_ratio_combine %>%
  dplyr::mutate(clinical_gsva = purrr::map2(GSVA,cancer_types,.f=function(.x,.y){
    print(.y)
    .x %>%
      dplyr::left_join(CTL_data,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,15)) %>%
      dplyr::left_join(TIMER_immunity,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,12)) %>%
      dplyr::left_join(TCGA_mutation_burden %>%
                         dplyr::mutate(barcode=substr(Tumor_Sample_Barcode,1,12)) %>%
                         dplyr::select(barcode,mutation_durden),by="barcode") %>%
      dplyr::inner_join(clinical_all_data,by=c("barcode","cancer_types")) -> .tmp
    if(nrow(.tmp)!=0){
      .tmp %>%
        dplyr::rename("status"="PFS","time"="PFS.time") %>%
        dplyr::select(-OS,-Status,-cancer_types,-Stage) %>%
        tidyr::gather(-barcode,-status,-time,key="Features",value="value")  %>%
        dplyr::filter(!is.na(value)) %>%
        dplyr::group_by(Features) %>%
        dplyr::mutate(group = ifelse(value >= quantile(value,0.5),"2high","1low")) %>%
        dplyr::mutate(n=n()) %>%
        dplyr::ungroup() %>%
        tidyr::nest(-Features,.key = group) %>%
        dplyr::mutate(group=purrr::map(group,.f=function(.x){
          .x %>%
            dplyr::mutate(nn=.x$group %>% unique() %>% length())
        })) %>%
        tidyr::unnest() %>%
        dplyr::filter(n>10,nn>=2) %>%
        dplyr::select(-n,-nn) %>%
        # rbind(stage.surv.data) %>%
        tidyr::nest(-Features) 
    } else{
      tibble::tibble()
    }
  })) %>%
  dplyr::select(-GSVA) %>%
  tidyr::unnest() %>%
  dplyr::mutate(surv_res = purrr::map2(data,Features,fn_survival_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS
GSVA.score.univarite.surv.PFS %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","GSVA.score.univarite.surv.PFS.tsv"))

# 3.3.2.OS ----
GSVA_exp_ratio_combine %>%
  dplyr::mutate(clinical_gsva = purrr::map2(GSVA,cancer_types,.f=function(.x,.y){
    print(.y)
    .x %>%
      dplyr::left_join(CTL_data,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,15)) %>%
      dplyr::left_join(TIMER_immunity,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,12)) %>%
      dplyr::left_join(TCGA_mutation_burden %>%
                         dplyr::mutate(barcode=substr(Tumor_Sample_Barcode,1,12)) %>%
                         dplyr::select(barcode,mutation_durden),by="barcode") %>%
      dplyr::inner_join(clinical_all_data,by=c("barcode","cancer_types")) -> .tmp
    if(nrow(.tmp)!=0){
      .tmp %>%
        dplyr::rename("status"="Status","time"="OS") %>%
        dplyr::select(-PFS,-PFS.time,-cancer_types,-Stage) %>%
        tidyr::gather(-barcode,-status,-time,key="Features",value="value")  %>%
        dplyr::filter(!is.na(value)) %>%
        dplyr::group_by(Features) %>%
        dplyr::mutate(group = ifelse(value >= quantile(value,0.5),"2high","1low")) %>%
        dplyr::mutate(n=n()) %>%
        dplyr::ungroup() %>%
        tidyr::nest(-Features,.key = group) %>%
        dplyr::mutate(group=purrr::map(group,.f=function(.x){
          .x %>%
            dplyr::mutate(nn=.x$group %>% unique() %>% length())
        })) %>%
        tidyr::unnest() %>%
        dplyr::filter(n>10,nn>=2) %>%
        dplyr::select(-n,-nn) %>%
        # rbind(stage.surv.data) %>%
        tidyr::nest(-Features) 
    } else{
      tibble::tibble()
    }
  })) %>%
  dplyr::select(-GSVA) %>%
  tidyr::unnest() %>%
  dplyr::mutate(surv_res = purrr::map2(data,Features,fn_survival_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS
GSVA.score.univarite.surv.OS %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","2.GSVA.score.univarite.surv.OS.tsv"))

# 3.4. multi-variate survival analysis ------
# 3.4.1.PFS ----
GSVA.score.univarite.surv.PFS %>%
  dplyr::filter(coxp<=0.05 | kmp <= 0.05) %>%
  dplyr::mutate(Features = gsub(" ","_",Features)) %>%
  dplyr::select(cancer_types,Features) %>%
  tidyr::nest(-cancer_types,.key="sig_features") -> GSVA.score.univarite.surv.PFS.sig


GSVA_exp_ratio_combine %>%
  dplyr::mutate(clinical_gsva = purrr::map2(GSVA,cancer_types,.f=function(.x,.y){
    print(.y)
    colnames(.x) <- gsub(" ","_",colnames(.x))
    .x %>%
      dplyr::left_join(CTL_data,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,15)) %>%
      dplyr::left_join(TIMER_immunity,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,12)) %>%
      dplyr::left_join(TCGA_mutation_burden %>%
                         dplyr::mutate(barcode=substr(Tumor_Sample_Barcode,1,12)) %>%
                         dplyr::select(barcode,mutation_durden),by="barcode") %>%
      dplyr::inner_join(clinical_all_data,by="barcode")  %>%
      dplyr::rename("status"="PFS","time"="PFS.time") 
  })) %>%
  dplyr::select(-GSVA) %>%
  dplyr::inner_join(GSVA.score.univarite.surv.PFS.sig,by="cancer_types") %>%
  dplyr::mutate(surv_res.multi = purrr::map2(clinical_gsva,sig_features,fn_survival_test.multiCox)) %>%
  dplyr::select(-clinical_gsva,-sig_features) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.multi
GSVA.score.univarite.surv.PFS.multi %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","GSVA.score.multi-varite.surv.PFS.tsv"))

# 3.4.2.OS ----
GSVA.score.univarite.surv.OS %>%
  dplyr::filter(coxp<=0.05 | kmp <= 0.05) %>%
  dplyr::mutate(Features = gsub(" ","_",Features)) %>%
  dplyr::select(cancer_types,Features) %>%
  tidyr::nest(-cancer_types,.key="sig_features") -> GSVA.score.univarite.surv.OS.sig

GSVA_exp_ratio_combine %>%
  dplyr::mutate(clinical_gsva = purrr::map2(GSVA,cancer_types,.f=function(.x,.y){
    print(.y)
    colnames(.x) <- gsub(" ","_",colnames(.x))
    .x %>%
      dplyr::left_join(CTL_data,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,15)) %>%
      dplyr::left_join(TIMER_immunity,by=c("barcode")) %>%
      dplyr::mutate(barcode = substr(barcode,1,12)) %>%
      dplyr::left_join(TCGA_mutation_burden %>%
                         dplyr::mutate(barcode=substr(Tumor_Sample_Barcode,1,12)) %>%
                         dplyr::select(barcode,mutation_durden),by="barcode") %>%
      dplyr::inner_join(clinical_all_data,by="barcode")  %>%
      dplyr::rename("status"="Status","time"="OS") 
  })) %>%
  dplyr::select(-GSVA) %>%
  dplyr::inner_join(GSVA.score.univarite.surv.OS.sig,by="cancer_types") %>%
  dplyr::mutate(surv_res.multi = purrr::map2(clinical_gsva,sig_features,fn_survival_test.multiCox)) %>%
  dplyr::select(-clinical_gsva,-sig_features) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS.multi
GSVA.score.univarite.surv.OS.multi %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","2.GSVA.score.multi-varite.surv.OS.tsv"))

### 3.5.1.OS cox plot -------
# 3.5.1.1. cancer groups 
# GSVA.score.univarite.surv.OS.multi %>%
#   dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
#   dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
#   dplyr::filter(abs(hr)<10) %>%
#   dplyr::mutate(hr_l=ifelse(hr_l<(-10),-10,hr_l),hr_h=ifelse(hr_h>10,10,hr_h)) %>%
#   tidyr::nest(-cancer_types) %>%
#   dplyr::mutate(Group = "Group1") %>%
#   tidyr::unnest() %>%
#   tidyr::nest(-Group) %>%
#   dplyr::mutate(filename = paste("2.OS.multi-variable.cox",Group,sep="_")) %>%
#   dplyr::mutate(res = purrr::map2(data,filename,fn_cox_plot_2,facet = "cancer_types~.",title="Multi-variable cox model of GSVA score of features (OS)",dir=file.path(res_path,"3.survival_with_GSVA_score.new"),w=10,h=12))

# 3.5.1.2. all cancer togather
label.data.CTL %>%
  dplyr::select(Features,label,n) %>%
  dplyr::mutate(name = ifelse(n<5,paste(Features,"(",label,")",sep=""),Features)) %>%
  dplyr::mutate(Features=gsub(" ","_",Features)) -> Features_name.detail
GSVA.score.univarite.surv.OS.multi %>%
  dplyr::inner_join(Features_name.detail,by="Features") %>%
  dplyr::mutate(Features=label) %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.05,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  dplyr::filter(abs(hr)<6 & coxp<0.05) %>%
  dplyr::mutate(hr_l=ifelse(hr_l<(-6),-6,hr_l),hr_h=ifelse(hr_h>6,6,hr_h)) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(Group = "All") %>%
  tidyr::unnest() %>%
  tidyr::nest(-Group) %>%
  dplyr::mutate(filename = paste("2.OS.multi-variable.cox",Group,sep="_")) %>%
  dplyr::mutate(res = purrr::map2(data,filename,fn_cox_plot_1,facet = "cancer_types~.",title="Multi-variable cox(OS)\nGSVA score of ICP features",dir=file.path(res_path,"3.survival_with_GSVA_score.new"),w=6,h=20))

### 3.5.2.PFS cox plot -------
# 3.5.2.1. cancer groups 
# GSVA.score.univarite.surv.PFS.multi %>%
#   dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
#   dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
#   dplyr::filter(abs(hr)<10) %>%
#   dplyr::mutate(hr_l=ifelse(hr_l<(-10),-10,hr_l),hr_h=ifelse(hr_h>10,10,hr_h)) %>%
#   tidyr::nest(-cancer_types) %>%
#   dplyr::mutate(Group = c(rep("Group1",7),rep("Group2",7),rep("Group3",7),rep("Group4",8))) %>%
#   tidyr::unnest() %>%
#   tidyr::nest(-Group) %>%
#   dplyr::mutate(filename = paste("1.PFS.multi-variable.cox",Group,sep="_")) %>%
#   dplyr::mutate(res = purrr::map2(data,filename,fn_cox_plot_2,facet = "cancer_types~.",title="Multi-variable cox model of GSVA score of features (PFS)",dir=file.path(res_path,"3.survival_with_GSVA_score.new"),w=10,h=12))

# 3.5.2.2. all cancer togather
GSVA.score.univarite.surv.PFS.multi %>%
  dplyr::inner_join(Features_name.detail,by="Features") %>%
  dplyr::mutate(Features=label) %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.05,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  dplyr::filter(abs(hr)<6 & coxp<0.05) %>%
  dplyr::mutate(hr_l=ifelse(hr_l<(-6),-6,hr_l),hr_h=ifelse(hr_h>6,6,hr_h)) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(Group = "All") %>%
  tidyr::unnest() %>%
  tidyr::nest(-Group) %>%
  dplyr::mutate(filename = paste("1.PFS.multi-variable.cox",Group,sep="_")) %>%
  dplyr::mutate(res = purrr::map2(data,filename,fn_cox_plot_1,facet = "cancer_types~.",title="Multi-variable cox(PFS)\nGSVA score of ICP features",dir=file.path(res_path,"3.survival_with_GSVA_score.new"),w=6,h=20))

# 5. score from cox model -------------------------------------------------
# 5.1.function--------
fn_cox_score <- function(GSVA,survival){
  GSVA %>%
    tidyr::gather(-barcode,key="Features",value="GSVAscore") %>%
    dplyr::filter(Features %in% survival$Features) %>%
    dplyr::inner_join(survival,by="Features") %>%
    dplyr::mutate(GSVA_muliple_coef = GSVAscore*coef) %>%
    dplyr::group_by(barcode) %>%
    dplyr::mutate(cox_score = sum(GSVA_muliple_coef)) %>%
    dplyr::ungroup() %>%
    dplyr::select(barcode,cox_score) %>%
    unique()
}
# 5.2. multi-cox model coef ------
# 5.2.1.PFS --------
GSVA.score.univarite.surv.PFS.multi %>%
  tidyr::nest(-cancer_types,.key="PFS") %>%
  dplyr::inner_join(GSVA.score.onlytumor,by="cancer_types") %>%
  dplyr::mutate(cox_score = purrr::map2(GSVA,PFS,fn_cox_score)) %>%
  dplyr::select(cancer_types,cox_score) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.multi.coxscore

GSVA.score.univarite.surv.PFS.multi.coxscore %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","PFS_cox_score(GSVA*coef_of_multi_cox).tsv"))

GSVA.score.univarite.surv.PFS.multi.coxscore %>%
  dplyr::mutate(barcode = substr(barcode,1,15)) %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$TIL,method = "spearman"))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.multi.coxscore.correlation.TIL

GSVA.score.univarite.surv.PFS.multi.coxscore %>%
  # dplyr::mutate(barcode = substr(barcode,1,15)) %>%
  dplyr::inner_join(CTL_data,by=c("barcode")) %>%
  tidyr::nest(-cancer_types.x) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$CTL,method = "spearman"))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.multi.coxscore.correlation.CTL

GSVA.score.univarite.surv.PFS.multi.coxscore %>%
  dplyr::mutate(sample = substr(barcode,1,15)) %>%
  dplyr::inner_join(TCGA_mutation_burden,
                    by=c("sample")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$mutation_durden,method = "spearman"))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.multi.coxscore.correlation.mutation_durden
# survival analysis -------
# PFS survival all samples ----
GSVA.score.univarite.surv.PFS.multi.coxscore %>%
  dplyr::filter(substr(barcode,14,15)=="01") %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_all_data,by="barcode") %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(cox_score > quantile(cox_score,0.5),"high","low")) -> GSVA.score.univarite.surv.PFS.multi.coxscore.for_surv

GSVA.score.univarite.surv.PFS.multi.coxscore.for_surv %>%
  dplyr::mutate(x = "all") %>%
  tidyr::nest(-x) %>%
  dplyr::mutate(surv_res = purrr::map2(data,x,fn_survival_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest()

sur_res_path <- file.path(res_path,"3.survival_with_GSVA_score.new")
tibble::tibble(color = c("red","blue"),
               group = c("high","low"))  -> color_list
sur_name <- paste("3.PFS_between_PFS-multi_coxscore_high-low")

GSVA.score.univarite.surv.PFS.multi.coxscore.for_surv %>%
  fn_durvival_plot(title="Pan-cancer\nProgression-free survival",color=color_list,filename=sur_name,out_path=sur_res_path,legend.pos="none",h=3,w=4)

# PFS survival cancer specific ----
GSVA.score.univarite.surv.PFS.multi.coxscore.for_surv %>%
  dplyr::group_by(cancer_types.y) %>%
  dplyr::mutate(group = ifelse(cox_score > quantile(cox_score,0.5),"2high","1low")) %>%
  dplyr::ungroup() %>%
  tidyr::nest(-cancer_types.y) %>%
  dplyr::mutate(surv_res = purrr::map2(data,cancer_types.y,fn_survival_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.multi.coxscore.for_surv.res

GSVA.score.univarite.surv.PFS.multi.coxscore.for_surv.res %>%
  dplyr::arrange(kmp) %>%
  .$cancer_types.y -> cancer_types.y

GSVA.score.univarite.surv.PFS.multi.coxscore.for_surv.res %>%
  dplyr::mutate(value = -log10(kmp)) %>%
  ggplot(aes(x=cancer_types.y,y=value)) +
  geom_bar(aes(fill=status),stat="identity",alpha = 0.5,
           position=position_dodge()) +
  scale_x_discrete(limits=cancer_types.y) +
  scale_fill_manual(
    values = c("red","blue"),
    name = "HR(High/Low)"
  ) +
  geom_hline(yintercept = 1.3) +
  labs(x="Cancer types",
       y="-log10(KMP)") +
  my_theme +
  coord_flip()


sur_res_path <- file.path(res_path,"3.survival_with_GSVA_score.new")
tibble::tibble(color = c("red","blue"),
               group = c("high","low"))  -> color_list
sur_name <- paste("3.PFS_between_OS_coxscore_high-low")

GSVA.score.univarite.surv.PFS.multi.coxscore.for_surv %>%
  dplyr::group_by(cancer_types.x) %>%
  dplyr::mutate(group = ifelse(cox_score > quantile(cox_score,0.5),"high","low")) %>%
  dplyr::ungroup() %>%
  tidyr::nest(-cancer_types.x) %>%
  dplyr::mutate(title = paste(cancer_types.x,"\nOverall survival")) %>%
  dplyr::mutate(filename = paste("3.PFS_between_OS_coxscore_high-low",cancer_types.x,sep=".")) %>%
  dplyr::mutate(res = purrr::pmap(list(data=data,title=title,filename=filename),fn_durvival_plot,color=color_list,out_path=file.path(sur_res_path,"3.PFS_cancer_specific.PFS_coxscore"),legend.pos="none",h=3,w=4))

# 5.2.2.OS --------

GSVA.score.univarite.surv.OS.multi %>%
  tidyr::nest(-cancer_types,.key="OS") %>%
  dplyr::inner_join(GSVA.score.onlytumor,by="cancer_types") %>%
  dplyr::mutate(cox_score = purrr::map2(GSVA,OS,fn_cox_score)) %>%
  dplyr::select(cancer_types,cox_score) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS.multi.coxscore
GSVA.score.univarite.surv.OS.multi.coxscore %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","OS_cox_score(GSVA*coef_of_multi_cox).tsv"))

GSVA.score.univarite.surv.OS.multi.coxscore %>%
  dplyr::mutate(barcode = substr(barcode,1,15)) %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$TIL,method = "spearman"))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS.multi.coxscore.correlation.TIL

GSVA.score.univarite.surv.OS.multi.coxscore %>%
  # dplyr::mutate(barcode = substr(barcode,1,15)) %>%
  dplyr::inner_join(CTL_data,by=c("barcode")) %>%
  tidyr::nest(-cancer_types.x) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$CTL,method = "spearman"))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS.multi.coxscore.correlation.CTL

GSVA.score.univarite.surv.OS.multi.coxscore %>%
  dplyr::mutate(sample = substr(barcode,1,15)) %>%
  dplyr::inner_join(TCGA_mutation_burden,
                    by=c("sample")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$mutation_durden,method = "spearman"))
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS.multi.coxscore.correlation.mutation_durden

# 5.3 uni-cox coef to get cox_score --------
# 5.2.1.PFS --------
GSVA.score.univarite.surv.PFS %>%
  tidyr::nest(-cancer_types,.key="PFS") %>%
  dplyr::inner_join(GSVA.score.onlytumor,by="cancer_types") %>%
  dplyr::mutate(cox_score = purrr::map2(GSVA,PFS,fn_cox_score)) %>%
  dplyr::select(cancer_types,cox_score) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.PFS.coxscore

GSVA.score.univarite.surv.PFS.coxscore %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","3.PFS_cox_score(GSVA*coef_of_multi_cox)-all_features.tsv"))

GSVA.score.univarite.surv.PFS.coxscore %>%
  dplyr::mutate(barcode = substr(barcode,1,15)) %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$TIL,method = "spearman"))
  })) %>%
  dplyr::select(-data)  -> GSVA.score.univarite.surv.PFS.coxscore.correlation.TIL

GSVA.score.univarite.surv.PFS.coxscore %>%
  dplyr::inner_join(CTL_data,by=c("barcode")) %>%
  tidyr::nest(-cancer_types.x) %>%
  dplyr::mutate(CTL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$CTL,method = "spearman"))
  })) %>%
  dplyr::select(-data)-> GSVA.score.univarite.surv.PFS.coxscore.correlation.CTL

GSVA.score.univarite.surv.PFS.coxscore %>%
  dplyr::mutate(sample = substr(barcode,1,15)) %>%
  dplyr::inner_join(TCGA_mutation_burden,
                    by=c("sample")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(MB_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$mutation_durden,method = "spearman"))
  })) %>%
  dplyr::select(-data) -> GSVA.score.univarite.surv.PFS.coxscore.correlation.MB

GSVA.score.univarite.surv.PFS.coxscore.correlation.TIL %>%
  dplyr::inner_join(GSVA.score.univarite.surv.PFS.coxscore.correlation.CTL,by="cancer_types") %>%
  dplyr::inner_join(GSVA.score.univarite.surv.PFS.coxscore.correlation.MB,by="cancer_types") -> GSVA.score.univarite.surv.PFS.coxscore.correlation.combine

GSVA.score.univarite.surv.PFS.coxscore.correlation.combine %>%
  readr::write_rds(file.path(res_path,"3.survival_with_GSVA_score.new","3.PFS_cox_score(GSVA*coef_of_multi_cox)-all_features-correlation.combine(TIL.CTL.MB).rds.gz"),compress = "gz")

# 5.2.2.OS --------
GSVA.score.univarite.surv.OS %>%
  tidyr::nest(-cancer_types,.key="PFS") %>%
  dplyr::inner_join(GSVA.score.onlytumor,by="cancer_types") %>%
  dplyr::mutate(cox_score = purrr::map2(GSVA,PFS,fn_cox_score)) %>%
  dplyr::select(cancer_types,cox_score) %>%
  tidyr::unnest() -> GSVA.score.univarite.surv.OS.coxscore

GSVA.score.univarite.surv.OS.coxscore %>%
  readr::write_tsv(file.path(res_path,"3.survival_with_GSVA_score.new","3.OS_cox_score(GSVA*coef_of_multi_cox)-all_features.tsv"))
# 5.2.2.1.correlation with TIl,CTL and MB --------
GSVA.score.univarite.surv.OS.coxscore %>%
  dplyr::mutate(barcode = substr(barcode,1,15)) %>%
  dplyr::inner_join(TIMER_immunity,by="barcode") %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(TIL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$TIL,method = "spearman"))
  })) %>%
  dplyr::select(-data)  -> GSVA.score.univarite.surv.OS.coxscore.correlation.TIL

GSVA.score.univarite.surv.OS.coxscore %>%
  dplyr::inner_join(CTL_data,by=c("barcode","cancer_types")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(CTL_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$CTL,method = "spearman"))
  })) %>%
  dplyr::select(-data)-> GSVA.score.univarite.surv.OS.coxscore.correlation.CTL

GSVA.score.univarite.surv.OS.coxscore %>%
  dplyr::mutate(sample = substr(barcode,1,15)) %>%
  dplyr::inner_join(TCGA_mutation_burden,
                    by=c("sample")) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(MB_cor = purrr::map(data,.f=function(.x){
    broom::tidy(cor.test(.x$cox_score,.x$mutation_durden,method = "spearman"))
  })) %>%
  dplyr::select(-data) -> GSVA.score.univarite.surv.OS.coxscore.correlation.MB

GSVA.score.univarite.surv.OS.coxscore.correlation.TIL %>%
  dplyr::inner_join(GSVA.score.univarite.surv.OS.coxscore.correlation.CTL,by="cancer_types") %>%
  dplyr::inner_join(GSVA.score.univarite.surv.OS.coxscore.correlation.MB,by="cancer_types") -> GSVA.score.univarite.surv.OS.coxscore.correlation.combine

GSVA.score.univarite.surv.OS.coxscore.correlation.combine %>%
  readr::write_rds(file.path(res_path,"3.survival_with_GSVA_score.new","3.OS_cox_score(GSVA*coef_of_multi_cox)-all_features-correlation.combine(TIL.CTL.MB).rds.gz"),compress = "gz")

# 5.2.2.2.survival between high and low cox score --------
# survival plot fucntion ---
fn_durvival_plot <- function(data,title="Pan-cancers progress-free survival",color,filename="Pan-can.PFS.survival",out_path,legend.pos,h=3,w=4){
  
  fit <- survfit(Surv(time, status) ~ group, data = data, na.action = na.exclude)
  tmp <- 0
  for (i in 1:length(fit$strata)) {
    tmp <- fit$strata[i] + tmp
    .group <- strsplit(names(fit$strata[i]),split = "=")[[1]][2]
    tibble::tibble(mid_time=fit$time,final_surv=fit$surv) %>%
      head(tmp) %>%
      tail(1) %>%
      dplyr::mutate(group = .group)-> res.tmp
    
    if(i==1){
      tail_position <- res.tmp
    } else{
      tail_position <- rbind(tail_position,res.tmp)
    }
  }
  diff <- survdiff(Surv(time, status) ~ group, data = data, na.action = na.exclude)
  kmp <- 1 - pchisq(diff$chisq, df = length(levels(as.factor(data$group))) - 1) %>% signif(2)
  color %>%
    dplyr::inner_join(data,by="group") %>%
    dplyr::inner_join(tail_position,by="group") %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(label = paste(group,", n=",n,sep="")) %>%
    dplyr::select(group,label,color,n,mid_time,final_surv) %>%
    dplyr::arrange(group) %>%
    unique() -> text_data
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        data = data,
                        # surv.median.line = "hv",
                        title = paste(title,", p=",kmp,sep=""), # change it when doing diff data
                        xlab = "Time (days)",
                        ylab = 'Probability of survival',
                        # legend.title = "Methyla group:",
                        legend= c(legend.pos),
                        # ggtheme = theme_survminer(),
                        ggtheme = theme(
                          panel.border = element_blank(), 
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), 
                          axis.line = element_line(colour = "black",size = 0.5),
                          panel.background = element_rect(fill = "white"),
                          legend.key = element_blank(),
                          legend.background = element_blank(),
                          legend.text = element_text(size = 8),
                          axis.text = element_text(size = 12,color = "black"),
                          legend.title = element_blank(),
                          axis.title = element_text(size = 12,color = "black"),
                          text = element_text(color = "black")
                        )
  )[[1]] +
    # geom_text(aes(x=time,y=surv,label=label),data=text_data) +
    geom_text_repel(aes(x=mid_time,y=final_surv,label=label,color=group),data=text_data) +
    scale_color_manual(
      values = c(text_data$color,text_data$color),
      labels = c(text_data$group,text_data$group)
    ) -> p;p
  ggsave(filename = paste(filename,kmp,"pdf",sep = "."), plot = p, path = out_path,device = "pdf",height = h,width = w)
  ggsave(filename = paste(filename,kmp,"png",sep = "."), plot = p, path = out_path,device = "png",height = h,width = w)
}

# PFS survival all samples ----
GSVA.score.univarite.surv.OS.coxscore %>%
  dplyr::filter(substr(barcode,14,15)=="01") %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(clinical_all_data,by="barcode") %>%
  dplyr::rename("time"="PFS.time","status"="PFS") %>%
  dplyr::mutate(group = ifelse(cox_score > quantile(cox_score,0.5),"high","low")) -> GSVA.score.univarite.surv.OS.coxscore.for_surv

GSVA.score.univarite.surv.OS.coxscore.for_surv %>%
  dplyr::mutate(x = "all") %>%
  tidyr::nest(-x) %>%
  dplyr::mutate(surv_res = purrr::map2(data,x,fn_survival_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest()

sur_res_path <- file.path(res_path,"3.survival_with_GSVA_score.new")
tibble::tibble(color = c("red","blue"),
               group = c("high","low"))  -> color_list
sur_name <- paste("3.PFS_between_OS_coxscore_high-low")

GSVA.score.univarite.surv.OS.coxscore.for_surv %>%
  fn_durvival_plot(title="Pan-cancer\nOverall survival",color=color_list,filename=sur_name,out_path=sur_res_path,legend.pos="none",h=3,w=4)

# PFS survival cancer specific ----
GSVA.score.univarite.surv.OS.coxscore.for_surv %>%
  dplyr::group_by(cancer_types.x) %>%
  dplyr::mutate(group = ifelse(cox_score > quantile(cox_score,0.5),"high","low")) %>%
  dplyr::ungroup() %>%
  tidyr::nest(-cancer_types.x) %>%
  dplyr::mutate(surv_res = purrr::map2(data,cancer_types.x,fn_survival_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest()

sur_res_path <- file.path(res_path,"3.survival_with_GSVA_score.new")
tibble::tibble(color = c("red","blue"),
               group = c("high","low"))  -> color_list
sur_name <- paste("3.PFS_between_OS_coxscore_high-low")

GSVA.score.univarite.surv.OS.coxscore.for_surv %>%
  dplyr::group_by(cancer_types.x) %>%
  dplyr::mutate(group = ifelse(cox_score > quantile(cox_score,0.5),"high","low")) %>%
  dplyr::ungroup() %>%
  tidyr::nest(-cancer_types.x) %>%
  dplyr::mutate(title = paste(cancer_types.x,"\nOverall survival")) %>%
  dplyr::mutate(filename = paste("3.PFS_between_OS_coxscore_high-low",cancer_types.x,sep=".")) %>%
  dplyr::mutate(res = purrr::pmap(list(data=data,title=title,filename=filename),fn_durvival_plot,color=color_list,out_path=file.path(sur_res_path,"3.PFS_cancer_specific.OS_coxscore"),legend.pos="none",h=3,w=4))

# 6. filter Features for each cancers -------------------------------------
GSVA.TIL.res %>%
  tidyr::unnest() %>%
  dplyr::filter(abs(estimate)>=0.3 & p.value<=0.05) %>%
  dplyr::filter(abs(`log2FC(High/Low)`)>=0.585 & p.value1<=0.05) %>%
  dplyr::group_by(Features) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>2) %>%
  dplyr::mutate(Analysis = "TIL_cor_DE") %>%
  dplyr::select(cancer_types,Features,Analysis) %>%
  unique() %>%
  dplyr::ungroup()-> GSVA.TIL.res.sig
GSVA.MB.res %>%
  tidyr::unnest() %>%
  dplyr::filter(abs(estimate)>=0.3 & p.value<=0.05) %>%
  dplyr::filter(abs(`log2FC(High/Low)`)>=0.585 & p.value1<=0.05) %>%
  dplyr::mutate(Analysis = "MB_cor_DE") %>%
  dplyr::select(cancer_types,Features,Analysis)-> GSVA.MB.res.sig
GSVA.CTL.res %>%
  tidyr::unnest() %>%
  dplyr::filter(abs(estimate)>=0.3 & p.value<=0.05) %>%
  dplyr::filter(abs(`log2FC(High/Low)`)>=1 & p.value1<=0.05)%>%
  dplyr::mutate(Analysis = "CTL_cor_DE") %>%
  dplyr::select(cancer_types,Features,Analysis) -> GSVA.CTL.res.sig
GSVA.score.univarite.surv.OS.multi %>%
  dplyr::filter(coxp<=0.1) %>%
  dplyr::mutate(Analysis = "Multi-cox,OS") %>%
  dplyr::select(cancer_types,Features,Analysis) -> GSVA.score.univarite.surv.OS.multi.sig
GSVA.score.univarite.surv.PFS.multi %>%
  dplyr::filter(coxp<=0.1) %>%
  dplyr::mutate(Analysis = "Multi-cox,PFS") %>%
  dplyr::select(cancer_types,Features,Analysis) -> GSVA.score.univarite.surv.PFS.multi.sig

GSVA.TIL.res.sig %>%
  rbind(GSVA.MB.res.sig) %>%
  rbind(GSVA.CTL.res.sig) %>%
  rbind(GSVA.score.univarite.surv.OS.multi.sig) %>%
  rbind(GSVA.score.univarite.surv.PFS.multi.sig) %>%
  dplyr::mutate(Features=gsub(" ","_",Features)) %>%
  dplyr::group_by(Features) %>%
  dplyr::mutate(Sig_n = n()) %>%
  dplyr::select(Features,Sig_n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Sig_n) %>%
  dplyr::rename("Counts_of_significant_analysis_in_all_cancers"="Sig_n") %>%
  dplyr::filter(!Features %in% c(colnames(TIMER_immunity),"CTL","Age")) %>%
  dplyr::inner_join(Features_name.detail,by="Features") %>%
  dplyr::mutate(color = ifelse(Counts_of_significant_analysis_in_all_cancers>=3,"red","black")) -> GSVA_score_res.counts

GSVA_score_res.counts %>%
  readr::write_tsv(file.path(res_path,"overall_feature_filter.tsv"))

GSVA_score_res.counts %>%
  # dplyr::filter(cancer_types %in% c("SKCM","STAD","LUSC","LUAD","KIRC","KICH","KIRP")) %>%
  ggplot(aes(x=label,y=Counts_of_significant_analysis_in_all_cancers)) +
  geom_bar(aes(color=color),stat="identity", fill="orange",alpha = 0.5,
           position=position_dodge()) +
  scale_x_discrete(limits=GSVA_score_res.counts$label) +
  scale_color_manual(values=c("black","red")) +
  coord_flip() +
  my_theme +
  theme(
    axis.text.y = element_text(color = GSVA_score_res.counts$color),
    legend.position = "none"
  ) +
  labs(y="Counts of significant analysis among cancers",x="Features") 
ggsave(file.path(res_path,"overall_feature_filter.pdf"),device = "pdf",height = 20,width =6)
ggsave(file.path(res_path,"overall_feature_filter.png"),device = "png",height = 20,width = 6)

# save.image --------------------------------------------------------------

# save.image(file.path(res_path,"TCGA.GSVA_score.ICP_features-all_cancer_togather.rda"))
