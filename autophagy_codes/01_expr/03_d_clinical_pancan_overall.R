###########################################################################
# survival situation of each cancers ------------------------------

library(magrittr)
library(ggplot2)

# data path config ---------------------------------------------------
basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
tcga_path <- file.path(basic_path,"/data/TCGA")
# gene_list_path <- file.path(immune_path,"checkpoint/20171021_checkpoint")
# expr_path <- file.path(basic_path,"immune_checkpoint/result_20171025/expr_rds")
out_path <- file.path(basic_path,"immune_checkpoint/result_20171025/c_3_survival")


# load data ---------------------------------------------------------------

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
  dplyr::rename("cancer_types"= "type")-> clinical_all_data


# draw survival plot ------------------------------------------------------
library(survival)
library(ggrepel)
library(ggplot2)
## function
fn_durvival_plot <- function(fit,clinical_all_data,title="Pan-cancers progress-free survival",text_data,filename="Pan-can.PFS.survival"){
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        data = clinical_all_data,
                        # surv.median.line = "hv",
                        title = paste(title), # change it when doing diff data
                        xlab = "Time (days)",
                        ylab = 'Probability of survival',
                        # legend.title = "Methyla group:",
                        legend= "bottom",
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
    geom_text_repel(aes(x=time,y=surv,label=label),data=text_data) +
    scale_color_manual(
      values = text_data$color,
      labels = text_data$cancer_types
    ) -> p;p
  ggsave(filename = paste(filename,"pdf",sep = "."), plot = p, path = out_path,device = "pdf",height = 8,width = 10)
  ggsave(filename = paste(filename,"png",sep = "."), plot = p, path = out_path,device = "png",height = 8,width = 10)
}
## PFS
fit <- survfit(Surv(PFS.time, PFS) ~ cancer_types, data = clinical_all_data, na.action = na.exclude)
fit.pfs <- fit
tmp <- 0
for (i in 1:length(fit$strata)) {
  tmp <- fit$strata[i] + tmp
  cancers <- strsplit(names(fit$strata[i]),split = "=")[[1]][2]
  tibble::tibble(time=fit$time,surv=fit$surv) %>%
    head(tmp) %>%
    tail(1) %>%
    dplyr::mutate(cancer_types = cancers)-> res.tmp
  
  if(i==1){
    tail_position <- res.tmp
  } else{
    tail_position <- rbind(tail_position,res.tmp)
  }
}

readr::read_tsv(file.path(tcga_path,"02_pcc.tsv")) %>%
  dplyr::inner_join(clinical_all_data,by="cancer_types") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(tail_position,by="cancer_types") %>%
  dplyr::mutate(label = paste(cancer_types,", n=",n,sep="")) %>%
  dplyr::select(cancer_types,label,color,n,time,surv) %>%
  dplyr::arrange(cancer_types) %>%
  unique() -> color_paired

fn_durvival_plot(fit,clinical_all_data,title = "Pan-cancers progress-free survival",text_data = color_paired,filename = "Pan-can.PFS.survival")

## OS

fit <- survfit(Surv(OS, Status) ~ cancer_types, data = clinical_all_data, na.action = na.exclude)
fit.os <- fit
tmp <- 0
for (i in 1:length(fit$strata)) {
  tmp <- fit$strata[i] + tmp
  cancers <- strsplit(names(fit$strata[i]),split = "=")[[1]][2]
  tibble::tibble(time=fit$time,surv=fit$surv) %>%
    head(tmp) %>%
    tail(1) %>%
    dplyr::mutate(cancer_types = cancers)-> res.tmp
  
  if(i==1){
    tail_position <- res.tmp
  } else{
    tail_position <- rbind(tail_position,res.tmp)
  }
}

readr::read_tsv(file.path(tcga_path,"02_pcc.tsv")) %>%
  dplyr::inner_join(clinical_all_data,by="cancer_types") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(tail_position,by="cancer_types") %>%
  dplyr::mutate(label = paste(cancer_types,", n=",n,sep="")) %>%
  dplyr::select(cancer_types,label,color,n,time,surv) %>%
  dplyr::arrange(cancer_types) %>%
  unique() -> color_paired

fn_durvival_plot(fit,clinical_all_data,title = "Pan-cancers overall survival",text_data = color_paired,filename = "Pan-can.OS.survival")

# middle survival
summary(fit.os)
# OS
fit.os.s <- summary(fit.os)
os.median <- fit.os.s$table[,"median"]

for (i in 1:length(names(fit.os$strata))) {
  name <- names(fit.os$strata)[i]
  surv_rate <- fit.os.s$surv[max(grep("TRUE",fit.os.s$strata==name))]
  cancers <- strsplit(name,"=")[[1]][2]
  res.tmp <- tibble::tibble(cancer_types = cancers, survival_rate = surv_rate)
  if(i==1){
    os.survrate <- res.tmp
  }else{
    os.survrate <- rbind(os.survrate,res.tmp)
  }
}
os.median %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::mutate(cancer_types = names(os.median)) %>%
  dplyr::mutate(cancer_types = purrr::map(cancer_types,.f=function(.x){
    strsplit(.x,"=")[[1]][2]
  })) %>%
  tidyr::unnest() %>%
  dplyr::rename("Median_survival_time" = ".") %>%
  dplyr::inner_join(os.survrate,by="cancer_types") %>%
  dplyr::arrange(survival_rate,Median_survival_time ) %>%
  dplyr::mutate(n=1:n())  %>%
  dplyr::mutate(Median_survival_time=ifelse(is.na(Median_survival_time),Inf,Median_survival_time)) -> OS.pancan.rank

# PFS
fit.pfs.s <- summary(fit.pfs)
pfs.median <- fit.pfs.s$table[,"median"]
for (i in 1:length(names(fit.pfs$strata))) {
  name <- names(fit.pfs$strata)[i]
  surv_rate <- fit.pfs.s$surv[max(grep("TRUE",fit.pfs.s$strata==name))]
  cancers <- strsplit(name,"=")[[1]][2]
  res.tmp <- tibble::tibble(cancer_types = cancers, survival_rate = surv_rate)
  if(i==1){
    pfs.survrate <- res.tmp
  }else{
    pfs.survrate <- rbind(os.survrate,res.tmp)
  }
}

pfs.median %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::mutate(cancer_types = names(pfs.median)) %>%
  dplyr::mutate(cancer_types = purrr::map(cancer_types,.f=function(.x){
    strsplit(.x,"=")[[1]][2]
  })) %>%
  tidyr::unnest() %>%
  dplyr::rename("Median_survival_time" = ".") %>%
  dplyr::inner_join(pfs.survrate,by="cancer_types") %>%
  dplyr::arrange(survival_rate,Median_survival_time ) %>%
  dplyr::mutate(n=1:n()) %>%
  dplyr::mutate(Median_survival_time=ifelse(is.na(Median_survival_time),Inf,Median_survival_time)) -> PFS.pancan.rank
  
  # dplyr::filter(!is.na(Median_PFS_survival_time)) 

# function
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

fn_histrom_line <- function(data,title,filename,w=6,h=5){
  data %>%
    dplyr::mutate(box="Median survival time",line="Ultimate survival rate") %>%
    ggplot() +
    geom_bar(aes(x=n,y=Median_survival_time,fill = box),stat="identity") +
    geom_line(aes(x=n, y=survival_rate*6000,color = line),stat="identity") +
    scale_y_continuous(sec.axis = sec_axis(~./6000, name = "Ultimate survival rate")) +
    scale_x_continuous(breaks = c(1:nrow(data)),
                       labels = data$cancer_types) +
    labs(title = title,y="Median survival time (days)") +
    scale_fill_manual(values="tan1") +
    scale_color_manual(values="blue") +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
      legend.title = element_blank(),
      legend.position = "top",
      axis.title.x = element_blank()
    )
  ggsave(file.path(out_path,paste(filename,"pdf",sep=".")),device = "pdf", width = w, height = h)
  ggsave(file.path(out_path,paste(filename,"png",sep=".")),device = "png", width = w, height = h)
}
fn_histrom_line(data=PFS.pancan.rank,title = "Progression-free survival",filename = "Pan-can.PFS.rank")
fn_histrom_line(data=OS.pancan.rank,title = "Overall survival",filename = "Pan-can.OS.rank")
