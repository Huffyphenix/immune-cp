# check data
rm(list=ls())
setwd("/extraspace/yye1/analysis/Hypoxia/PSM")
cancerNames <- c("BLCA","BRCA","CESC","ESCA","GBM","HNSC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","SKCM","PRAD","SARC","STAD","THCA","TGCT","THYM","UCEC")
folder <- "cnv.lesion"
if (!file.exists(folder)) { dir.create(folder) }
scripts.dir <- paste("/extraspace/yye1/analysis/Hypoxia/PSM/",folder,sep="")
setwd(scripts.dir)
analysis="myclusters" #Oxygen_Content
sum.cnvAll <- data.frame()
####must exist stum  clinical, stratification, cnv files.
cnv.dirs.names <- list.dirs(path="/extraspace/yye1/share_data/TCGA_CNV_GISTIC/")
cnv.dirs.names <- cnv.dirs.names[2:length(cnv.dirs.names)]
cnv.files.Abs <- gsub("/extraspace/yye1/share_data/TCGA_CNV_GISTIC//GDAC_","",cnv.dirs.names)
clinical.files.names <- list.files(path="/extraspace/TCGA/TCGA_clinical/",pattern=".txt")
clinical.files.Abs <- gsub("_clinical_clean.txt","",clinical.files.names)
stratification.files.names <- list.files(path="/extraspace/yye1/analysis/Hypoxia/Stratification.tumor/",pattern="Hypoxia.stratification.txt")
stratification.files.Abs <- gsub(".Hypoxia.stratification.txt","",stratification.files.names)
SampleCount <- read.delim("/extraspace/yye1/analysis/Hypoxia/Stratification.tumor/Hypoxia and normoxia sample size large than 30.txt",header=T)
stratification.files.Abs <- intersect(SampleCount$type,intersect(clinical.files.Abs,cnv.files.Abs))
stratification.files.Abs <- intersect(stratification.files.Abs,cancerNames)

for(cancer in stratification.files.Abs){
  clinical.file <- paste("/extraspace/TCGA/TCGA_clinical/",cancer,"_clinical_clean.txt", sep="")
  clinical.raw <- read.table(clinical.file, header=TRUE, check.names=FALSE, sep="\t",quote = "", comment.char="")
  print(paste("Total clinical samples for", cancer, ":", nrow(clinical.raw)))
  ###add hypoxic and normoxic stratification data
  stratification <- read.delim(paste("/extraspace/yye1/analysis/Hypoxia/Stratification.tumor/",cancer,".Hypoxia.stratification.txt",sep=""),header=T)
  stratification$barcode <- gsub("\\.","\\-",substr(stratification$SampleID,1,12))
  stratification <- stratification[which(stratification$myclusters %in% c("hypoxic","normoxic")),c("myclusters","barcode")]
  clinical.raw <- merge(stratification,clinical.raw,by="barcode")
  if(nrow(clinical.raw) >=50){
    end.col <- which(colnames(clinical.raw)=="os_days")-1
    #str(clinical.raw[,1:end.col])
    #summary(clinical.raw[,1:end.col])
    # apply(clinical.raw[1:end.col], 2, function(x) length(which(is.na(x))))
    
    # process data
    data <- clinical.raw[,1:end.col]
    ## rm race_ASIAN and race_NA
    consider_factors <- c("barcode","myclusters","age_at_initial_pathologic_diagnosis","gender","race","pathologic_stage")
    consider_factors <- intersect(consider_factors,colnames(data))
    if(length(table(data$pathologic_stage))==0){
      consider_factors <- consider_factors[which(consider_factors != "pathologic_stage")]
      data <- data[,consider_factors]
    }else{
      data <- data[,consider_factors]
    }
    
    if(length(which(is.na(data$race))) >0){
      data <- data[-which(data$race=="ASIAN" | is.na(data$race)),]
    }
    if(length(which(is.na(data$pathologic_stage))) > 0){
      data <- data[-which(is.na(data$pathologic_stage)),]
    }
    
    data$race <- as.vector(data$race)
    #data$race[which(is.na(data$race))] <- "UNKNOWN"
    data$race <- factor(data$race)
    data <- unique(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    
    analysis <- "myclusters"
    if(analysis=="myclusters"){
      # convert hypoxic and normoxic to numeric 1,0 to suppress the warning message in lm
      data$myclusters <- ifelse(data$myclusters=="hypoxic",1,0)
      colnames(data)[which(colnames(data)=="myclusters")] <- "Z"
    }
    
    # convert to dummy
    library(dummies)
    dummy.feature <- setdiff(colnames(data),c("Z","age_at_initial_pathologic_diagnosis"))#,"pathologic_stage"))
    data.dum <- dummy.data.frame(data, names=dummy.feature)
    dummy.list <- attr(data.dum,"dummies")
    rm.col <- c()
    for (i in 1:length(dummy.list))
    {
      rm.col <- c(rm.col, dummy.list[[i]][length(dummy.list[[i]])])
    }
    data.dum <- data.dum[,-rm.col]
    data.dum$X0 <- rep(1, nrow(data.dum))
    #form <- as.formula("Z~.") # should exclude X0
    exclude.col <- match(c("Z","X0"), colnames(data.dum))
    colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
    form <- as.formula(paste("Z~",paste(colnames(data.dum)[-exclude.col],collapse="+"),sep=""))
    # perform calculation
    source("~/code/cal.R")
    library(doMC)
    library(foreach)
    registerDoMC(15)
    # cnv.exp  /extraspace/yye1/share_data/TCGA_CNV_GISTIC/
    cnv <- read.delim(paste("/extraspace/yye1/share_data/TCGA_CNV_GISTIC/GDAC_",cancer,"/all_lesions.conf_99.txt",sep=""),header=T)
    cnv <- cnv[grep("CN values",cnv$Unique.Name),]
    cnv$Unique.Name <- gsub("Peak|- CN values|Peak","",cnv$Unique.Name)
    cnv$Unique.Name <- gsub(" ","",cnv$Unique.Name)
    cnv <- cnv[,c("Unique.Name",grep("TCGA",colnames(cnv),value=T))]
    
    if(length(duplicated(substr(colnames(cnv),1,12))[duplicated(substr(colnames(cnv),1,12))==TRUE]) >0){
      cnv <- cnv[,colnames(cnv)[-which(colnames(cnv) %in% colnames(cnv)[duplicated(substr(colnames(cnv),1,12))])]]
    }
    cnv[,2:ncol(cnv)] <- apply(cnv[,2:ncol(cnv)],2,function(x){2^(1+x)})
    cnv.pri <- data.frame(t(cnv[,2:ncol(cnv)]))
    colnames(cnv.pri) <- as.vector(cnv[,1])
    rownames(cnv.pri) <- substr(rownames(cnv.pri),1,12)
    cnv.pri <- rm.zero.col(cnv.pri)
    
    folder <- paste(cancer,"_",analysis,sep="")
    if (!file.exists(folder)) { dir.create(folder) }
    cnv.result <- weight.test(data.dum, form, cnv.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, data.type= "cnv", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=FALSE)
    sum.cnv <- summarize.fdr(cnv.pri, cnv.result, print=TRUE)
    summarize.fdr(cnv.pri, cnv.result, print=TRUE, cutoff=0.05)
    write.summary(sum.cnv, cancer, analysis,"cnv")
    write.result(cnv.result, cancer, analysis,"cnv")
    save(cnv.result, file=paste(cancer,"_", analysis,"_result.RData",sep=""))
    if(length(which(cnv.result$fdr < 0.05)) > 0){
      sum.cnv <- data.frame(sum.cnv)
      sum.cnv$class <- rep(cancer,times=nrow(sum.cnv))
      if(nrow(sum.cnvAll) == 0){
        sum.cnvAll <- sum.cnv
      }else{
        sum.cnvAll <- rbind(sum.cnvAll,sum.cnv)
      }
      perm.result <- foreach (seed =1:100, .combine='rbind') %dopar%
      {
        ## cnv, only for KIRC        
        perm.cnv.result <- weight.test(data.dum, form, cnv.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, "cnv", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=TRUE, seed=seed)
        perm.sum.cnv <- summarize.fdr(cnv.pri, perm.cnv.result)
        
        write(c(seed, perm.sum.cnv$n.sig), file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_count_",seed,".txt", sep=""), ncolumns=7)
        save(seed,perm.cnv.result, file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
        
      }
      
      
      cutoff <- 0.05
      seedV <- 1:100
      perm.cal(cancer, analysis, "cnv", cnv.pri, cutoff=cutoff, seedV=seedV)
      
      if(FALSE)
      {
        cnv.ttest <- myttest(data.dum, cnv.pri, cancer,"cnv")
        sum.cnv <- summarize.fdr(cnv.pri, cnv.ttest)
        save(cnv.ttest, file=paste(cancer,"_cnv_ttest.RData", sep=""))
      }
    }
  }
  
}

write.table(sum.cnvAll,file="cnv.genes.across.cancer.types.txt",quote = F,sep="\t",row.names = F)

cnvsamples <- list.files(path="/extraspace/yye1/analysis/Hypoxia/PSM/cnv.lesion/",pattern="_cnv_samples.txt")
cnvsamplenames <- gsub("_cnv_samples.txt","",cnvsamples)
for(cancer in intersect(cnvsamplenames,cancerNames)){
  cnvsample <- read.delim(paste("/extraspace/yye1/analysis/Hypoxia/PSM/cnv.lesion/",cancer,"_cnv_samples.txt",sep=""),header=T)
  cnvsample$class <- rep(cancer,nrow(cnvsample))
  if(cancer == intersect(cnvsamplenames,cancerNames)){
    cnvsampleAll <- cnvsample
  }else{
    cnvsampleAll <- rbind(cnvsampleAll,cnvsample)
  }
}

cnvCount <- t(sapply(split(cnvsampleAll[,"V2"],cnvsampleAll$class),function(x){table(factor(x,levels=c("hypoxic","normoxic")))}))
cnvCount <- data.frame(cnvCount)
cnvCount$class <- rownames(cnvCount)
write.table(cnvCount,file = "/extraspace/yye1/analysis/Hypoxia/PSM/cnv.lesion/cnvSampleCount.txt")

for(cancer in cancerNames){
  cnv <- read.delim(paste("/extraspace/yye1/share_data/TCGA_CNV_GISTIC/GDAC_",cancer,"/all_lesions.conf_99.txt",sep=""),header=T)
  cnv <- cnv[grep("CN values",cnv$Unique.Name),]
  cnv$Unique.Name <- gsub("Peak|- CN values|Peak","",cnv$Unique.Name)
  cnv$Unique.Name <- gsub(" ","",cnv$Unique.Name)
  cnv <- cnv[,c("Unique.Name",grep("TCGA",colnames(cnv),value=T))]
  SCNA <- data.frame(cancer,Amp=length(grep("Amp",cnv$Unique.Name)),Del=length(grep("Del",cnv$Unique.Name)))
  if(cancer==cancerNames[1]){
    SCNAs <- SCNA
  }else{
    SCNAs <- rbind(SCNAs,SCNA)
  }
 
}
SCNAs <- SCNAs[order(SCNAs$Amp+SCNAs$Del),]
sign.cnv.PSM <- read.delim("/extraspace/yye1/analysis/Hypoxia/PSM/cnv.lesion/cnv.genes.across.cancer.types.txt",header=T)
sign.cnv.PSM <- sign.cnv.PSM[which(sign.cnv.PSM$class %in% cancerNames),]
sign.cnvCount <- data.frame(t(sapply(split(sign.cnv.PSM[,c("feature.sig","coef.sig")],sign.cnv.PSM$class),function(x){c(nrow(x[which(substr(x$feature.sig,1,3)=="Amp" & x$coef.sig > 0),]),nrow(x[which(substr(x$feature.sig,1,3)=="Amp" & x$coef.sig < 0),]),nrow(x[which(substr(x$feature.sig,1,3)=="Del" & x$coef.sig > 0),]),nrow(x[which(substr(x$feature.sig,1,3)=="Del" & x$coef.sig < 0),]))})))
colnames(sign.cnvCount) <- c("Hypo.Amp","Norm.Amp","Hypo.Del","Hypo.Amp")
sign.cnvCount$cancer <- rownames(sign.cnvCount)
sign.cnvCount <- merge(SCNAs,sign.cnvCount,by="cancer")
sign.cnvCount <- sign.cnvCount[!(sign.cnvCount$cancer %in% c("ESCA","PRAD")),]
write.csv(sign.cnvCount,file="/extraspace/yye1/analysis/Hypoxia/PSM/cnv.lesion/results/cnv hypoxia-associated lesion and background.csv",quote = F,row.names = F)
