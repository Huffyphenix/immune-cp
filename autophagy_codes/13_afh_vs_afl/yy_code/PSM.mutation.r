# check data
folder <- "Mutation"
if (!file.exists(folder)) { dir.create(folder) }
scripts.dir <- paste("/extraspace/yye1/analysis/Hypoxia/PSM/",folder,sep="")
setwd(scripts.dir)
analysis="myclusters" #Oxygen_Content
sum.mutAll <- data.frame()
####must exist stum  clinical, stratification, mutation files.
mut.files.names <- list.files(path="/extraspace/TCGA/Mutation/",pattern=".txt")
mut.files.Abs <- gsub("mutation_|.txt","",mut.files.names)
clinical.files.names <- list.files(path="/extraspace/TCGA/TCGA_clinical/",pattern=".txt")
clinical.files.Abs <- gsub("_clinical_clean.txt","",clinical.files.names)
stratification.files.names <- list.files(path="/extraspace/yye1/analysis/Hypoxia/Stratification.tumor/",pattern="Hypoxia.stratification.txt")
SampleCount <- read.delim("/extraspace/yye1/analysis/Hypoxia/Stratification.tumor/Hypoxia and normoxia sample size large than 30.txt",header=T)
stratification.files.Abs <- intersect(SampleCount$type,intersect(clinical.files.Abs,mut.files.Abs))
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
    # mutation
    mutation <- read.delim(paste("/extraspace/TCGA/Mutation/mutation_",cancer,".txt",sep=""),header=T)
    mutation.pri <- data.frame(t(mutation[,2:ncol(mutation)]))
    colnames(mutation.pri) <- as.vector(mutation[,1])
    # only focuse on the highly mutated genes
    # if the top 100 mutated gene mutation frequency less than 0.05, keep top 100 mutated genes, or keep all genes which matutation frequency large than 0.05
    top <- 100
    mut.cutoff <- 0.05
    mut.rate <- colSums(mutation.pri)/nrow(mutation.pri)
    mut.rate.sort <- sort(mut.rate, decreasing=T, index.return=T)
    if(length(which(mut.rate > mut.cutoff)) > 100){
      keep.index <- mut.rate.sort$ix[match( names(which(mut.rate.sort$x > mut.cutoff)),names(mut.rate.sort$x))]
    }else{
      top.index <-mut.rate.sort$ix[1:top]
      top.cutoff <- mut.rate.sort$x[top] # keep ties
      keep.index <- which(mut.rate>= top.cutoff)
    }
    #folder <- "sample_list"
    #if (!file.exists(folder)) { dir.create(folder) }
    # For LUSC: gender, ATT/IPW could not get balanced weight. OVERLAP and MW get balanced results. Double check!
    folder <- paste(cancer,"_",analysis,sep="")
    if (!file.exists(folder)) { dir.create(folder) }
    mut.result <- weight.test(data.dum, form, mutation.pri[,keep.index], is.continuous=FALSE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, "mut", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=FALSE)
    sum.mut <- summarize.fdr(mutation.pri[,keep.index], mut.result, print=TRUE)
    summarize.fdr(mutation.pri[,keep.index], mut.result, print=TRUE, cutoff=0.05)
    write.summary(sum.mut, cancer, analysis,"mut")
    write.result(mut.result, cancer, analysis,"mut")
    if(length(which(mut.result$fdr < 0.05)) > 0){
      sum.mut <- data.frame(sum.mut)
      sum.mut$class <- rep(cancer,times=nrow(sum.mut))
      if(nrow(sum.mutAll) == 0){
        sum.mutAll <- sum.mut
      }else{
        sum.mutAll <- rbind(sum.mutAll,sum.mut)
      }
      perm.result <- foreach (seed =1:100, .combine='rbind') %dopar%
      {
        ## mutation, only for KIRC        
        perm.mut.result <- weight.test(data.dum, form, mutation.pri[,keep.index], is.continuous=FALSE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=TRUE, cancer, "mut", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=TRUE, seed=seed)
        perm.sum.mut <- summarize.fdr(mutation.pri[,keep.index], perm.mut.result)
        
        write(c(seed, perm.sum.mut$n.sig), file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_count_",seed,".txt", sep=""), ncolumns=7)
        save(seed,perm.mut.result, file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
        
      }
      
      
      cutoff <- 0.05
      seedV <- 1:100
      perm.cal(cancer, analysis, "mut", mutation.pri[,keep.index], cutoff=cutoff, seedV=seedV)
      
      if(FALSE)
      {
        mut.ttest <- myttest(data.dum, mutation.pri, cancer,"mut")
        sum.mut <- summarize.fdr(mutation.pri, mut.ttest)
        save(mRNAseq.ttest, file=paste(cancer,"_mut_ttest.RData", sep=""))
      }
    }
  }
 
}

write.table(sum.mutAll,file="Mutation.genes.across.cancer.types.txt",quote = F,sep="\t",row.names = F)


####
MutationCountAll <- c("Type","Mutation.Num")
for( cancer in stratification.files.Abs){
  mutation <- read.delim(paste("/extraspace/TCGA/Mutation/mutation_",cancer,".txt",sep=""),header=T)
  mutation.pri <- data.frame(t(mutation[,2:ncol(mutation)]))
  colnames(mutation.pri) <- as.vector(mutation[,1])
  mut.rate <- colSums(mutation.pri)/nrow(mutation.pri)
  mutation.count <- c(cancer,length(mut.rate[mut.rate > 0.05]))
  MutationCountAll <- rbind(MutationCountAll,mutation.count)
}
write.csv(MutationCountAll,file="/extraspace/yye1/analysis/Hypoxia/PSM/Mutation/MutationCountAll.csv",quote = F,row.names = F)




