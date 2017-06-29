###############
#read in data
###############
read.table("mut_rate.filter_cancer.txt",sep="\t",header = T,row.names = 1)->mut_rate.filter_cancer.txt
read.table("number_of_samples_in_cancers",sep="\t",row.names = 1)->number_of_samples_in_cancers
mut_num.filter.txt[rownames(mut_rate.filter_cancer.txt),]->mut_num.filter.txt
mut_num.filter.txt[,colnames(mut_rate.filter_cancer.txt)]->mut_num.filter.txt
cbind(colnames(mut_rate.filter_cancer.txt),number_of_samples_in_cancers[colnames(mut_rate.filter_cancer.txt),"V2"])->sample.n

###################################
#sample names change as cancer(n=..)
###################################
sample.n1<-c()
for(i in 1:nrow(sample.n)){
  sample.n1<-c(sample.n1,paste(sample.n[i,1],"(n=",sample.n[i,2],")",sep=""))
}
colnames(mut_rate.filter_cancer.txt)=sample.n1

#############
#draw heatmap
#############
library(pheatmap)
mut_rate.filter_cancer.txt<-mut_rate.filter_cancer.txt*100
pdf(file = "mut_rate.pdf")
pheatmap(mut_rate.filter_cancer.txt,cluster_row = FALSE,cluster_cols = FALSE,
         color =colorRampPalette(c("white","red"))(20),border_color="grey60",
         legend_breaks=seq(0,15,5),legend_labels = c("0","5","10","14"),
         display_numbers = mut_num.filter.txt,
         main = c("Mutation rate of check-point genes in 32 cancers."))
dev.off()
