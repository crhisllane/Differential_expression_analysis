##########################################
## Build count matrix with HTSeq counts ##
##########################################
# Define the names of counts files for all samples

direc = "/home/crhisllane/Doutorado/Suzane/htseq_count";
setwd(direc)
count.files<-c("C1_htseq_result.out", "C2_htseq_result.out", "C3_htseq_result.out", "L4_htseq_result.out", "L5_htseq_result.out", "L6_htseq_result.out")

# Build the count matrix
count_matrix<-lapply(1:length(count.files), function(x){
  count.vec<-read.delim(file=paste(count.files[x]), header=FALSE, stringsAsFactor=FALSE)
  return(count.vec)
})
count_matrix<-do.call(cbind, count_matrix)

##take a look on data
head(count_matrix)
tail(count_matrix)

#V1      V2                     V1      V2                     V1
#19866          CPIJ040898-RA       0          CPIJ040898-RA       2          CPIJ040898-RA
#19867           __no_feature 1697998           __no_feature 1806531           __no_feature
#19868            __ambiguous    9007            __ambiguous   10112            __ambiguous
#19869        __too_low_aQual       0        __too_low_aQual       0        __too_low_aQual
#19870          __not_aligned 3116073          __not_aligned  971522          __not_aligned
#19871 __alignment_not_unique  909549 __alignment_not_unique  718848 __alignment_not_unique
#V2                     V1      V2                     V1      V2
#19866       0          CPIJ040898-RA       0          CPIJ040898-RA       0
#19867 2820138           __no_feature 1600830           __no_feature 1319851
#19868   16285            __ambiguous    6090            __ambiguous    6750
#19869       0        __too_low_aQual       0        __too_low_aQual       0
#19870  464710          __not_aligned  316727          __not_aligned  547068
#19871 1087766 __alignment_not_unique 1003293 __alignment_not_unique  649457
#V1      V2
#19866          CPIJ040898-RA       0
#19867           __no_feature 1596838
#19868            __ambiguous    6926
#19869        __too_low_aQual       0
#19870          __not_aligned  389022
#19871 __alignment_not_unique 1214950



# exclude the last 5 rows (they are aligment informaion)
count_matrix<-count_matrix[5:nrow(count_matrix),]
head(count_matrix)
#the first column of each count.file has the gene name
#set gene names as row names
row.names(count_matrix)<-count_matrix[,1]
# keep counts only
count_matrix<-count_matrix[,c(2,4,6,8,10,12)]
#set names of count files as column names
colnames(count_matrix)<-count.files

#explore the count matrix
head(count_matrix)
##how meny trasncript do we have?
dim(count_matrix)
#[1] 19867     6

# build 2 data frames. The first one will contain sample information and the second gene annotation
sample_data<-data.frame(sample=names(count_matrix), condition=c(rep("Field", times=3), rep("Lab", times=3)))


# We must ensure that all genes have 10 counts in at least one condition
summary(rowSums(count_matrix[,1:3]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       1      11     878      66 6324667 
summary(rowSums(count_matrix[,4:6]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0       0       2     536      23 4517519

# build a index of genes with good counts
#in this case good counts are defined as those satisfying count >= 10
index_good_counts<-sapply(1:nrow(count_matrix), function(x){
  idx<-all(count_matrix[x,1:3] > rep(9, 3)) | all(count_matrix[x,4:6] > rep(9,3)) 
  return(idx)
})
table(index_good_counts)
#index_good_counts
#FALSE  TRUE 
#13953  5914

# keep those genes that have good counts
count_matrix<-count_matrix[index_good_counts,]
dim(count_matrix)
# [1] 5914    6

# save filtered count matrix, gene data and sample data 
dados<-list(count_matrix, sample_data)
save(dados, file="filteredCountMatrix.RData", compress="xz")

######################
## Data exploration ##
######################

head(sample_data$condition)
#[1] Field Field Field Lab   Lab   Lab  
#Levels: Field Lab

library(DESeq2)

# We use count_matrix as countData, sample_data as colData and we choose condition column of sample_data as principal factor to the Binomial Negative Model
y_DESeq<-DESeqDataSetFromMatrix(countData=count_matrix, colData=sample_data, design= ~condition)
#set size factor equal 1 (no library size effect) to avoind rlogTransformation does for you
sizeFactors(y_DESeq)<-rep(1,6)

# To the first exploration, we make a log2 conversion over count data
#this data does not have any factor correction
matrix_rlog_DESeq<-rlogTransformation(y_DESeq, blind=TRUE)

head(count_matrix_rlog_DESeq<-assays(matrix_rlog_DESeq)[[1]])
names(count_matrix_rlog_DESeq)<-names(count_matrix)
head(count_matrix_rlog_DESeq)
#takes a look on data distributions

###################################
##PCA over log2 transformed data ##
###################################
# The first exploration must be done over our count_matrix in order to determine if our samples are clustered togheter using count values
# Principal Component Analysis: explore samples separability  
pr_comp_y<-prcomp(t((count_matrix_rlog_DESeq)))
resumen <- summary(pr_comp_y)
labX <- signif((resumen$importance[2,1])*100, 3)
labY <- signif((resumen$importance[2,2])*100, 3)
png(file="PCA_DESeq_counts.png", width=1024, height=840)
par(mfrow=c(1,2))
plot(x=pr_comp_y$x[,1], y=pr_comp_y$x[,2], col=c(rep("red",3), rep("blue",3)), xlab=paste("PC1 (",labX,"%)", sep=""), ylab=paste("PC2 (",labY,"%)", sep=""))
abline(h=0)
abline(v=0)
title(main="PCA DESeq counts")
biplot(pr_comp_y,cex=c(1,0.001),xlabs=c("Field1","Field2","Field3","Lab1","Lab2","Lab3"),ylabs=rep("",nrow(count_matrix_rlog_DESeq)),var.axes=FALSE)
title(main="Biplot DESeq counts")
dev.off()
#######################################
##Boxplot over log2 transformed data ##
#######################################
# The second exploration is over all counts for each sample in order to determine if samples are comparable
# Boxplot of the raw and normalized read counts
png(file="comparisson_of_boxplot_distributions.png", width=2000, height=1200)
par(mfrow=c(1,2))
boxplot(count_matrix_rlog_DESeq, col=colors()[c(137, 134, 59,30,132,125)], main="DESeq", names=c("Field1","Field2","Field3","Lab1","Lab2","Lab3"))
FDP_DESeq<-sapply(1:ncol(count_matrix_rlog_DESeq), function(x){
  fdp<-density(count_matrix_rlog_DESeq[,x])
  return(fdp)})
FDP_DESeq[,1]$call<-"log2counts Densities"
plot(FDP_DESeq[,1], col=colors()[137], xlab="log2DESeqCounts", ylab="Density", xlim=c(0,16), ylim=c(0,0.250),type="l", lwd=4)
lines(FDP_DESeq[,2],col=colors()[134], lwd=4)
lines(FDP_DESeq[,3],col=colors()[59], lwd=4)
lines(FDP_DESeq[,4],col=colors()[30], lwd=4)
lines(FDP_DESeq[,5],col=colors()[132], lwd=4)
lines(FDP_DESeq[,6],col=colors()[125], lwd=4)
legend(legend=c("Field1","Field2","Field3","Lab1","Lab2","Lab3"),x="topright", fill = colors()[c(137, 134, 59,30,132,125)])
dev.off()

################################
## Library size normalization ##
################################
# There are some diferences between density functions. We should remember that existing differences in library size
# In order to estimate those diferences, we use estimateSizeFactors function
y_DESeq<-estimateSizeFactors(y_DESeq)
sizeFactors(y_DESeq)
#C1_htseq_result.out C2_htseq_result.out C3_htseq_result.out L4_htseq_result.out 
#1.2580496           1.3347290           2.2918593           0.6380649 
#L5_htseq_result.out L6_htseq_result.out 
#0.7597429           0.6195442 
norm_rlog_DESeqData<-rlogTransformation(y_DESeq, blind=TRUE)
count_matrix_norm_rlog_DESeqData<-assays(norm_rlog_DESeqData)[[1]]
colnames(count_matrix_norm_rlog_DESeqData)<-sample_data$sample
head(count_matrix_norm_rlog_DESeqData)

png(file="comparisson_of_boxplot_distributions_afterNorm.png", width=2000, height=1200)
par(mfrow=c(2,2))
boxplot(count_matrix_rlog_DESeq, col=colors()[c(137,  59,101,30,125,128)], main="Rlog Counts", names=c("Field1","Field2","Field3","Lab1","Lab2","Lab3"))

FDP<-sapply(1:ncol(count_matrix_rlog_DESeq), function(x){
  fdp<-density(count_matrix_rlog_DESeq[,x])
  return(fdp)})

FDP_DESeq<-sapply(1:ncol(count_matrix_norm_rlog_DESeqData), function(x){
  fdp<-density(count_matrix_norm_rlog_DESeqData[,x])
  return(fdp)})
FDP[,1]$call<-"log2counts Densities"
plot(FDP[,1], col=colors()[137], xlab="log2DESeqcounts", ylab="Density", xlim=c(0,16), ylim=c(0,0.30),type="l", lwd=4)
lines(FDP[,2],col=colors()[59], lwd=4)
lines(FDP[,3],col=colors()[101], lwd=4)
lines(FDP[,4],col=colors()[30], lwd=4)
lines(FDP[,5],col=colors()[125], lwd=4)
lines(FDP[,6],col=colors()[128], lwd=4)
legend(legend=c("Field1","Field2","Field3","Lab1","Lab2","Lab3"),x="topright", fill = colors()[c(137,  59,101,30,125,128)])

boxplot(count_matrix_norm_rlog_DESeqData, col=colors()[c(137,  59,101,30,125,128)], main="Normalized Rlog Counts", names=c("Field1","Field2","Field3","Lab1","Lab2","Lab3"))
FDP_DESeq[,1]$call<-"log2counts Densities"
plot(FDP_DESeq[,1], col=colors()[137], xlab="log2DESeqcounts", ylab="Density", xlim=c(0,16), ylim=c(0,0.30),type="l", lwd=4)
lines(FDP_DESeq[,2],col=colors()[59], lwd=4)
lines(FDP_DESeq[,3],col=colors()[101], lwd=4)
lines(FDP_DESeq[,4],col=colors()[30], lwd=4)
lines(FDP_DESeq[,5],col=colors()[125], lwd=4)
lines(FDP_DESeq[,6],col=colors()[128], lwd=4)
legend(legend=c("Field1","Field2","Field3","Lab1","Lab2","Lab3"),x="topright", fill = colors()[c(137, 59,101,30,125,128)])
dev.off()

#########################################
## Binomial Negative Model and DE Test ##
#########################################
# Before DE testing is necessary estimate dispersion of the model
y_DESeq<-estimateDispersions(y_DESeq, fitType="local")
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates

#Test de DE
et_DESeq<-nbinomWaldTest(y_DESeq)
results.DESeq <- results(et_DESeq)
head(results.DESeq)
head(results(et_DESeq)[order(results(et_DESeq)$padj),])

# Outlier's analysis. For outlier detection, we should determine which of all genes has its pvalue as NA
table(is.na(results(et_DESeq)$pvalue))
#FALSE 
#5914

cond1<-"Field"
cond2<-"Lab"
dds=et_DESeq
decide<-function(dss, p.adj=TRUE, value=0.05, FC_teste=1,cond1,cond2){
  if(class(dds) != "DESeqDataSet") stop("dds must be a DESeqDataSet")
  var<-"padj"
  if(!p.adj) var<-"pvalue"
  dds<-as.data.frame(results(dds,contrast = c("condition",cond1,cond2)))
  de_dds<-rep(0, nrow(dds))
  de_dds<-sapply(1:nrow(dds), function(x){
    
    if(!is.na(dds[x, var]) & dds[x, var] <= value & dds[x,"log2FoldChange"] >= FC_teste ) {
      aux<-1} else{
        if(!is.na(dds[x, var]) & dds[x, var] <= value & dds[x,"log2FoldChange"] <= -(FC_teste)) {
          aux<-(-1)
        }else {
          if(is.na(dds[x, var])) {
            aux<-NA}else aux<-0
        }}
    return(aux)
  })
  return(de_dds)
}
# We consider DE genes those that have a adjusted pavalue lower than 0.01
# table(de_DESeq <-decide(et_DESeq_clean, p.adj=TRUE, value=0.01 ))
table(de_DESeq <-decide(dss, p.adj=TRUE, value=0.05,FC_teste = 1,"Field","Lab" ))
#-1    0    1 
#294 5234  386

DESeqDEResults<-results(et_DESeq,contrast = c("condition","Field","Lab"))
write.table(DESeqDEResults,"Field-vs-Lab_DEanalysis.tab")


############ Explore DE results ############ 


## 1)  check if samples are clustered together using counts of DE genes
library(gplots)
select<-which(de_DESeq ==1 | de_DESeq== -1)
png(file="HeatMap_log2counts_UP_DOWN_genes.png",width=1024, height=860)
hm<-heatmap.2(scale(t(as.matrix(log2(counts(et_DESeq[select,], normalized=TRUE)+1))), center=TRUE, scale=FALSE),dendrogram="col", scale="none", col="redgreen",trace="none", labRow=c("Field1","Field2","Field3","Lab1","Lab2","Lab3"),labCol="",margins = c(10, 7))
dev.off()

## 2) check expression bias. 
png(file="MA_DE_genes_DESeq.png", width=1200, height=900)
par(mfrow=c(1,1))
plot(rowMeans(log2(counts(et_DESeq, normalized=TRUE))),DESeqDEResults$log2FoldChange, xlab="Mean log2normalizedCounts", ylab="log2FoldChange", col=8,cex=0.5)
points(rowMeans(log2(counts(et_DESeq, normalized=TRUE)))[de_DESeq==1], DESeqDEResults$log2FoldChange[de_DESeq==1], col=3,cex=0.8)
points(rowMeans(log2(counts(et_DESeq, normalized=TRUE)))[de_DESeq==-1], DESeqDEResults$log2FoldChange[de_DESeq==-1], col=2,cex=0.8)
title(main="MA plot logCPM DESeq_counts")
legend(legend=c("UP Genes","DOWN Genes"),x="topright", fill = c(3,2))
abline(h=0, col=1, lwd=2)
dev.off()

FDP_RLE<-density(DESeqDEResults$log2FoldChange[!is.na(DESeqDEResults$log2FoldChange)])
FDP_RLE$call<-"log2FC Density"
png(file="log2FC_DESeq_distribution.png",width=1024, height=860)
par(mfrow=c(1,2))
boxplot(DESeqDEResults$log2FoldChange, col=c(3))
plot(FDP_RLE, col=3, lwd=2)
dev.off()
