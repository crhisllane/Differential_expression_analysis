
##########################################
## Build count matrix with HTSeq counts ##
##########################################
# Define the names of counts files for all samples

direc = "/home/crhisllane/Doutorado/Suzane";
setwd(direc)
count.files<-c("v2_C1_ReadsPerGene_reverse.out.tab", "v2_C2_ReadsPerGene_reverse.out.tab", "v2_C3_ReadsPerGene_reverse.out.tab", "v2_L4_ReadsPerGene_reverse.out.tab", "v2_L5_ReadsPerGene_reverse.out.tab", "v2_L6_ReadsPerGene_reverse.out.tab")


# Build the count matrix
count_matrix<-lapply(1:length(count.files), function(x){
  count.vec<-read.delim(file=paste(count.files[x]), header=FALSE, stringsAsFactor=FALSE)
  return(count.vec)
})
count_matrix<-do.call(cbind, count_matrix)

##take a look on data
head(count_matrix)
# V1   V2            V1   V2            V1    V2            V1    V2            V1   V2            V1    V2
# 1 CPIJ012138-RA    0 CPIJ012138-RA    0 CPIJ012138-RA     0 CPIJ012138-RA     0 CPIJ012138-RA    0 CPIJ012138-RA     0
# 2 CPIJ012139-RA    1 CPIJ012139-RA    3 CPIJ012139-RA     4 CPIJ012139-RA     0 CPIJ012139-RA    0 CPIJ012139-RA     0
# 3 CPIJ012140-RA    2 CPIJ012140-RA    4 CPIJ012140-RA     5 CPIJ012140-RA     1 CPIJ012140-RA    3 CPIJ012140-RA     3
# 4 CPIJ012141-RA   28 CPIJ012141-RA   25 CPIJ012141-RA    28 CPIJ012141-RA     9 CPIJ012141-RA    9 CPIJ012141-RA     7
# 5 CPIJ012142-RA    3 CPIJ012142-RA    3 CPIJ012142-RA     9 CPIJ012142-RA     3 CPIJ012142-RA    2 CPIJ012142-RA     0
# 6 CPIJ012143-RA 8134 CPIJ012143-RA 7315 CPIJ012143-RA 11488 CPIJ012143-RA 13692 CPIJ012143-RA 7358 CPIJ012143-RA 11234

tail(count_matrix)
# V1 V2            V1 V2            V1 V2            V1 V2            V1 V2            V1 V2
# 19861 CPIJ019929-RA 29 CPIJ019929-RA 28 CPIJ019929-RA 51 CPIJ019929-RA 21 CPIJ019929-RA 12 CPIJ019929-RA  5
# 19862 CPIJ020103-RA  0 CPIJ020103-RA  0 CPIJ020103-RA  0 CPIJ020103-RA  0 CPIJ020103-RA  0 CPIJ020103-RA  0
# 19863 CPIJ020273-RA  0 CPIJ020273-RA  0 CPIJ020273-RA  0 CPIJ020273-RA  0 CPIJ020273-RA  0 CPIJ020273-RA  0
# 19864 CPIJ018847-RA  0 CPIJ018847-RA  0 CPIJ018847-RA  0 CPIJ018847-RA  0 CPIJ018847-RA  0 CPIJ018847-RA  0
# 19865 CPIJ019515-RA  0 CPIJ019515-RA  0 CPIJ019515-RA  0 CPIJ019515-RA  0 CPIJ019515-RA  0 CPIJ019515-RA  0
# 19866 CPIJ020086-RA  1 CPIJ020086-RA  0 CPIJ020086-RA  1 CPIJ020086-RA  0 CPIJ020086-RA  0 CPIJ020086-RA  0


## the first column of each count.file has the gene name
## set gene names as row names
row.names(count_matrix)<-count_matrix[,1]
# keep counts only
count_matrix<-count_matrix[,c(2,4,6,8,10,12)]
#set names of count files as column names
colnames(count_matrix)<-count.files

#explore the count matrix
head(count_matrix)
# v2_C1_ReadsPerGene_reverse.out.tab v2_C2_ReadsPerGene_reverse.out.tab v2_C3_ReadsPerGene_reverse.out.tab
# CPIJ012138-RA                                  0                                  0                                  0
# CPIJ012139-RA                                  1                                  3                                  4
# CPIJ012140-RA                                  2                                  4                                  5
# CPIJ012141-RA                                 28                                 25                                 28
# CPIJ012142-RA                                  3                                  3                                  9
# CPIJ012143-RA                               8134                               7315                              11488
# v2_L4_ReadsPerGene_reverse.out.tab v2_L5_ReadsPerGene_reverse.out.tab v2_L6_ReadsPerGene_reverse.out.tab
# CPIJ012138-RA                                  0                                  0                                  0
# CPIJ012139-RA                                  0                                  0                                  0
# CPIJ012140-RA                                  1                                  3                                  3
# CPIJ012141-RA                                  9                                  9                                  7
# CPIJ012142-RA                                  3                                  2                                  0
# CPIJ012143-RA                              13692                               7358                              11234


##how many trasncript do we have?
dim(count_matrix)
# [1] 19866     6

# build 2 data frames. The first one will contain sample information and the second gene annotation
sample_data<-data.frame(sample=names(count_matrix), condition=c(rep("Field", times=3), rep("Lab", times=3)))

# We must ensure that all genes have 10 counts in at least one condition
summary(rowSums(count_matrix[,1:3]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     0.0     7.0   107.6    42.0 66750.0 

summary(rowSums(count_matrix[,4:6]))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00     0.00     2.00    56.15    16.00 40770.00


# build a index of genes with good counts
#in this case good counts are defined as those satisfying count >= 10
index_good_counts<-sapply(1:nrow(count_matrix), function(x){
  idx<-all(count_matrix[x,1:3] > rep(9, 3)) | all(count_matrix[x,4:6] > rep(9,3)) 
  return(idx)
})
table(index_good_counts)
# index_good_counts
# FALSE  TRUE 
# 15086  4780 


# keep those genes that have good counts
count_matrix<-count_matrix[index_good_counts,]
dim(count_matrix)
# [1] 4780    6


# save filtered count matrix, gene data and sample data 
dados<-list(count_matrix, sample_data)
save(dados, file="filteredCountMatrix.RData", compress="xz")

######################
## Data exploration ##
######################

head(sample_data$condition)
# [1] Field Field Field Lab   Lab   Lab  
# Levels: Field Lab


library(DESeq2)

# We use count_matrix as countData, sample_data as colData and we choose condition column of sample_data as principal factor to the Binomial Negative Model
y_DESeq<-DESeqDataSetFromMatrix(countData=count_matrix, colData=sample_data, design= ~condition)
#set size factor equal 1 (no library size effect) to avoind rlogTransformation does for you
sizeFactors(y_DESeq)<-rep(1,6)

# To the first exploration, we make a log2 conversion over count data
#this data does not have any factor correction
matrix_rlog_DESeq<-rlogTransformation(y_DESeq, blind=TRUE)

head(count_matrix_rlog_DESeq<-assays(matrix_rlog_DESeq)[[1]])
# v2_C1_ReadsPerGene_reverse.out.tab v2_C2_ReadsPerGene_reverse.out.tab v2_C3_ReadsPerGene_reverse.out.tab
# CPIJ012141-RA                           4.402173                           4.314787                           4.402173
# CPIJ012143-RA                          13.085714                          12.996301                          13.391514
# CPIJ007879-RA                           4.212466                           3.893633                           4.212466
# CPIJ007889-RA                           3.962913                           4.319828                           4.729265
# CPIJ007892-RA                           5.336209                           6.020603                           6.210744
# CPIJ007894-RA                           4.039105                           4.431734                           4.535559
# v2_L4_ReadsPerGene_reverse.out.tab v2_L5_ReadsPerGene_reverse.out.tab v2_L6_ReadsPerGene_reverse.out.tab
# CPIJ012141-RA                           3.659157                           3.659157                           3.537552
# CPIJ012143-RA                          13.550374                          13.001161                          13.371304
# CPIJ007879-RA                           3.552090                           3.361964                           3.608402
# CPIJ007889-RA                           4.009291                           4.319828                           4.176268
# CPIJ007892-RA                           5.357050                           5.377477                           5.126607
# CPIJ007894-RA                           3.658761                           3.998572                           3.320453

names(count_matrix_rlog_DESeq)<-names(count_matrix)
head(count_matrix_rlog_DESeq)
# v2_C1_ReadsPerGene_reverse.out.tab v2_C2_ReadsPerGene_reverse.out.tab v2_C3_ReadsPerGene_reverse.out.tab
# CPIJ012141-RA                           4.402173                           4.314787                           4.402173
# CPIJ012143-RA                          13.085714                          12.996301                          13.391514
# CPIJ007879-RA                           4.212466                           3.893633                           4.212466
# CPIJ007889-RA                           3.962913                           4.319828                           4.729265
# CPIJ007892-RA                           5.336209                           6.020603                           6.210744
# CPIJ007894-RA                           4.039105                           4.431734                           4.535559
# v2_L4_ReadsPerGene_reverse.out.tab v2_L5_ReadsPerGene_reverse.out.tab v2_L6_ReadsPerGene_reverse.out.tab
# CPIJ012141-RA                           3.659157                           3.659157                           3.537552
# CPIJ012143-RA                          13.550374                          13.001161                          13.371304
# CPIJ007879-RA                           3.552090                           3.361964                           3.608402
# CPIJ007889-RA                           4.009291                           4.319828                           4.176268
# CPIJ007892-RA                           5.357050                           5.377477                           5.126607
# CPIJ007894-RA                           3.658761                           3.998572                           3.320453




# takes a look on data distributions

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
plot(FDP_DESeq[,1], col=colors()[137], xlab="log2DESeqCounts", ylab="Density", xlim=c(0,16), ylim=c(0,0.400),type="l", lwd=4)
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
# v2_C1_ReadsPerGene_reverse.out.tab v2_C2_ReadsPerGene_reverse.out.tab v2_C3_ReadsPerGene_reverse.out.tab 
# 1.2114707                          1.3045608                          2.1610338 
# v2_L4_ReadsPerGene_reverse.out.tab v2_L5_ReadsPerGene_reverse.out.tab v2_L6_ReadsPerGene_reverse.out.tab 
# 0.6527135                          0.7709879                          0.6415473


norm_rlog_DESeqData<-rlogTransformation(y_DESeq, blind=TRUE)
count_matrix_norm_rlog_DESeqData<-assays(norm_rlog_DESeqData)[[1]]
colnames(count_matrix_norm_rlog_DESeqData)<-sample_data$sample
head(count_matrix_norm_rlog_DESeqData)
# v2_C1_ReadsPerGene_reverse.out.tab v2_C2_ReadsPerGene_reverse.out.tab v2_C3_ReadsPerGene_reverse.out.tab
# CPIJ012141-RA                           4.176752                           4.057928                           3.829990
# CPIJ012143-RA                          12.922771                          12.774771                          12.731812
# CPIJ007879-RA                           3.994693                           3.695047                           3.660788
# CPIJ007889-RA                           3.894998                           4.133818                           4.149897
# CPIJ007892-RA                           5.203498                           5.755195                           5.527103
# CPIJ007894-RA                           3.879145                           4.154831                           3.925499
# v2_L4_ReadsPerGene_reverse.out.tab v2_L5_ReadsPerGene_reverse.out.tab v2_L6_ReadsPerGene_reverse.out.tab
# CPIJ012141-RA                           3.873941                           3.795465                           3.773431
# CPIJ012143-RA                          13.959583                          13.225297                          13.784956
# CPIJ007879-RA                           3.754553                           3.517977                           3.814699
# CPIJ007889-RA                           4.255958                           4.452039                           4.427327
# CPIJ007892-RA                           5.654220                           5.553004                           5.435790
# CPIJ007894-RA                           3.875605                           4.096226                           3.591727

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
plot(FDP[,1], col=colors()[137], xlab="log2DESeqcounts", ylab="Density", xlim=c(0,16), ylim=c(0,0.40),type="l", lwd=4)
lines(FDP[,2],col=colors()[59], lwd=4)
lines(FDP[,3],col=colors()[101], lwd=4)
lines(FDP[,4],col=colors()[30], lwd=4)
lines(FDP[,5],col=colors()[125], lwd=4)
lines(FDP[,6],col=colors()[128], lwd=4)
legend(legend=c("Field1","Field2","Field3","Lab1","Lab2","Lab3"),x="topright", fill = colors()[c(137,  59,101,30,125,128)])

boxplot(count_matrix_norm_rlog_DESeqData, col=colors()[c(137,  59,101,30,125,128)], main="Normalized Rlog Counts", names=c("Field1","Field2","Field3","Lab1","Lab2","Lab3"))
FDP_DESeq[,1]$call<-"log2counts Densities"
plot(FDP_DESeq[,1], col=colors()[137], xlab="log2DESeqcounts", ylab="Density", xlim=c(0,16), ylim=c(0,0.40),type="l", lwd=4)
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
# log2 fold change (MLE): condition Lab vs Field 
# Wald test p-value: condition Lab vs Field 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE       stat       pvalue        padj
# <numeric>      <numeric> <numeric>  <numeric>    <numeric>   <numeric>
#   CPIJ012141-RA    15.26762    -0.57606372 0.4512284 -1.2766566 0.2017235025 0.464245711
# CPIJ012143-RA 10944.80204     1.44534476 0.3772841  3.8309190 0.0001276655 0.003684142
# CPIJ007879-RA    13.60457    -0.25388412 0.4701962 -0.5399536 0.5892290337 0.802197317
# CPIJ007889-RA    19.38574     0.76922454 0.3924551  1.9600320 0.0499920527 0.204941691
# CPIJ007892-RA    46.82105     0.06297557 0.3334994  0.1888326 0.8502240013 0.945316261
# CPIJ007894-RA    15.38480    -0.30455449 0.4549231 -0.6694637 0.5031997521 0.741688195

head(results(et_DESeq)[order(results(et_DESeq)$padj),])
# log2 fold change (MLE): condition Lab vs Field 
# Wald test p-value: condition Lab vs Field 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE       stat       pvalue         padj
# <numeric>      <numeric> <numeric>  <numeric>    <numeric>    <numeric>
#   CPIJ013879-RA 361.65043      -4.106832 0.3040383 -13.507613 1.410200e-41 6.740754e-38
# CPIJ003554-RA 102.17859       2.709234 0.3479642   7.785955 6.918852e-15 1.209920e-11
# CPIJ016316-RA  42.37781       2.531877 0.3256776   7.774182 7.593639e-15 1.209920e-11
# CPIJ017875-RA  46.55938      -5.216048 0.6939778  -7.516160 5.640838e-14 6.740801e-11
# CPIJ008801-RA 116.08847       2.580591 0.3709662   6.956405 3.490665e-12 3.337076e-09
# CPIJ017845-RA 155.57356       2.653391 0.3859783   6.874456 6.222696e-12 4.957414e-09


# Outlier's analysis. For outlier detection, we should determine which of all genes has its pvalue as NA
table(is.na(results(et_DESeq)$pvalue))
# FALSE 
# 4780 

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

# We consider DE genes those that have a adjusted pavalue lower than 0.05
table(de_DESeq <-decide(dss, p.adj=TRUE, value=0.05,FC_teste = 1,"Field","Lab" ))
# -1    0    1 
# 159 4432  189 

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
