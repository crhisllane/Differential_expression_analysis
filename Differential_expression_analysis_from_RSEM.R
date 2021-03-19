##########################################
## Build count matrix with RSEM counts ##
##########################################
# Define the names of counts files for all samples


library(tximport)
direc = "/home/crhisllane/Doutorado/Suzane/RSEM_tximport";
setwd(direc)
count.files<-c("1C.genes.results", "2C.genes.results", "3C.genes.results", "4L.genes.results", "5L.genes.results", "6L.genes.results")

# importing data
txi.rsem <- tximport(count.files, type = "rsem", txIn = FALSE, txOut = FALSE)

#etapa necessaria devido ao erro:
#-using counts and average transcript lengths from tximport
#-Error: all(lengths > 0) is not TRUE
#consultei essa saída em https://support.bioconductor.org/p/92763/
txi.rsem$length[txi.rsem$length == 0] <- 1

#take a look on data
head(txi.rsem$counts)
# [,1] [,2] [,3] [,4] [,5] [,6]
# CPIJ000001-RA    0    0    0    0    0    0
# CPIJ000002-RA    7    4    0    0    2    0
# CPIJ000003-RA    0    0    0    0    0    0
# CPIJ000004-RA   12   10   19    3   10   16
# CPIJ000005-RA    8    6   14    2    3    2
# CPIJ000006-RA   15   14   32    5    4    7

tail(txi.rsem$counts)
# [,1] [,2] [,3] [,4] [,5] [,6]
# CPIJ040893-RA    0    0    0    0    0    0
# CPIJ040894-RA    1    0    0    0    0    0
# CPIJ040895-RA   11   13   23   15   16    6
# CPIJ040896-RA    1    1    0    0    1    1
# CPIJ040897-RA    5    2    8    1    3    0
# CPIJ040898-RA    0    1    0    0    0    0

#how many trasncript do we have?
dim(txi.rsem$counts) 
# [1] 19866     6

#build a dataframe with the sample information
sample_data<-data.frame(condition=c(rep("Field", times=3), rep("Lab", times=3)))

# We must ensure that all genes have 10 counts in at least one condition
summary(rowSums(txi.rsem$counts[,1:3]))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.0      0.0      3.0    183.7     20.0 497500.0 

summary(rowSums(txi.rsem$counts[,4:6]))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.0      0.0      1.0     92.1      6.1 906200.0 


# build a index of genes with good counts
#in this case good counts are defined as those satisfying count >= 10
index_good_counts<-sapply(1:nrow(txi.rsem$counts), function(x){
  idx<-all(txi.rsem$counts[x,1:3] > rep(9, 3)) | all(txi.rsem$counts[x,4:6] > rep(9,3)) 
  return(idx)
})
table(index_good_counts)
# index_good_counts
# FALSE  TRUE 
# 17099  2767


# keep those genes that have good counts
txi.rsem$counts<-txi.rsem$counts[index_good_counts,]
dim(txi.rsem$counts)
# [1] 2767    6
txi.rsem$length <-txi.rsem$length[index_good_counts,]
dim(txi.rsem$length)
# [1] 2767    6
txi.rsem$abundance<-txi.rsem$abundance[index_good_counts,]
dim(txi.rsem$abundance)
# [1] 2767    6



# save filtered count matrix, gene data and sample data 
dados<-list(txi.rsem, sample_data)
save(dados, file="filteredCountMatrix.RData", compress="xz")

library(DESeq2)
rownames(sample_data) <- colnames(txi.rsem$counts)
y_DESeq<-DESeqDataSetFromTximport(txi.rsem, colData=sample_data, design= ~condition)
# using counts and average transcript lengths from tximport

#########################################
#Test de DE
#########################################

y_DESeq <- DESeq(y_DESeq)
et_DESeq<-nbinomWaldTest(y_DESeq)
results.DESeq <- results(et_DESeq)
head(results.DESeq)
# log2 fold change (MLE): condition Lab vs Field 
# Wald test p-value: condition Lab vs Field 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE       stat     pvalue      padj
# <numeric>      <numeric> <numeric>  <numeric>  <numeric> <numeric>
#   CPIJ000004-RA  11.75497     0.74534919 0.5754262  1.2952995 0.19521694 0.4389592
# CPIJ000006-RA  10.40816    -0.63617394 0.5477474 -1.1614367 0.24546433 0.4972180
# CPIJ000011-RA  22.76117     0.11657681 0.4087117  0.2852299 0.77546800 0.8854523
# CPIJ000016-RA  47.12154     0.04167215 0.3927520  0.1061030 0.91550066 0.9616811
# CPIJ000019-RA  18.09786     0.47201580 0.4399821  1.0728067 0.28335785 0.5362867
# CPIJ000023-RA  11.57636    -1.20361965 0.5762020 -2.0888850 0.03671807 0.1704680

head(results(et_DESeq)[order(results(et_DESeq)$padj),])
# log2 fold change (MLE): condition Lab vs Field 
# Wald test p-value: condition Lab vs Field 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat       pvalue         padj
# <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
#   CPIJ015947-RA  17.54917       5.415348 0.8535047  6.344837 2.226608e-10 6.161026e-07
# CPIJ003554-RA  48.74001       3.136093 0.5048606  6.211799 5.238141e-10 7.246969e-07
# CPIJ009101-RA  33.88239      -5.012639 0.8264257 -6.065445 1.315890e-09 1.068173e-06
# CPIJ018582-RA  24.60324      -7.447073 1.2330236 -6.039684 1.544161e-09 1.068173e-06
# CPIJ008785-RA  85.06854       2.385044 0.4023482  5.927811 3.069989e-09 1.698932e-06
# CPIJ018557-RA 583.25766      -1.658998 0.2894787 -5.730985 9.984896e-09 4.604701e-06



#MOSTRAR ISSO A ANTONIO
table(is.na(results(et_DESeq)$pvalue))
#Sem o corte anterior ficava assim:
# FALSE  TRUE 
# 13285  6581
#com o corte anterior:
# FALSE 
# 2767
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
# 122 2553   92 

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

