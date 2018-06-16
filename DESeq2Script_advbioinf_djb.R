### This is an R script that implements the program DESeq2 for gene expression analysis.
### Much more information on the program and specific function (particularly for checking quality) 
### can be found here: https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
### The elements of this script were written by Melissa Pespeni and Daniel Barshis.

#Only need to do this the first time to install the package
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

#setwd("~/dansstuff/Projeks/ODU/CoursesTaught/18Sp_AdvBioinf/data/counts")  # The drag and drop from finder works in R, too.
setwd("~/ODU_MS/ODU_MS/MolecularWork/RNASeq_Analysis")

library(DESeq2)

#useful functions
#head() - prints out the top 6 lines
#dim() - prints the dimensions of a variable
#nrow() - returns the number of rows in a vector or matrix
# ?[functionName] - opens documentation describing the function

#read in your data to make counts table
countsTable <- read.delim('ALLLANESFullCounts_summed.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

#read in a table with conditions of each individual (e.g. "VA" and "RI")  There should be the same number of conditions described as there are samples in your data file, and in the same order.
#NOTE: It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. DESeq2 will not make guesses as to which column of the count matrix belongs to which row of the column data, these must be provided to DESeq2 already in consistent order.
#Here's an example
# sample		origin	symstate		temp
# RI_B_06_18		RI	brown	18
# RI_B_07_14		RI	brown	14
# RI_B_07_18		RI	brown	18
# RI_B_07_22		RI	brown	22
# RI_W_06_18		RI	white	18
# RI_W_07_14		RI	white	14
# RI_W_07_18		RI	white	18
# RI_W_07_22		RI	white	22
# VA_B_06_18		VA	brown	18
# VA_B_07_14		VA	brown	14
# VA_B_07_18		VA	brown	18
# VA_B_07_22		VA	brown	22
# VA_W_06_18		VA	white	18
# VA_W_07_14		VA	white	14
# VA_W_07_18		VA	white	18
# VA_W_07_22		VA	white	22
conds <- read.delim('ALLLANES_conditions.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
conds$temp<-factor(conds$temp)

#make count data sets
dds <- DESeqDataSetFromMatrix(countData=countsTable, colData=conds, design=~ origin + symstate + temp)
dds <- DESeq(dds)

dim(dds)
#prefilter to keep only rows with at least 10 counts to reduce memory consumption
keep <- rowSums(counts(dds)) >= 10
dds10 <- dds[keep,]
dim(dds10)
keep <- rowSums(counts(dds)) >= 5
dds05 <- dds[keep,]
dim(dds05)

###Figure out which contrast you want to examine (i.e. which two groups do you want to compare)
res <- results(dds)
head(res)
summary(res)
resSymState <- results(dds, contrast=c("symstate", "brown", "white"))
head(resSymState)
summary(resSymState)
resOrigin <- results(dds, contrast=c("origin", "VA", "RI"))
head(resOrigin)
summary(resOrigin)
resOrdered <- res[order(res$padj),]
head(resOrdered)
#count the number of significantly differentially expressed genes
sum(res$padj < 0.3, na.rm =T)
sum(res$padj < 0.2, na.rm =T)
sum(res$padj < 0.1, na.rm =T)
sum(res$pvalue < 0.05, na.rm =T)
											 
#filter for contigs with average(baseMean) >5
res5<-res[res$baseMean>5, ]
dim(res)
dim(res5)  # number of genes that have >5 counts

#p-value readjustment Benjamini and Hochberg after >5 filtering
res5$padj <- p.adjust(res5$pvalue, method="BH")

sum(res5$padj < 0.3, na.rm =T)
sum(res5$padj < 0.2, na.rm =T)
sum(res5$padj < 0.1, na.rm =T)

#### Now filter for average counts and variance

# Make a counts table that is scaled by the size factors
temp = t(sizeFactors(dds))
sizematrix<-matrix(data=temp, nrow=nrow(countsTable), ncol=ncol(temp), byrow=TRUE)
scaledcounts = countsTable/sizematrix
head(scaledcounts)

#building heat map data
head(scaledcounts)
genes4heatmap<-res5[res5$pvalue <0.05 & !is.na(res5$pvalue),]
names(genes4heatmap)
head(genes4heatmap)
dim(genes4heatmap)
data4heatmap<-scaledcounts[row.names(scaledcounts)%in%row.names(genes4heatmap),]
dim(data4heatmap)
head(data4heatmap)

temp = as.matrix(rowMeans(data4heatmap))
head(temp)
scaledmatrix<-matrix(data=temp, nrow=nrow(data4heatmap), ncol=ncol(data4heatmap), byrow=FALSE)
data4heatmapscaled = data4heatmap/scaledmatrix
head(data4heatmapscaled)

dim(data4heatmapscaled)

library(gplots)
pairs.breaks <- seq(0, 3.0, by=0.1)
length(pairs.breaks)
mycol <- colorpanel(n=30, low="black", high="yellow") 

pdf(file="XXXXXX_byrow.pdf",7,7)
heatmap.2(data.matrix(data4heatmapscaled), Rowv=T, Colv=F, dendrogram = c("row"), scale="none", keysize=1, breaks=pairs.breaks, col=mycol, trace = "none", symkey = F, density.info = "density", colsep=c(24), sepcolor=c("white"), sepwidth=c(.1,.1), margins=c(10,10), labRow=F)
dev.off()

pdf(file="XXXXXX_bycolumn.pdf",7,7)
heatmap.2(data.matrix(data4heatmapscaled), Rowv=T, Colv=T, dendrogram = c("col"), scale="none", keysize=1, breaks=pairs.breaks, col=mycol, trace = "none", symkey = F, density.info = "density", colsep=c(24), sepcolor=c("white"), sepwidth=c(.1,.1), margins=c(10,10), labRow=F)
dev.off()
