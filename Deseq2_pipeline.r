#Deseq2
#Install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
install.packages("BiocVersion")
library('BiocManager')
install.packages("BiocManager")
BiocManager::install("DESeq2")
library('DESeq2')
library("RColorBrewer")
library("gplots")
library("ggplot2")
library( "genefilter" )
# Grep "counts" format file
sampleFiles <- grep("counts",list.files(directory),value=TRUE)
sampleFiles
# Set the group with control group and treatment group
sampleCondition<- c("Control","Control","Control","Treat","Treat","Treat" )
sampleTable<-data.frame(samTleName= sampleFiles, fileName=sampleFiles, condition=sampleCondition)
#Conduct DESeq package
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
# Set the p value lower than 0.05
resSig <- res[ which(res$padj < 0.05 ), ]
view(res)
# output DESeq
write.csv(as.data.frame(resSig), file=paste(directory,"MLAF9-25um.resSig.csv"))
MLAF925um<- read.csv(paste(directory,"MLAF9-25um.resSig.csv"))
summary(res)
#normalize counts result
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
write.csv(as.data.frame(normalized_counts), file=paste(directory,"normalized_counts.csv"))
normalized_counts_output <- read.csv(paste(directory,"normalized_counts.csv"))
# draw a histogram with p value
hist(res$pvalue, breaks=20, col="white") 
# add linear regression and draw out scatter plot
aa <- c(1,2,3,4,5,6,7,8,9,10)
rld <- rlog( dds )
assay(rld)[ 1:4 , 1:4]
plot( log2( 1+counts(dds, normalized=TRUE)[,1:3] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 7:8], col="#00000020", pch=20, cex=0.3 ) + abline(lsfit(1:10,aa))
# draw heatmap
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
write.csv(as.data.frame(sampleDistMatrix), file=paste(directory,"sampleDistMatrix.csv"))
sampleDistMatrixoutput <- read.csv(paste(directory,"sampleDistMatrix.csv"))
colours = colorRampPalette( rev(brewer.pal(9, "YlOrRd")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours,margins = c(12,12),cexRow = 1,cexCol = 0.8)
# draw MA-plot
plotMA(dds,ylim=c(-3,3),main="DESeq2")
d <- plotCounts(dds, gene="Myc", intgroup="condition", returnData=TRUE)
plotdata <- plotCounts(dds, gene = "Myc", intgroup=c("condition"),returnData = T)
library(ggplot2)
ggplot(plotdata,aes(x=condition,y=count,col=condition))+
  geom_boxplot()+
  theme_bw()
# adjust logFC and plot MA plot
contrast <- c("condition", "Control", "Treat")
dd1 <- results(dds, contrast=contrast, alpha=0.01)
plotMA(dd1, ylim=c(-2.5,2.5))
write.csv(as.data.frame(dd1), file=paste(directory,"dd1.csv"))
dd1_output<- read.csv(paste(directory,"dd1.csv"))
# plot Volcano plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2.5)))
with(subset(res, pvalue<.3 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, pvalue<.1 & abs(log2FoldChange)>0.34), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
vsdata <- varianceStabilizingTransformation(dds, blind=FALSE)
plotPCA(vsdata)
# plot 25 different gene with heatmap
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 25 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
