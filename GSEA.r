#GSEA and metabolism 
# install necessary package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(dplyr)
library(tibble)
library("AnnotationDbi")
library(clusterProfiler)
library("org.Mm.eg.db")
library("org.Hs.eg.db")
require(DOSE)
DOSE::dotplot
library(enrichplot)
library(pathview)
# adjust deviation 
dd2 <- lfcShrink(dds, contrast=contrast, res=dd1,type="apeglm")
plotMA(dd2, ylim=c(-2,2))
# Deviation analysis
summary(dd1, alpha = 1)
# Filter Na value
res <- res %>% 
  data.frame() %>% 
# Gene Id add to column
rownames_to_column("gene_id")
# plot specific feature between control and treat
ggplot2(data_plot, aes(condition,"ENSG00000206503",fill=condition)) +
  geom_boxplot() +
  geom_point(size=2, alpha=0.5) +
  geom_line(aes(condition=pairinfo), colour="black", linetype="11") +
  xlab("") +
  ylab(paste("Expression of ","ENSG00000206503"))+
  theme_classic()+
  theme(legend.position = "none") 
# add symbol name to column
rownames_to_column("symbol")
# Exchange Gene ID to ensemble name and entrez 
res$ENSEMBL <- mapIds(org.Mm.eg.db,
                     keys=res$symbol,
                     column="ENSEMBL",
                     keytype="SYMBOL",
                     multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=res$symbol,
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")
# Filter off Na value after exchange 
gene_df <- res %>% 
dplyr::select(symbol,log2FoldChange,ENSEMBL,entrez) %>% 
filter(entrez!="NA") %>% 
distinct(entrez,.keep_all = T)
# Log foldchange with gene name
geneList <- gene_df$log2FoldChange
view(geneList)
# Name as entraz and ensemble
names(geneList) = gene_df$entrez
names(geneList) = gene_df$ENSEMBL
# Sorting with gen list
geneList = sort(geneList, decreasing = TRUE)
head(geneList)
# Write ouput result in csv 
write.csv(as.data.frame(gene_df), file=paste(directory,"gene_df.csv"))
gene_df_output<- read.csv(paste(directory,"gene_df.csv"))
write.csv(as.data.frame(geneList), file=paste(directory,"geneList.csv"))
geneList_output<- read.csv(paste(directory,"geneList.csv"))
# Set reference species in Mouse
keytypes(org.Mm.eg.db)
# Run out KEGG pathway result
gseaKEGG <- gseKEGG(geneList     = geneList,
                    organism     = 'mmu',
                    nPerm        = 1000,
                    minGSSize    = 20,
                    pvalueCutoff = 0.1,
                    verbose      = FALSE)
view(gseaKEGG)
geneList_organism = "mmu"
kk2 <- gseKEGG(geneList     = geneList,
               organism     = geneList_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.99,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

geneList_organism = "mmu"
kk2 <- gseKEGG(organism     = geneList_organism,
               geneList     = geneList,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none")
# Run out gene ontology result
kk2 <- gseGO(geneList=geneList, 
             ont ="ALL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none")

gse <- gseGO(geneList=geneList, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none")
# Run out pathway
gsepathway <- gseKEGG(geneList=geneList, 
             organism = "mmu",
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.99, 
             verbose = TRUE, 
             pAdjustMethod = "none")
# dotplot show top 50 pathway(Up and down)
dotplot(kk2,showCategory=50,split=".sign")+facet_grid(~.sign)
# Emission plot show top 50 pathway
emapplot(kk2, showCategory = 50)
# Central plot show top 20
cnetplot(kk2, showCategory = 20)
dotplot(gse,showCategory=20,split=".sign")+facet_grid(~.sign)
emapplot(gse, showCategory = 50)
cnetplot(gse, showCategory = 20) 
dotplot(gseaKEGG,showCategory=20,split=".sign")+facet_grid(~.sign)
# Write dataframe in csv 
gseaKEGG_results <- gseaKEGG@result
write.csv(as.data.frame(gseaKEGG_results), file=paste(directory,"gseaKEGG_results.csv"))
gseaKEGG_output<- read.csv(paste(directory,"gseaKEGG_results.csv"))
kk2_results <- kk2@result
write.csv(as.data.frame(kk2_results), file=paste(directory,"kk2_results.csv"))
kk2_output<- read.csv(paste(directory,"kk2_results.csv"))
gsepathway_results <- gsepathway@result
write.csv(as.data.frame(gsepathway_results), file=paste(directory,"gsepathway_results.csv"))
gsepathway_output<- read.csv(paste(directory,"gsepathway_results.csv"))
gse_results <- gse@result
write.csv(as.data.frame(gse_results), file=paste(directory,"gse_results.csv"))
gse_output<- read.csv(paste(directory,"gse_results.csv"))
# Plot specific pathway GSEA
pathway.id = "mmu04215"
gseaplot2(gsepathway,color = "red",geneSetID = pathway.id,pvalue_table = T)
pathway.id = "mmu04215"
# Show KEGG pathway (upregulation gene and downregulation gene in pathway) 
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "mmu")

pathway.id = "mmu04146"
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "mmu",
                   kegg.native = F)
