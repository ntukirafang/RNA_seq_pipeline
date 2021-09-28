#GSEA and metabolism
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")
library(dplyr)
library(tibble)
library("AnnotationDbi")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library("AnnotationDbi")
library("org.Mm.eg.db")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
keytypes(org.Mm.eg.db)
require(DOSE)
DOSE::dotplot
library(enrichplot)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")
library(pathview)
dd2 <- lfcShrink(dds, contrast=contrast, res=dd1,type="apeglm")
plotMA(dd2, ylim=c(-2,2))
summary(dd1, alpha = 1)
res <- dd2 %>% 
  data.frame() %>% 
  rownames_to_column("gene_id")
ggplot2(data_plot, aes(condition,"ENSG00000206503",fill=condition)) +
  geom_boxplot() +
  geom_point(size=2, alpha=0.5) +
  geom_line(aes(condition=pairinfo), colour="black", linetype="11") +
  xlab("") +
  ylab(paste("Expression of ","ENSG00000206503"))+
  theme_classic()+
  theme(legend.position = "none")
  res <- res %>% 
  data.frame() %>% 
  rownames_to_column("symbol")
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
gene_df <- res %>% 
dplyr::select(symbol,log2FoldChange,ENSEMBL,entrez) %>% 
filter(entrez!="NA") %>% 
distinct(entrez,.keep_all = T)
geneList <- gene_df$log2FoldChange
view(geneList)
names(geneList) = gene_df$entrez
names(geneList) = gene_df$ENSEMBL
geneList = sort(geneList, decreasing = TRUE)
head(geneList)
write.csv(as.data.frame(gene_df), file=paste(directory,"gene_df.csv"))
gene_df_output<- read.csv(paste(directory,"gene_df.csv"))
write.csv(as.data.frame(geneList), file=paste(directory,"geneList.csv"))
geneList_output<- read.csv(paste(directory,"geneList.csv"))
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

gsepathway <- gseKEGG(geneList=geneList, 
             organism = "mmu",
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.99, 
             verbose = TRUE, 
             pAdjustMethod = "none")
dotplot(kk2,showCategory=50,split=".sign")+facet_grid(~.sign)
emapplot(kk2, showCategory = 50)
cnetplot(kk2, showCategory = 20)
dotplot(gse,showCategory=20,split=".sign")+facet_grid(~.sign)
emapplot(gse, showCategory = 50)
cnetplot(gse, showCategory = 20) 
dotplot(gseaKEGG,showCategory=20,split=".sign")+facet_grid(~.sign)
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
pathway.id = "mmu04215"
gseaplot2(gsepathway,color = "red",geneSetID = pathway.id,pvalue_table = T)
pathway.id = "mmu04215"
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "mmu")
library(pathview)
pathway.id = "mmu04146"
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "mmu",
                   kegg.native = F)
