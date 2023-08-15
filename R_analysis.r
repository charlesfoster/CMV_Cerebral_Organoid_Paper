## Author: Dr Susan Corley through bioinformatics consultaton
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("statmod")
BiocManager::install("DESeq2")
BiocManager::install("gplots")
BiocManager::install("biomaRt")
BiocManager::install("RColorBrewer")
BiocManager::install("tidyverse")
BiocManager::install("pheatmap")

library("limma")
library("edgeR")
library("DESeq2")
library("gplots")
library("biomaRt")
library("RColorBrewer")
library("tidyverse")
library("ggplot2")
library("pheatmap")
library("gplots")



### Read in data 
Data<-read.delim("featureCounts.txt", header=T) 

### remove unnecessary columns from Data
Data.ed<-Data[, -c(2:6)]
rownames(Data.ed)<-Data$Geneid
colnames(Data.ed)

### read in sample information
sample_info <- read.delim("Sample_info.txt", header=T)
head(sample_info)
sample_names<-as.character(sample_info$Sample_name)

colnames(Data.ed)[c(2:9)] <- sample_names
colnames(Data.ed)
dim(Data.ed) #
Data.ed<-Data.ed[rowSums(Data.ed[,2:9])>0,] # remove rows if all values zero
dim(Data.ed) #

## Add in annotations
library(biomaRt)

ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "asia")


listMarts(host="https://www.ensembl.org")


#attributes<-listAttributes(human)
attributes=c("ensembl_gene_id", "entrezgene_id","hgnc_symbol","gene_biotype")

Rseq_GeneID<-Data.ed$Geneid
head(Rseq_GeneID)
Genemap.108<-getBM(attributes, filters="ensembl_gene_id", values=Rseq_GeneID, mart=ensembl)
idx <-match(Rseq_GeneID, Genemap.108$ensembl_gene_id)
Data.ed$Entrez<-Genemap.108$entrezgene [ idx ]
Data.ed$Ensembl<-Genemap.108$ensembl_gene_id [ idx ]
Data.ed$Symbol<-Genemap.108$hgnc_symbol [ idx ]
Data.ed$gene_biotype<-Genemap.108$gene_biotype [ idx] 
colnames(Data.ed)
dim(Data.ed) #46379 13

###
Annot<-Data.ed[,c(1,10:13)]
head(Annot)
Annot[Annot==""] <- NA # if there is a blank insert NA


head(sample_info)
condition=factor(sample_info$Condition, levels=c("Con", "Merlin"))
design<-model.matrix(~condition)
print(design)
colnames(design)<- c("Intercept","Merlin")

################# create DGE objects using edgeR

Rseq_y<-DGEList(counts=Data.ed[,2:9], group=condition, genes=Annot)
dim(Rseq_y) #

### filter to remove lowly expressed genes 

Rseq_keep<-rowSums(edgeR::cpm(Rseq_y)>0.5) >=4
Rseq_y<-Rseq_y[Rseq_keep, keep.lib.size =FALSE]

length(unique(Rseq_y$genes$Symbol))

### normalisation
Rseq_y <-calcNormFactors(Rseq_y, method="TMM")


### estimate dispersion
Rseq_y <- estimateGLMCommonDisp(Rseq_y, design, verbose=T) 
Rseq_y <- estimateGLMTrendedDisp(Rseq_y, design)
Rseq_y <- estimateGLMTagwiseDisp(Rseq_y,design)


alpha=0.05

edgeR_fit<-glmQLFit(Rseq_y, design, robust=T) #Fit a quasi-likelihood negative binomial generalized log-linear model to count data 

print(design)
head( edgeR_fit$coefficients)

Merlin.res<- glmQLFTest(edgeR_fit)
Merlin.res.df<-topTags(Merlin.res, n=nrow(Merlin.res), adjust.method="BH")$table
head(Merlin.res.df)

Merlin_fdr_sig<-Merlin.res.df[Merlin.res.df$FDR<alpha,]
head(Merlin_fdr_sig)
dim(Merlin_fdr_sig) #

### capture log counts
rm(logCPMcounts)
logCPMcounts<-edgeR::cpm(Rseq_y, log=T)
colnames(logCPMcounts)
logCPMcounts.df<-as.data.frame(logCPMcounts)
logCPMcounts.df$Entrez<-Rseq_y$genes$Entrez
logCPMcounts.df$Ensembl<-Rseq_y$genes$Ensembl
logCPMcounts.df$Symbol<-Rseq_y$genes$Symbol
logCPMcounts.df$gene_biotype<-Rseq_y$genes$gene_biotype

### Functional analysis

### KEGG pathways

head(Merlin.res)
Merlin_Kegg<-kegga(Merlin.res,coef=2, geneid=Merlin.res$genes$Entrez, FDR=0.05, species = "Hs", plot=F )
Merlin_Kegg_up<- topKEGG(Merlin_Kegg, sort = "up",number = 50L)
Merlin_Kegg_down<- topKEGG(Merlin_Kegg, sort = "down",number = 50L)

### GO terms
Merlin_GO<-goana(Merlin.res,coef=2, geneid=Merlin.res$genes$Entrez, FDR=0.05, species = "Hs", plot=F )
Merlin_GO_up<- topGO(Merlin_GO, sort = "up",number = 100L,ontology=c("BP", "CC", "MF"))
Merlin_GO_up.sig<- Merlin_GO_up[Merlin_GO_up$P.Up<0.05,]
Merlin_GO_down<- topGO(Merlin_GO, sort = "down",number = 100L,ontology=c("BP", "CC", "MF"))
Merlin_GO_down.sig<- Merlin_GO_down[Merlin_GO_up$P.Down<0.05,]

Merlin_GO_BP_up<- topGO(Merlin_GO, sort = "up",number = 300L,ontology=c("BP"))
Merlin_GO_BP_up.sig<- Merlin_GO_BP_up[Merlin_GO_BP_up$P.Up<0.05]
Merlin_GO_BP_down<- topGO(Merlin_GO, sort = "down",number = 300L,ontology=c("BP"))
Merlin_GO_BP_down.sig<- Merlin_GO_BP_down[Merlin_GO_BP_up$P.Down<0.05]

## prepare data for making kegg heatmaps

logCPMcounts.DE<-logCPMcounts.df%>%dplyr::filter(Ensembl %in% Merlin_fdr_sig$Ensembl)

###
GK <- limma::getGeneKEGGLinks(species.KEGG="hsa")
dim(GK) #35642     2

### Example Pathways of interest

p1<-c("path:hsa04360") #Axon guidance
n1<-paste("Axon_guidance", p1, sep="_")

p2<-c("path:hsa04728") #Dopaminergic synapse
n2<-paste("Dopaminergic_synapse", p2, sep="_")

###KEGG pathways

Axon_guidance <- GK[ which(GK$PathwayID==p1), ] #
dim(Axon_guidance) #182 2

Dopaminergic_synapse <- GK[ which(GK$PathwayID==p2), ] #
dim(Dopaminergic_synapse) #132 2

dat1<- logCPMcounts.DE %>%filter(Entrez %in% Axon_guidance$GeneID) #41 12
dat2<- logCPMcounts.DE %>%filter(Entrez %in% Dopaminergic_synapse$GeneID) #35 12

###  Heatmap using library(pheatmap)

set.seed(123)

matrix<-as.matrix(dat1[,1:8])
rownames(matrix)<-dat1$Symbol
nr=nrow(matrix)
cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
rows.cor <- cor(t(matrix), use = "pairwise.complete.obs", method = "pearson")
dcols = as.dist(1 - cols.cor)
drows = as.dist(1 - rows.cor)

# Plot the heatmap
Axon_pheatmap<-pheatmap(matrix, scale="row", clustering_distance_rows=drows, cluster_row = T,
                        clustering_distance_cols = dcols, cluster_cols = F,
                        clustering_method = "complete",show_rownames = T,
                        fontsize_row=6,fontsize_col=12,cellwidth = 16, cellheight=6, main=paste0(n1," DE"))
                        
save_pheatmap_png <- function(x, filename, width=900, height=1400, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(Huntington.pheatmap, "Axon_pheatmap.png")





