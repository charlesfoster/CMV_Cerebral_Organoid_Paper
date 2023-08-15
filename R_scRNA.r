## Author: Dr Susan Corley through bioinformatics consultaton
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install('tidyverse')


library("tidyverse")
library("Seurat")


sessionInfo()

###

### Read in sample information

metadata_cells<-read.csv("metadata_human_cells.tsv", sep = "\t", row.names=1)

colnames(metadata_cells)

unique(metadata_cells$PredCellType)
unique(metadata_cells$Sample)

### create subsets of the human organoid data
org.65<- metadata_cells %>%filter(Sample=="SC102A1_65d_org1" | Sample=="SC102A1_65d_org2")
dim(org.65)  #organoid at 65d  #10408 17
head(org.65)

## Look at how many of each cell Type in the subsets

org.65.cellType.tab<-table(org.65$PredCellType)
org.65.cellType.tab<-as.data.frame(org.65.cellType.tab)
org.65.cellType.tab<-filter(org.65.cellType.tab,Freq>10 )
org.65.cellType.tab<-org.65.cellType.tab[order(-org.65.cellType.tab$Freq),]
write.table(org.65.cellType.tab, "org.65.cellType.tab.txt", sep="\t", col.names=NA, quote=F)


### create subsets by cell Type in the org.65 cells

EN.org.65<-org.65 %>%filter(PredCellType=="EN")
dim(EN.org.65)

IN.org.65<-org.65 %>%filter(PredCellType=="IN")
dim(IN.org.65)

IPC.org.65<-org.65 %>%filter(PredCellType=="IPC")
dim(IPC.org.65)

RG.org.65<-org.65 %>%filter(PredCellType=="RG")
dim(RG.org.65)

Glyc.org.65<-org.65 %>%filter(PredCellType=="Glyc")
dim(Glyc.org.65)

Astrocyte.org.65<-org.65 %>%filter(PredCellType=="Astrocyte")
dim(Astrocyte.org.65)

Microglia.org.65<-org.65 %>%filter(PredCellType=="Microglia")
dim(Microglia.org.65)


######## Create Seurat object ######

### Read in Cellranger data - 3 files in cellRanger directory matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz

Org.data <- Read10X(data.dir = "/Data/cellRanger")

str(Org.data)

str(Org.data$protein_coding)
seurat_object = CreateSeuratObject(counts = Org.data$protein_coding)

#19826 features across 90576 samples within 1 assay 

seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

seurat_filt <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &percent.mt < 5)
19826 features across 80653 samples within 1 assay 

##### create seurat object for org65

Org.65_seurat<-subset(seurat_filt, cells=org.65$Barcode)
#19826 features across 10202 samples within 1 assay 

#normalisation

Org.65_seurat<- NormalizeData(Org.65_seurat,scale.factor = 1e6) # normalize to CPM

Org.65_seurat <- FindVariableFeatures(Org.65_seurat, selection.method = "vst", nfeatures = 2000)

### select celltypes by barcode
colnames(EN)
EN.org.65_seurat<-subset(Org.65_seurat, cells=EN.org.65$Barcode)

IN.org.65_seurat<-subset(Org.65_seurat, cells=IN.org.65$Barcode)

IPC.org.65_seurat<-subset(Org.65_seurat, cells=IPC.org.65$Barcode)

RG.org.65_seurat<-subset(Org.65_seurat, cells=RG.org.65$Barcode)

Glyc.org.65_seurat<-subset(Org.65_seurat, cells=Glyc.org.65$Barcode)

Astrocyte.org.65_seurat<-subset(Org.65_seurat, cells=Astrocyte.org.65$Barcode)

Microglia.org.65_seurat<-subset(Org.65_seurat, cells=Microglia.org.65$Barcode)

### Get data

EN.org.65.data<-GetAssayData(object = EN.org.65_seurat[["RNA"]], slot = "data") 
EN.org.65.data.df<-as.data.frame(EN.org.65.data[,1:2])
EN.org.65.data.df<-mutate(EN.org.65.data.df, EN.total=rowSums(EN.org.65.data))
head(EN.org.65.data.df)
dim(EN.org.65.data.df)
EN.org.65.data.df$cell_count <- rowSums(EN.org.65.data!=0, na.rm=T)
EN.org.65.data.df$ncol <- ncol(EN.org.65.data)
EN.org.65.data.df$EN.cellProp <- EN.org.65.data.df$cell_count/EN.org.65.data.df$ncol
EN.org.65.data.df$EN.mean <- EN.org.65.data.df$EN.total/EN.org.65.data.df$ncol
head(EN.org.65.data.df)

###
IN.org.65.data<-GetAssayData(object = IN.org.65_seurat[["RNA"]], slot = "data") 
IN.org.65.data.df<-as.data.frame(IN.org.65.data[,1:2])
IN.org.65.data.df<-mutate(IN.org.65.data.df, IN.total=rowSums(IN.org.65.data))
head(IN.org.65.data.df)
dim(IN.org.65.data.df)
IN.org.65.data.df$cell_count <- rowSums(IN.org.65.data!=0, na.rm=T)
IN.org.65.data.df$ncol <- ncol(IN.org.65.data)
IN.org.65.data.df$IN.cellProp <- IN.org.65.data.df$cell_count/IN.org.65.data.df$ncol
IN.org.65.data.df$IN.mean <- IN.org.65.data.df$IN.total/IN.org.65.data.df$ncol
head(IN.org.65.data.df)

###

IPC.org.65.data<-GetAssayData(object = IPC.org.65_seurat[["RNA"]], slot = "data") 
IPC.org.65.data.df<-as.data.frame(IPC.org.65.data[,1:2])
IPC.org.65.data.df<-mutate(IPC.org.65.data.df, IPC.total=rowSums(IPC.org.65.data))
head(IPC.org.65.data.df)
dim(IPC.org.65.data.df)
IPC.org.65.data.df$cell_count <- rowSums(IPC.org.65.data!=0, na.rm=T)
IPC.org.65.data.df$ncol <- ncol(IPC.org.65.data)
IPC.org.65.data.df$IPC.cellProp <- IPC.org.65.data.df$cell_count/IPC.org.65.data.df$ncol
IPC.org.65.data.df$IPC.mean <- IPC.org.65.data.df$IPC.total/IPC.org.65.data.df$ncol
head(IPC.org.65.data.df)

###
rm(RG.org.65.data)
RG.org.65.data<-GetAssayData(object = RG.org.65_seurat[["RNA"]], slot = "data") 
RG.org.65.data.df<-as.data.frame(RG.org.65.data[,1:2])
RG.org.65.data.df<-mutate(RG.org.65.data.df, RG.total=rowSums(RG.org.65.data))
head(RG.org.65.data.df)
dim(RG.org.65.data.df)
RG.org.65.data.df$cell_count <- rowSums(RG.org.65.data!=0, na.rm=T)
RG.org.65.data.df$ncol <- ncol(RG.org.65.data)
RG.org.65.data.df$RG.cellProp <- RG.org.65.data.df$cell_count/RG.org.65.data.df$ncol
RG.org.65.data.df$RG.mean <- RG.org.65.data.df$RG.total/RG.org.65.data.df$ncol
head(RG.org.65.data.df)
#
###
rm(Glyc.org.65.data)
Glyc.org.65.data<-GetAssayData(object = Glyc.org.65_seurat[["RNA"]], slot = "data") 
Glyc.org.65.data.df<-as.data.frame(Glyc.org.65.data[,1:2])
Glyc.org.65.data.df<-mutate(Glyc.org.65.data.df, Glyc.total=rowSums(Glyc.org.65.data))
head(Glyc.org.65.data.df)
dim(Glyc.org.65.data.df)
Glyc.org.65.data.df$cell_count <- rowSums(Glyc.org.65.data!=0, na.rm=T)
Glyc.org.65.data.df$ncol <- ncol(Glyc.org.65.data)
Glyc.org.65.data.df$Glyc.cellProp <- Glyc.org.65.data.df$cell_count/Glyc.org.65.data.df$ncol
Glyc.org.65.data.df$Glyc.mean <- Glyc.org.65.data.df$Glyc.total/Glyc.org.65.data.df$ncol
head(Glyc.org.65.data.df)

###
rm(Astrocyte.org.65.data)
Astrocyte.org.65.data<-GetAssayData(object = Astrocyte.org.65_seurat[["RNA"]], slot = "data") 
Astrocyte.org.65.data.df<-as.data.frame(Astrocyte.org.65.data[,1:2])
Astrocyte.org.65.data.df<-mutate(Astrocyte.org.65.data.df, Astrocyte.total=rowSums(Astrocyte.org.65.data))
head(Astrocyte.org.65.data.df)
dim(Astrocyte.org.65.data.df)
Astrocyte.org.65.data.df$cell_count <- rowSums(Astrocyte.org.65.data!=0, na.rm=T)
Astrocyte.org.65.data.df$ncol <- ncol(Astrocyte.org.65.data)
Astrocyte.org.65.data.df$Astrocyte.cellProp <- Astrocyte.org.65.data.df$cell_count/Astrocyte.org.65.data.df$ncol
Astrocyte.org.65.data.df$Astrocyte.mean <- Astrocyte.org.65.data.df$Astrocyte.total/Astrocyte.org.65.data.df$ncol
head(Astrocyte.org.65.data.df)

###
rm(Microglia.org.65.data)
Microglia.org.65.data<-GetAssayData(object = Microglia.org.65_seurat[["RNA"]], slot = "data") 
Microglia.org.65.data.df<-as.data.frame(Microglia.org.65.data[,1:2])
Microglia.org.65.data.df<-mutate(Microglia.org.65.data.df, Microglia.total=rowSums(Microglia.org.65.data))
head(Microglia.org.65.data.df)
dim(Microglia.org.65.data.df)
Microglia.org.65.data.df$cell_count <- rowSums(Microglia.org.65.data!=0, na.rm=T)
Microglia.org.65.data.df$ncol <- ncol(Microglia.org.65.data)
Microglia.org.65.data.df$Microglia.cellProp <- Microglia.org.65.data.df$cell_count/Microglia.org.65.data.df$ncol
Microglia.org.65.data.df$Microglia.mean <- Microglia.org.65.data.df$Microglia.total/Microglia.org.65.data.df$ncol
head(Microglia.org.65.data.df)

### create signature for org.65
rm(signature_org.65)
signature_org.65<-subset(EN.org.65.data.df, select=c(EN.total, EN.cellProp,EN.mean))
signature_org.65<-cbind(signature_org.65,IN.org.65.data.df[c("IN.total", "IN.cellProp", "IN.mean")])
signature_org.65<-cbind(signature_org.65,IPC.org.65.data.df[c("IPC.total", "IPC.cellProp", "IPC.mean")])
signature_org.65<-cbind(signature_org.65,RG.org.65.data.df[c("RG.total", "RG.cellProp", "RG.mean")])
signature_org.65<-cbind(signature_org.65,Glyc.org.65.data.df[c("Glyc.total", "Glyc.cellProp", "Glyc.mean")])
signature_org.65<-cbind(signature_org.65,Astrocyte.org.65.data.df[c("Astrocyte.total", "Astrocyte.cellProp", "Astrocyte.mean")])
signature_org.65<-cbind(signature_org.65,Microglia.org.65.data.df[c("Microglia.total", "Microglia.cellProp", "Microglia.mean")])

signature_org.65<-signature_org.65 %>% mutate(mean.total = EN.mean + IN.mean + IPC.mean + RG.mean + Glyc.mean + Astrocyte.mean + Microglia.mean)
head(signature_org.65)
dim(signature_org.65) #19826

signature_org.65.filt<-filter(signature_org.65, mean.total>0)
dim(signature_org.65.filt) #16923    22
colnames(signature_org.65.filt)

signature_org.65.filt2<-filter(signature_org.65.filt, EN.mean>1 |IN.mean>1|IPC.mean>1 |RG.mean>1| Glyc.mean>1|Astrocyte.mean>1| Microglia.mean>1)
dim(signature_org.65.filt2) # 5474   15

signature_org.65.filt3<-filter(signature_org.65.filt, EN.mean>1 &IN.mean>1&IPC.mean>1 &RG.mean>1& Glyc.mean>1&Astrocyte.mean>1& Microglia.mean>1)
dim(signature_org.65.filt3) #1576 15

colnames(signature_org.65.filt2)

signature_org.65_1.mat<-signature_org.65.filt[,c(9:15)]
signature_org.65.mat<-signature_org.65.filt2[,c(9:15)]

signature_org.65_2.mat<-signature_org.65.filt3[,c(9:15)]

write.table(signature_org.65.filt, "signature_org.65.txt", sep="\t",col.names=NA, quote=F)

write.table(signature_org.65_2.mat, "signature_org.65_2.mat.txt", sep="\t",col.names=NA, quote=F)










