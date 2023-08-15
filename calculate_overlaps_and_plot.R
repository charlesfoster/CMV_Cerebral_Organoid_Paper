suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(org.Hs.eg.db))

args <- commandArgs(trailingOnly = TRUE)

if(length(args) <= 3){
  message("Usage: Rscript calculate_and_plot_overlaps.R differential_expression_results.csv <de_column_name> <comparison_column_name> <dataset 1> <dataset 2> ...")
  quit()
}

message("\nReading files...")
### set up variable with path to ASD SFARI genes
datasets <- args[-c(1,2,3)]

### Comparing 
keep = args[3]

interesting.genes <- unique(unlist(lapply(datasets, function(x){
  d <- data.table::fread(x) %>% 
    dplyr::pull(all_of(keep)) %>% 
    toupper(.)
})))

### read in differential expression results: table with DE genes & table with all genes
### note: set the paths based on the relevant supplementary tables with results
all_genes_table <- data.table::fread(args[1])
de_genes_table <- all_genes_table %>% 
  dplyr::filter(FDR < 0.05)

### define gene sets (and numbers) using Ensembl IDs 
de.genes <- unique(toupper(de_genes_table[[args[2]]]))
de.genes <- de.genes[de.genes!=""]
number_de_genes <- length(de.genes)

all.exp.genes <- unique(toupper(all_genes_table[[args[2]]]))
all.exp.genes <- all.exp.genes[all.exp.genes!=""]
number_total_genes <- length(all.exp.genes)

### dplyr::filter ASD genes list to only include those detected at all in the present study
### this defines the background set of genes for hypergeometric/Fisher tests
interesting.genes <- intersect(interesting.genes, all.exp.genes)
number_interesting_genes <- length(interesting.genes)

### find the overlap, and numberoverlapping, between de genes and asd genes
overlap <- intersect(de.genes, interesting.genes)
number_overlapping <- length(overlap)

### Significance tests: 
## lower.tail=FALSE is set to test for over-representation (also need to take 1 away from overlap)
## lower.tail=TRUE is set to test for under-representation
## run with phyper and fisher.test to show equivalence
message("Calculating overlap significance...")

# phyper results: overrepresented
over.results.phyper <- phyper(q = number_overlapping-1, 
                                  m = number_interesting_genes, 
                                  n = number_total_genes - number_interesting_genes, 
                                  k = number_de_genes, 
                                  lower.tail=FALSE)

# fisher results: overrepresented
over.results.fisher <- fisher.test(matrix(c(number_overlapping, 
                                            number_interesting_genes - number_overlapping, 
                                            number_de_genes - number_overlapping, 
                                            number_total_genes - number_interesting_genes - number_de_genes + number_overlapping), 
                                          2, 
                                          2), 
                                   alternative='greater')$p.value

# optional sanity check: should return true if run
# over.results.phyper == over.results.fisher

# phyper results: underrepresented
under.results.phyper <- phyper(q = number_overlapping, 
                              m = number_interesting_genes, 
                              n = number_total_genes - number_interesting_genes, 
                              k = number_de_genes, 
                              lower.tail=TRUE)

# fisher results: underrepresented
under.results.fisher <- fisher.test(matrix(c(number_overlapping, 
                                             number_interesting_genes - number_overlapping, 
                                            number_de_genes - number_overlapping, 
                                            number_total_genes - number_interesting_genes - number_de_genes + number_overlapping), 
                                          2, 
                                          2), 
                                   alternative='less')$p.value

# optional sanity check: should return true if run
# under.results.phyper == under.results.fisher

### Generate results summary table
summary <- data.table::data.table(total_genes = number_total_genes,
                     number_DE = number_de_genes,
                     number_overlapping_background = number_interesting_genes,
                     number_overlapping_DEGs = number_overlapping,
                     p.value.over = over.results.phyper,
                     p.value.under = under.results.phyper,
                     overrepresented = fifelse(over.results.phyper < 0.05, TRUE, FALSE),
                     underrepresented = fifelse(under.results.phyper < 0.05, TRUE, FALSE))


### write to file
fname <- "hypergeometric_results_summary.csv"
data.table::fwrite(summary, fname)
knitr::kable(summary)

overlapping_dataset <- all_genes_table %>% 
  filter(get(args[2]) %in% interesting.genes) 

fname <- "overlapping_genes.csv"
data.table::fwrite(overlapping_dataset, fname)

### plot some of the overlap data
### note: to get the desired plots, some unusual modifications to the clusterProfiler plotting objects/classes is necessary
message("Plotting...")
# select top 25 up and down DE genes (by FDR)
DE_overlap.up.top25 <- de_genes_table %>% 
  dplyr::filter(get(args[2]) %in% interesting.genes) %>% 
  dplyr::filter(logFC > 0) %>% 
  dplyr::arrange(desc(logFC)) %>% 
  dplyr::slice_head(n=25)

DE_overlap.down.top25 <- de_genes_table %>% 
  dplyr::filter(get(args[2]) %in% interesting.genes) %>% 
  dplyr::filter(logFC < 0) %>% 
  dplyr::arrange(logFC) %>% 
  dplyr::slice_head(n=25)

DE_overlap <- rbind(DE_overlap.up.top25,DE_overlap.down.top25)
if(args[2] == toupper("ENSEMBL")){
  # convert symbols using clusterProfiler functions for compatibility
  converted <- suppressMessages(bitr(DE_overlap$Ensembl, 
                      fromType = "ENSEMBL", 
                      toType=c("ENTREZID","SYMBOL"), 
                      OrgDb="org.Hs.eg.db"))

  plot_data <- dplyr::full_join(converted, DE_overlap, by = c("ENSEMBL" = "Ensembl"))

  # do temporary GO enrichment to get an object in the right format for clusterProfiler plotting
  # note: these enrichment results are not used for anything else
  tmp <- enrichGO(plot_data$Entrez, OrgDb="org.Hs.eg.db", pvalueCutoff=1, qvalueCutoff=1)
  ego <- suppressMessages(setReadable(tmp, OrgDb = org.Hs.eg.db))

  # reformatted DE results for plotting: all overlapping with ASD
  reformatted <- data.frame(ID = DE_overlap$Ensembl,
                      Description = DE_overlap$Symbol,
                      GeneRatio = DE_overlap$logFC,
                      BgRatio = rep("1/1", length(DE_overlap$Ensembl)),
                      pvalue = DE_overlap$PValue,
                      p.adjust = DE_overlap$FDR,
                      qvalue = DE_overlap$FDR,
                      geneID = DE_overlap$Entrez,
                      Count = rep(1, length(DE_overlap$Ensembl)))
} else {
  plot_data <- DE_overlap
  tmp <- enrichGO(plot_data$Entrez, OrgDb="org.Hs.eg.db", pvalueCutoff=1, qvalueCutoff=1)
  ego <- suppressMessages(setReadable(tmp, OrgDb = org.Hs.eg.db))
  reformatted <- data.frame(ID = DE_overlap$Symbol,
                      Description = DE_overlap$Symbol,
                      GeneRatio = DE_overlap$logFC,
                      BgRatio = rep("1/1", length(DE_overlap$Symbol)),
                      pvalue = DE_overlap$PValue,
                      p.adjust = DE_overlap$FDR,
                      qvalue = DE_overlap$FDR,
                      geneID = DE_overlap$Entrez,
                      Count = rep(1, length(DE_overlap$Symbol)))
}

rownames(reformatted) <- reformatted$ID
ego@result <- reformatted
ego@result$LFC <- DE_overlap$logFC

plot1 <- barplot(ego, x="LFC",showCategory = nrow(reformatted),xlab="log2 fold change")+
  xlab("log2 fold change")

ggsave("top25Up_top25Down_by_FDR.pdf", plot=plot1, width=6,height=8,units="in")

# reformatted DE results for plotting: upregulated overlapping with ASD
up <- de_genes_table %>% 
  filter(logFC > 0)

up$Symbol[which(up$Symbol == "")] <- up$Ensembl[which(up$Symbol == "")]

reformatted <- data.frame(ID = up$Ensembl,
                    Description = up$Symbol,
                    GeneRatio = up$logFC,
                    BgRatio = rep("1/1", length(up$Ensembl)),
                    pvalue = up$PValue,
                    p.adjust = up$FDR,
                    qvalue = up$FDR,
                    geneID = up$Entrez,
                    Count = rep(1, length(up$Ensembl)))
reformatted$ID[is.na(reformatted$ID)] <- "NA"
rownames(reformatted) <- make.unique(reformatted$ID)
ego@result <- reformatted
ego@result$LFC <- up$logFC

plot2 <- barplot(ego, x="LFC",showCategory = 25) +
  xlab("log2 fold change")

ggsave("de_up_barplot.pdf", plot=plot2, width=6,height=8)

# reformatted DE results for plotting: downregulated overlapping with ASD
down <- de_genes_table %>% 
  filter(logFC < 0)

down$Symbol[which(down$Symbol == "")] <- down$Ensembl[which(down$Symbol == "")]

reformatted <- data.frame(ID = down$Ensembl,
                             Description = down$Symbol,
                             GeneRatio = down$logFC,
                             BgRatio = rep("1/1", length(down$Ensembl)),
                             pvalue = down$PValue,
                             p.adjust = down$FDR,
                             qvalue = down$FDR,
                             geneID = down$Entrez,
                             Count = rep(1, length(down$Ensembl)))
reformatted$ID[is.na(reformatted$ID)] <- "NA"
rownames(reformatted) <- make.unique(reformatted$ID)
ego@result <- reformatted
ego@result$LFC <- down$logFC

plot3 <- barplot(ego, x="LFC",showCategory = 25) +
  xlab("log2 fold change")

ggsave("de_down_barplot.pdf", plot=plot3, width=6,height=8)

# reformatted DE results for plotting: DE top 25 total overlapping with ASD by FDR
all <- de_genes_table
all$Symbol[which(all$Symbol == "")] <- all$Ensembl[which(all$Symbol == "")]
reformatted <- data.frame(ID = all$Ensembl,
                    Description = all$Symbol,
                    GeneRatio = all$logFC,
                    BgRatio = rep("1/1", length(all$Ensembl)),
                    pvalue = all$PValue,
                    p.adjust = all$FDR,
                    qvalue = all$FDR,
                    geneID = all$Entrez,
                    Count = rep(1, length(all$Ensembl)))

ego@result <- reformatted
ego@result$LFC <- all$logFC

plot4 <- barplot(ego, x="LFC",showCategory = 50) +
  xlab("log2 fold change")

ggsave("top50_ranked_by_FDR.pdf", plot=plot4, width=6,height=8,units="in")

message("Done!")

