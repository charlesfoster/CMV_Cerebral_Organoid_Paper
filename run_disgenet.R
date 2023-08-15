library(data.table)
library(disgenet2r)
library(dplyr)

# read in and subset data
dat <- data.table::fread("Merlin.res.df.csv") 
sig <- dat %>% 
    dplyr::filter(FDR < 0.05)

# remove genes with no gene symbol
universe_genes <- unique(dat$Symbol[dat$Symbol != ""])
de_genes <- unique(sig$Symbol[sig$Symbol != ""])

# run enrichment analysis
all <- disease_enrichment(
  entities = de_genes,
  universe = "CUSTOM",
  custom_universe = universe_genes,
  vocabulary = "HGNC",
  verbose = TRUE,
  database = "ALL",
  warnings = TRUE
)

# extract significant associations
enrichment_results <- extract(all) %>% 
  arrange(FDR) %>% 
  filter(FDR < 0.05)

# plot all enriched diseases with a count of at least 25 disease-associated genes in the DE genes list and an FDR of <=0.05
eplot <- plot( all, class = "Enrichment", count = 25,  cutoff = 0.05, nchars=70)

# write plot and results to file
ggplot2::ggsave(eplot, filename= "top_enriched_diseases.pdf", h=11, w=9)
fwrite(extract(all),"disgenet_enrichment_ALL.csv")
fwrite(enrichment_results,"disgenet_enrichment_ALL_significant.csv")
