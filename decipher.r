#clean slate
rm(list=ls(all.names = TRUE))
gc()
set.seed(123)

library(DECIPHER); packageVersion("DECIPHER")
library(dada2)
library(Biostrings)

path <- "~/"
ASV_tab.nochim <- readRDS(paste(path,"/asv/ASVBtab_file",sep=""))
dna <- readDNAStringSet(paste(path, "/fa/ASV_ba.fasta" ,sep =""))
load("~/ref/16s/SILVA_SSU_r138_2019.RData") #
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(ASV_tab.nochim)

saveRDS(taxid, paste(path,"/fa/taxa_silva138.rds",sep=""))
quit()
y
