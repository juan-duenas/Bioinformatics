## Bioinformatic processing of runs targeting 18s of Eukaryota with cercozoa specific primers.
# Juan F. Dueñas, PhD ##

# Primer coctail used to produce libraries
#(S616F_Cerco, S616F_Eocer, S963R_Cerco, S947R_Cerco, following Fiore-Donno et al. 2018. Molecular Ecology Resoursces)
# Expected 18s fragment will be about 300 bp within the v4 region of 18s.

## Dada2 bioinformatic pipeline based in tutorial version 1.6 ## 
## https://benjjneb.github.io/dada2/tutorial.html ##

library(dada2)
library(ShortRead)
library(Biostrings)


######## Adapter, primer removal and filtering steps ########
path <- "path/to/fastq.gz"

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#Remove ambiguous bases -> "N" 
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = 4) # on windows, set multithread = FALSE

#sanity check
test <-letterFrequency(sread(readFastq(fnFs[1])), letters = "N")
sum(test)

test <-letterFrequency(sread(readFastq(fnFs.filtN[1])), letters = "N")
sum(test)


# Adapter sequences and orientations
p5 <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
p7 <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"


allOrients_ad <- function(adapter) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(adapter)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
p5.orients.ad <- allOrients_ad(p5) #all possible orientations of p5
p7.orients.ad <- allOrients_ad(p7) #all possible orientations of p7


# find out whether adapters are present among sequences
adapterHits <- function(adapter, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(adapter, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(p5.ForwardReads = sapply(p5.orients.ad, adapterHits, fn = fnFs.filtN[[1]]), # change number within brackets to examine other samples
      p5.ReverseReads = sapply(p5.orients.ad, adapterHits, fn = fnRs.filtN[[1]]), 
      p7.ForwardReads = sapply(p7.orients.ad, adapterHits, fn = fnFs.filtN[[1]]), 
      p7.ReverseReads = sapply(p7.orients.ad, adapterHits, fn = fnRs.filtN[[1]]))

# Reverse compliments of adapters are present

# Check if Cutadapt can be run form R
cutadapt <- "cutadapt"
system2(cutadapt, args = "--version")

# Path to files with adapter free reads
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnf.cut <- file.path(path.cut, basename(fnFs.filtN))
fnr.cut <- file.path(path.cut, basename(fnRs.filtN))

p7.RC <- dada2:::rc(p7)
p5.RC <- dada2:::rc(p5)

# Arguments to trim the reverse-complement of p7 off of R1 (Fowards)
R1.flags.ad <- paste("-b", p7.RC, "--minimum-length 10", "-O 4") #sequences with less than 10 bp were removed
# Arguments to trim the reverse-complement of p5 off of R2 (Reverse)
R2.flags.ad <- paste("-B", p5.RC, "--minimum-length 10", "-O 4") #sequences with less than 10 bp were removed


# Trim adapters with Cutadapt 
# If running on windows, set multiple cores are not supported. So set multiple core option as FALSE

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags.ad, R2.flags.ad, "-n", 1, # -n 1 required to remove adapter only once
                             "-o", fnf.cut[i], "-p", fnr.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

## check for adapter presence in files
rbind(p5.ForwardReads = sapply(p5.orients.ad, adapterHits, fn = fnf.cut[[1]]), 
      p5.ReverseReads = sapply(p5.orients.ad, adapterHits, fn = fnr.cut[[1]]), 
      p7.ForwardReads = sapply(p7.orients.ad, adapterHits, fn = fnf.cut[[1]]), 
      p7.ReverseReads = sapply(p7.orients.ad, adapterHits, fn = fnr.cut[[1]]))

# Now remove primers
S616F_Cerco <- "TTAAAAACGTCGTAGTTG"
S616F_Eocer <- "TTAAAAAGCGCGTAGTTG"
S963R_Cerco <- "CAACTTTCGTTCTTGATTAAA"
S947R_Cerco <- "AAGAAGACATCCTTGGTG"
                
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD1.orients <- allOrients(S616F_Cerco) #all possible orientations of forward prim
FWD2.orients <- allOrients(S616F_Eocer) #all possible orientations of forward prim
REV1.orients <- allOrients(S963R_Cerco) #all possible orientations of reverse prim
REV2.orients <- allOrients(S947R_Cerco) #all possible orientations of reverse prim
FWD1.orients

#find out whether primer sequences are present in sequences
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD1.ForwardReads = sapply(FWD1.orients, primerHits, fn = fnf.cut[[1]]),
      FWD1.ReverseReads = sapply(FWD1.orients, primerHits, fn = fnr.cut[[1]]), 
      REV1.ForwardReads = sapply(REV1.orients, primerHits, fn = fnf.cut[[1]]), 
      REV1.ReverseReads = sapply(REV1.orients, primerHits, fn = fnr.cut[[1]]),
      FWD2.ForwardReads = sapply(FWD2.orients, primerHits, fn = fnf.cut[[1]]),
      FWD2.ReverseReads = sapply(FWD2.orients, primerHits, fn = fnr.cut[[1]]), 
      REV2.ForwardReads = sapply(REV2.orients, primerHits, fn = fnf.cut[[1]]), 
      REV2.ReverseReads = sapply(REV2.orients, primerHits, fn = fnr.cut[[1]]))

# Path to files with primer free reads
path.cut1 <- file.path(path, "cutadaptc")
if(!dir.exists(path.cut1)) dir.create(path.cut1)
fnfs.cut <- file.path(path.cut1, basename(fnf.cut))
fnrs.cut <- file.path(path.cut1, basename(fnr.cut))

# get orientations of interest
FWD2.RC <- dada2:::rc(S616F_Eocer)
REV2.RC <- dada2:::rc(S947R_Cerco)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", S616F_Eocer, "-a", REV2.RC, "--minimum-length 10") #sequences with less than 10 bp were removed
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", S963R_Cerco , "-G", S947R_Cerco , "-A", FWD2.RC, "--minimum-length 10") #sequences with less than 10 bp were removed

# Run Cutadapt 
for(i in seq_along(fnf.cut)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnfs.cut[i], "-p", fnrs.cut[i], # output files
                             fnf.cut[i], fnr.cut[i])) # input files
}

#check for primer presence in files
rbind(FWD1.ForwardReads = sapply(FWD1.orients, primerHits, fn = fnf.cut[[1]]),
      FWD1.ReverseReads = sapply(FWD1.orients, primerHits, fn = fnr.cut[[1]]), 
      REV1.ForwardReads = sapply(REV1.orients, primerHits, fn = fnf.cut[[1]]), 
      REV1.ReverseReads = sapply(REV1.orients, primerHits, fn = fnr.cut[[1]]),
      FWD2.ForwardReads = sapply(FWD2.orients, primerHits, fn = fnf.cut[[1]]),
      FWD2.ReverseReads = sapply(FWD2.orients, primerHits, fn = fnr.cut[[1]]), 
      REV2.ForwardReads = sapply(REV2.orients, primerHits, fn = fnf.cut[[1]]), 
      REV2.ReverseReads = sapply(REV2.orients, primerHits, fn = fnr.cut[[1]]))

#concatenate files containing primer and adapter free reads into R objects

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut1, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut1, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names for latter:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Quality of samples
library(ggplot2)

#inspect quality base calling of reads within samples. Changing numbers within brackets choses other samples
plotQualityProfile(cutFs[1:4]) + scale_x_continuous(n.breaks = 7) #sequences generally good quality (30-40). Quality starts to decline at about base 250 

plotQualityProfile(cutRs[1:4]) + scale_x_continuous(n.breaks = 7) #sequences worst quality than forward reads. Quality starts to decline at about base 200.

### from here on the process can continue at a computer cluster, moving trimmed files to the cluster and executing the script below ####
#clean slate
rm(list=ls(all.names = TRUE))
gc()

# load libraries
library(dada2)
library(ShortRead)
library(Biostrings)

#paths
path.s <- "/path/to/trimmed_files/"
path.h <- "~"
path.cut1 <-file.path(paste(path.s, "folder/", sep=""))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut1, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut1, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names for latter:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#create folder "filtered" within cutadapt folder
filtFs <- file.path(path.cut1, "filtered", basename(cutFs))
filtRs <- file.path(path.cut1, "filtered", basename(cutRs))

#filter command
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(250,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=2) # specify the # of cores in linux. No more than 2 as this command uses mcparallel

# Use machine learning algorithm to "learn" errors and then see how error rate relates to consensus quality score
errF <- learnErrors(filtFs, multithread = T) 
errR <- learnErrors(filtRs, multithread = T)

# Dereplication - find unique reads for both forward and reverse reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
names(derepFs) <- sample.names
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepRs) <- sample.names

# Find 'real' variants with DADA2 algorithm = otherwise known as 'sample inference'
dadaFs <- dada(derepFs, err = errF, multithread = T)
dadaRs <- dada(derepRs, err = errR, multithread = T)

# see how merging works
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=F)

#Construct sequence table
ASV_tab <- makeSequenceTable(mergers)
dim(ASV_tab)

save.image(paste(path.h,"/image.RData", sep=""))

#### from here the script analyzes the outputs of the previous section and removes chimeras ####
#load results
load("~/path/to/image.RData")


# table tracking the loss of samples during initial filtering steps
head(out)
tail(out) 

# check results of dada
dadaFs[[1]]

# see how merging worked
head(mergers[[1]]) # inspect merger results
summary(mergers[[1]]$nmatch)  


# inspect distribution of ASV lengths - consider removing size outliers if they do not fall within the range size expected
table(nchar(getSequences(ASV_tab))) 
hist(nchar(getSequences(ASV_tab)), xlab="Lenght of reads (bp)", main = "ASV length distribution")

#ASV_tab2 <- ASV_tab[,nchar(colnames(ASV_tab)) %in% 000:000] # optional line to trim outliers 000:000 should indicate the expected range

# remove chimeras
ASV_tab.nochim <- removeBimeraDenovo(ASV_tab, method="consensus", multithread=T, verbose=TRUE)

#check how many reads remained after chimera removal
sum(ASV_tab.nochim)/sum(ASV_tab2)

#check size distribution again
hist(nchar(getSequences(ASV_tab.nochim)), xlab="Lenght of reads (bp)", main = "ASV length distribution")

#new path.h
path.h <- "/path/to/deposit/final/results"

#get FASTA file to be able to blast to NCBIs db just in case taxonomic assignation with RDP classifier is not enough
uniquesToFasta(ASV_tab.nochim, paste(path.h, "/ASV_ce.fasta", sep=""), ids=paste0("ASV", seq(length(getUniques(ASV_tab.nochim)))))

# Save ASV table
saveRDS(ASV_tab.nochim, paste(path.h, "/ASVCtab_file", sep=""))

# inspect the number of reads removed in every step 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(ASV_tab.nochim))
colnames(track) <- c("input", "filtered","F_denoised", "R_denoised", "merged", "nonchim")
rownames(track) <- sample.names

library(dplyr)
library(ggplot2)
# track loss of reads
trackl <- track %>% as_tibble(rownames="ID") %>% # keep rownames as column
  tidyr::pivot_longer(!ID, names_to = "step", values_to = "count") %>% # cast in long format to check how many reads were lost per step 
  dplyr::mutate(step=forcats::fct_relevel(step, c("input", "filtered","F_denoised", "R_denoised", "merged", "nonchim"))) # as factor


#plot
ggplot(data = trackl, aes(x=step, y=count, fill=step)) +
  geom_point()+
  geom_boxplot()+
  labs(x="Bioinfo Step", y="No. Reads")+
  theme_bw()

# Save all outputs in final image
save.image(paste(path.h,"/new.image.RData", sep=""))