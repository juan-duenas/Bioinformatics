## Bioinformatic processing of runs targeting ITS2 of Fungi (fITS7-ITS4)
# Juan F. Dueñas, PhD ##

# The fragment of interest should overlap with the ITS2 domain (starting in the 5.8s flanking and ending in LSU)  
# Resulting amplicons are expected to be very variable in size: 120-240 bp (see Tedersoo et al. 2015 Mycokeys) 
# Get the fragment size of primers
#width(fITS7)
#[1] 19
#width(ITS4)
#[1] 20


## Bioinfo ##
## Dada2 bioinformatic pipeline for ITS sequences based in tutorial version 1.8 ## 
## https://benjjneb.github.io/dada2/ITS_workflow.html ##

library(dada2)
library(ShortRead)
library(Biostrings)


######## Adapter, primer removal and filtering steps ########
path <- "path/to/fastq.gz"
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#Remove ambiguous bases -> "N" 

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = 4)

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
R1.flags.ad <- paste("-a", p7.RC, "--minimum-length 10", "-O 4", "-j 4") #sequences with less than 10 bp were removed, min overlap is 4, 4 cores
# Arguments to trim the reverse-complement of p5 off of R2 (Reverse)
R2.flags.ad <- paste("-A", p5.RC, "--minimum-length 10", "-O 4", "-j 4") #sequences with less than 10 bp were removed, min overlap is 4, 4 cores


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
fITS7 <- "GTGARTCATCGAATCTTTG"
ITS4 <- "TCCTCCGCTTATTGATATGC"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(fITS7) #all possible orientations of forward prim
REV.orients <- allOrients(ITS4) #all possible orientations of reverse prim
FWD.orients

#find out whether primer sequences are present in sequences
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnf.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnr.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnf.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnr.cut[[1]]))

# Path to files with primer free reads
path.cut1 <- file.path(path, "cutadapt1")
if(!dir.exists(path.cut1)) dir.create(path.cut1)
fnfs.cut <- file.path(path.cut1, basename(fnf.cut))
fnrs.cut <- file.path(path.cut1, basename(fnr.cut))

# get orientations of interest
FWD.RC <- dada2:::rc(fITS7)
REV.RC <- dada2:::rc(ITS4)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", fITS7, "-g", ITS4, "-a", REV.RC, "-a", FWD.RC,  "--minimum-length 10", "-j 4") #sequences with less than 10 bp were removed, min overlap is 4, 4 cores
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", ITS4, "-G", fITS7, "-A", FWD.RC, "-A", REV.RC, "--minimum-length 10", "-j 4") #sequences with less than 10 bp were removed, min overlap is 4, 4 cores

# Run Cutadapt 
for(i in seq_along(fnf.cut)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnfs.cut[i], "-p", fnrs.cut[i], # output files
                             fnf.cut[i], fnr.cut[i])) # input files
}

#check for primer presence in files
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnfs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnrs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnfs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnrs.cut[[1]]))

# All clear
#concatenate files containing primer and adapter free reads into R objects

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut1, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut1, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names for latter:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
saveRDS(sample.names, paste(path, "sample_names.rds", sep=""))


# Quality of samples
library(ggplot2)
#inspect quality base calling of reads within samples. Changing numbers within brackets choses other samples
plotQualityProfile(cutFs[3:4]) + scale_x_continuous(n.breaks = 7) #sequences generally good quality (30-40). Quality starts to decline at about base 240 

plotQualityProfile(cutRs[3:4]) + scale_x_continuous(n.breaks = 7) #sequences worst quality than forward reads. Quality starts to decline at about base 170. 


#### from here on the process can continue at a computer cluster, moving trimmed files to the cluster and executing the script below ####
# test of CPU efficiency. The basic idea is that the more cores used, the more time needed to set coms between cores (overhead) and not processing a task
# thus it is advice to do a performance experiment with 1, 2, 4, 8 CPUs. Total memory needed for 20 samples was 5 Gb
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
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxN=0, maxEE=c(2,2), minLen = 50, 
                     truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=3) # specify the # of cores in linux.

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
#### from here the script only analyzes the outputs of the previous section code ####
#load results
load("~/path/to/image.RData")

#plot errors
plotErrors(errF, nominalQ = TRUE)

# table tracking the loss of samples during initial filtering steps
head(out) # quite a hefty loss of reads
tail(out) 

# check results of dada
dadaFs[[1]] # 

# see how merging works
head(mergers[[1]]) # inspect merger results
summary(mergers[[1]]$nmatch)  

# inspect distribution of ASV lengths - consider removing size outliers if they do not fall within the range size expected
table(nchar(getSequences(ASV_tab))) 
hist(nchar(getSequences(ASV_tab)), xlab="Lenght of reads (bp)", main = "ASV length distribution")

#remove chimeras
ASV_tab.nochim <- removeBimeraDenovo(ASV_tab, method="consensus", multithread=T, verbose=TRUE)

#check how many reads remained after chimera removal
sum(ASV_tab.nochim)/sum(ASV_tab)

#new path.h
path.h <- "/path/to/deposit/final/results"

#get FASTA file to be able to blast to NCBIs db just in case taxonomic assignation with RDP classifier is not enough
uniquesToFasta(ASV_tab.nochim, paste(path.h, "/fa/ASV_fu.fasta", sep=""), ids=paste0("ASV", seq(length(getUniques(ASV_tab.nochim)))))

# Save ASV table
saveRDS(ASV_tab.nochim, paste(path.h, "/asv/ASVFtab_file", sep=""))

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

# Save results in a new image file
save.image(paste(path.h,"/new.image.RData", sep=""))

