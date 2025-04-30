# load packages
pkgs <- c("tidyverse") 

vapply(pkgs, FUN = library, FUN.VALUE = logical(1L), logical.return = TRUE, character.only = TRUE)

# custom function
'%notin%'<-Negate('%in%')

# Runs any project ####
# path
path <- "/path/to/where/accession/files/are/" # on linux

#this script requires you to have registered a project and samples from a study at the European Nucleotide Database ENA https://www.ebi.ac.uk/ena/submit/webin/login
#once you have the project and sample accessions, store them in a folder and provide the path to them at the beginning
#move all the R1 and R2 files to a folder and calculate a md5sum of each of them -> then concatenate all and store these into a file: e.g. "checklist.chk"
#see bash script xx.sh for more details

# the following code is an example
#load md5sums file, remove unnecessary characters in file names, extract sense of read and unique sample id from file names

mdsums <- read.delim(paste(path,"checklist.chk",sep=""), sep="", header = F)%>%
          rename(mdsum=V1, file_name=V2)%>% 
          mutate(file_name=str_remove(file_name, "^\\./"))%>%
          filter(file_name%notin%c("checklist.chk"))%>%
          mutate(id=str_extract(file_name, "^[^_]+"))%>%
          mutate(Temp=str_extract(file_name, "_R\\d(?=_)"))%>%
          mutate(Rsense=str_remove(Temp, "^_"))%>%
          select(-Temp)%>%
          mutate(id=str_remove(id, "^[^-]*-"))%>%
          mutate(id=str_remove(id, "^[^-]*-"))

# pivot file names and corresponding md5sums wide
mdsums_w <- mdsums %>% pivot_wider(id_cols = id, names_from = Rsense, values_from = c(file_name, mdsum)) %>%
            rename(forward_file_name=file_name_R1, reverse_file_name=file_name_R2,
            forward_file_md5=mdsum_R1, reverse_file_md5=mdsum_R2)#%>%

# create library_names and ALIAS fields to relate to the sample accessions at ENA
mdsums_f <- mdsums_w %>%
            mutate(Temp = case_when(
            str_starts(id, "b") ~ "16s amplicon of Plot/sample/assay",
            str_starts(id, "f") ~ "ITS2 amplicon of Plot/sample/assay",
            str_starts(id, "c") ~ "18s amplicons of Plot/sample/assay",
           TRUE ~ "other")) %>%
           filter(Temp %notin% c("other")) %>%
           mutate(
           p_n = str_extract(id, "(?<=-)[0-9]+"),
           library_name = paste(Temp, p_n)
           ) %>%
           select(-Temp) %>%
          { 
              row_count <- nrow(.)
              . %>%
              mutate(
              year = rep("-0000", row_count),
              pref = rep("xyz-", row_count)
              )
          } %>%
          mutate(ALIAS = paste(pref, p_n, year, sep = "")) %>%
          select(-pref, -year, -p_n)

# load ENA sample accessions
acc <- read.delim(paste(path,"Webin-accessions-sample.txt",sep=""), sep="", header = T)

# Combine in one file and add additional variables
fq2 <- left_join(mdsums_f, acc, by="ALIAS")%>%
       rename(sample=ACCESSION)%>%
       select(-c(TYPE, ALIAS))%>%
        { 
        row_count <- nrow(.)
        . %>%
        mutate(study=rep("PRJXXXXXXXX", row_count),
               instrument_model=rep("Illumina MiSeq", row_count),
               library_source=rep("METAGENOMIC", row_count),
               library_selection=rep("PCR", row_count),
               library_strategy=rep("AMPLICON", row_count),
               library_layout=rep("PAIRED", row_count)
               )
        }

#reorder columns and save output
fq2.1 <- fq2[, c(7,8,9,6,10,11,12,13,2,4,3,5,1)]
write_tsv(fq2.1, paste(path,"fastq_md5sum_xyz_0000.tsv"))
