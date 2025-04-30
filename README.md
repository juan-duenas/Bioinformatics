# Bioinformatics
Annotated scripts for different bioinformatic tasks used in different projects

### B_dada_script.r
Script to process amplicon libraries targeting the 16s of Procariotes. The script is based on the *DADA2* tutorial version [1.6](https://benjjneb.github.io/dada2/tutorial.html)
Libraries were generated targeting a region of about 400 bp spanning the V4/V5 domains in 16s using primers [515f and 926r](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-38)

### F_dada_script.r
Script to process amplicon libraries targeting the ITS2 domain (starting in the 5.8s flanking and ending in LSU). The script is based on the *DADA2* tutorial for ITS version [1.8](https://benjjneb.github.io/dada2/ITS_workflow.html)
Libraries were generated targeting a fragment of variable size using primers [fITS7 and ITS4](https://doi.org/10.1111/j.1574-6941.2012.01437.x)

### C_dada_script.r
Script to process amplicon libraries targeting the 18s of Eukaryota with cercozoa specific primers. The script is based on the *DADA2* tutorial version [1.6](https://benjjneb.github.io/dada2/tutorial.html)
Libraries were generated targeting a fragment within the v4 region of 18s using a coctail of primers [S616F_Cerco, S616F_Eocer, S963R_Cerco, S947R_Cerco](https://doi.org/10.1111/1755-0998.12729)

### decipher.r
Example **R** script to taxonomically annotate operational taxonomic units or amplicon sequence variants with the [IDTaxa](https://doi.org/10.1186/s40168-018-0521-5) algorithm. This example uses the Silva 16s database [v. 138.1](https://academic.oup.com/nar/article/41/D1/D590/1069277) as a reference for taxonomic annotation.

### taxid_slurm.sh
Bash script that executes an **R** script remotely by taking advantage of a high performance computer cluster administered throuth [SLURM](https://slurm.schedmd.com/overview.html)

### format_tvs_arch_ENA.r
**R** script to produce a tvs file to archive runs at ENA the European Nucleotide Database [ENA](https://www.ebi.ac.uk/ena/submit/webin/login). The script relates the run files with the accession numbers provided by ENA.

### ena_upload_script.sh
**Shell** commands to archive Illumina run files genereted during a study. The files can be archived at the European Nucleotide Database [ENA](https://www.ebi.ac.uk/ena/submit/webin/login)

### cat_fastq_files_script.sh
**Shell** commands and loop functions to concatenate files repeated during sequencing
