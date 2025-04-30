# Shell functions to archive files in ENA Juan F. Dueñas juanfduenas@proton.me
# this script requires you to have registered a project and samples from a study at the European Nucleotide Database ENA https://www.ebi.ac.uk/ena/submit/webin/login
#
# [optional] go to <path_to_directory> and compress files if they are not compressed
#gzip *.fastq

# [optional] These lines add a prefix to every file name to identify whether the run is 16s, 18s or ITS (the code skips folders) 

for file in *; do
  if [ -f "$file" ]; then
    mv "$file" "B-$file"
  fi
done

for file in *; do
  if [ -f "$file" ]; then
    mv "$file" "F-$file"
  fi
done

#copy (cp) or move (mv) all forward and reverse files to a folder - one must be in the directory where the files are to run this
mv *_R1_001.fastq.gz <path_to_directory>
mv *_R2_001.fastq.gz <path_to_directory>

# create MD5 file for all compressed fastq files
find -type f -exec md5sum "{}" + > checklist.chk

# this will also produce a md5sum of the newly file created and will include a "./" at the beggining of each file. I removed these using an R script

#connect with ENA account through lftp -- You will need a Webin user id and password
lftp webin2.ebi.ac.uk -u Webin-0000

#use command mput to upload
mput *.fastq.gz

