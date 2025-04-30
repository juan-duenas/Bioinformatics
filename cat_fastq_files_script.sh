# Shell commands and loop functions to concatenate files repeated during sequencing
# Juan F. Dueñas juanfduenas@proton.me

# first remove prefixes from file names added by Sequencing facility
# remove the sample identifier in the middle of each file

# test #run this to satisfy that everything is ok.
for file in * ; do
    echo mv -v "$file" "${file/_S+([[:digit:]])_/_}"
done

#if it is, remove echo from command and it will rename files as you want.
for file in * ; do
    mv -v "$file" "${file/_S+([[:digit:]])_/_}"
done

# then concatenate files based on names

# get into the folder where you want to place concatenated files - ot should be a folder within the folder containing the fastq.gz files to be concatenated
cd <namefile>

# then run this to check you are getting what you want - you might modify the code below to achieve what you want.
for name in ../*.fastq.gz; do
rnum=${name#*_};
Sample=${name%%_*};
Sample=${Sample#*-};
printf "%s\n" "$Sample" "$rnum";
done

# once the function above has been proofed, then run this function to append files in one.
for name in ../*.fastq.gz; do
    rnum=${name#*_};
    sample=${name%%_*};
    sample=${sample#*-};
    cat "$name" >> "${sample}_$rnum"
done

# Sanity check - I create a temp folder named <filename>. This is where all the concatenated reads go. I unzip the files and check how many total reads there are after the concatenation. I cross check these values with the table provided by the sequencing facility to be sure I am getting all the reads per sample.

cd <namefile>

#unzip
gunzip *

# count number of reads per file and concatenate into a table
for file in *; do
 printf "%s\t" "$file" ;
 expr $(cat "$file" | wc -l) / 4
done

# Once all is checked compress the folder containing all the concatenated files

tar -zcvf xyz.tar.gz <filename>/

