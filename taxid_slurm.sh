#!/bin/bash
#SBATCH -N1 -n4 --mem-per-cpu=1024M --qos=prio -t01:30:00
module add R/4.2.2-foss-2022b
Rscript decipher.r
