#!/bin/bash
# @ job_name = contam
# @ initialdir = ~/HuConTest
# @ total_tasks = 1
# @ cpus_per_task = 1
# @ wall_clock_limit = 1:00:00
# @ class = normal
# @ array = 1

## these modules might be different depending on the system
module load bcftools/1.18 --auto
module load R/4.3.0 --auto

## this should be executed directly in the directory of the HuConTest script
# cd HuConTest/

# CAUTION 1: Different versions of SAMTOOLS work differently. (Same might be true for R, but the commands are fairly standard in this case.)
# For example, SAMTOOLS 1.0 does not take CRAM files, so it will not give any output. SAMTOOLS 1.9 can use such files. Only these two versions were tested.
# It is advisable to use BAM files, if possible, since this goes faster.
## IMPORTANT UPGRADE on this matter: samtools mpileup is deprecated in more recent versions, instead bcftools mpileup is used now.

# The script will take a comma-separated list of THREE arguments.
# The FIRST is the path to the reference genome to use. Only paths including the strings "hg19" or "hg38" are allowed, and they should point to exactly one of these genomes.
# Obviously, it should be the same genome your sequences were mapped to. 
# The SECOND is the species to test, only three inputs are supported: "pan" (chimpanzee or bonobo), "orang" or "gorilla".
# The THIRD argument is the full path to a BAM file (or CRAM, if necessary), containing sequence alignments to the human reference genome.

# CAUTION 2: Mapped reads should have aligned sequences/inserts of >= 35bp. Spurious alignments (for example, from bacteria) will be identical to the reference genome and appear as "contaminant".

# In case you run many samples as an array, store all file names in a list and get them like that:
ID=$SLURM_ARRAY_TASK_ID
file=$(sed -n ${ID}p /your/file/names.list)

# Or just use the command line, like in this example:
# cd ~/test/
# file="/your/test.bam"

# The script itself is using R to call samtools mpileup for species-specific diagnostic sites, count reads that look like human or the species of choice, and calculate contamination & standard deviation.
# The calculation is: (sum of reads with human allele)/(sum of all reads with any allele).
# The output is FOUR numbers: point estimate of % contaminating reads, point estimate minus one standard deviation across 22 chromosomes, point estimate plus one SD, and the total number of sites covered.
# The script spits it out on the command line.

argum="/your/assembly/hg19.fa,pan,$file"
Rscript --vanilla contamination_test_ape.R ${argum}

exit
