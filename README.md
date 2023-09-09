# HuConTest
A Human Contamination Test for Great Ape Samples (bcftools mpileup version)

## Description
This is a script for testing human contamination in great ape samples. It uses data mapped to the human genome, and counts the overlap with diagnostic positions where the other great apes differ from humans. This can be useful for any type of sample, including tissue, fecal or ancient DNA from from any of the great ape species (Bonobos, Chimpanzees, Gorillas and Orangutans).

More detail is provided in a technical note, published in Genome Biology & Evolution 2021: https://doi.org/10.1093/gbe/evab117. Please cite if you use!

## UPDATE SEPTEMBER 2023
Samtools mpileup was deprecated in more recent versions, instead bcftools mpileup was recommended. This is now implemented in the code. The old version with samtools can still be found [here](https://github.com/kuhlwilm/HuConTest/tree/8638ad7a27bd8fd499760288e4245b07c55acd54). Furthermore, in this version target files (only containing the position) are used instead of bed files.


## Preparation & Requirements
First, you need to have installations of the following programs:

R (>=3.2) (https://www.r-project.org/)

BCFTOOLS (http://www.htslib.org/)

No specific R libraries are required.


Second, you need to download the script and the data. The data consists of tables with the diagnostic sites to be retrieved, and the great ape allele at these sites, for each species complex.
Create a directory from which to run the test script. Retrieve the "contamination_test_ape.R" script from here, and store it in this directory.
If you are working on a cluster, which is likely if you use this kind of sequencing data, you may want to use arrays jobs for processing several samples in parallel. Then consider the "contamination_ape.sh" wrapper as well.

Create a subdirectory named "data".

Retrieve these files from FigShare and store them in the "data" subdirectory : https://dx.doi.org/10.6084/m9.figshare.24112359.

*Note: This is the version with position files instead of bed files.*

Then you should be ready.

## How to apply it
The script will take a comma-separated string of exactly THREE arguments.

The FIRST is the path to the reference genome to use. 

The SECOND is the species to test, only three inputs are supported here: "pan" (chimpanzee or bonobo), "orang" or "gorilla".

The THIRD argument is the full path to a BAM file (or CRAM, if necessary), containing sequence alignments to the human reference genome.

The full input string should resemble this:

```
argum="/your/assembly/hg19.fa,pan,/your/test/file.bam"
```

The script itself is using R and will call bcftools mpileup for species-specific diagnostic sites, count reads that look like human or the species of choice, and calculate contamination & standard deviation across chromosomes.

The calculation is simple: (sum of reads with human allele)/(sum of all reads with any allele).



The output is FOUR numbers: The point estimate of % contaminating reads, point estimate minus one standard deviation across 22 chromosomes, point estimate plus one SD, and the total number of sites covered.

The script spits it out on the command line. If you want it to be written into a separate file, add something like "> your/file.txt" to the command line.


In order to apply the script, just run:

```
Rscript --vanilla contamination_test_ape.R ${argum}
```

## Considerations
- For the reference genome, only paths including the strings "hg19" or "hg38" are allowed, and they should point to exactly one of these versions of the human reference genome (hg19 or hg38/GRCh38). Patches should not matter regarding the coordinates. The reference genome should be a version with "chr" in the chromosome names. Obviously, it has to be the same genome your sequences were mapped to.
- Different versions of SAMTOOLS may work differently. SAMTOOLS 1.0 does not process files with CRAM extension. Files with BAM extension work faster in either version.
- Some filtering of the raw reads is advisable to get appropriate results. Mapped reads should have aligned sequences/inserts of >= 35bp. Spurious alignments (for example, from bacteria in fecal or ancient samples) will be identical to the reference genome and appear as "contaminant".
- If you get very high numbers (40%, 80%), you should check if you are testing the right species. Using the diagnostic sites of orangutans does not work perfectly with gorilla, but gives intermediate estimates. The same is true if your sequencing data comes from another primate species. If you get almost 100%, it might be either a human library, or spurious alignments (see previous point).
- In principle, the test will be sensitive with at least ~500 read pairs mapped to the reference genome, better more. Less than that may result in absence of observations.
