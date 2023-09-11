#!/usr/bin/r
position=(commandArgs(TRUE))
options("scipen"=100)
sourceDir <- getSrcDirectory(function() {})

if(length(position)==0) { cat('This script needs some input! That is a comma-separated string of characters of the form:\nreferencegenome.fa,species,testfile.bam"\nHave a look at the documentation.\n');q() }

infor=unlist(strsplit(as.character(position),split=","))
## infor[1] is the fasta file, infor[2] is the species, infor[3] the test bam file

if (!infor[2]%in%c("pan","orang","gorilla")) { print('Species must be "pan", "orang", or "gorilla"!');q() }
if (length(grep("hg19|hg38|GRCh38",infor[1]))==0) { print('Reference genome must contain "hg19", or "hg38"/"GRCh38"!');q() }
spec=infor[2]
refg=ifelse(length(grep("hg19",infor[1]))==1,"hg19",ifelse(length(grep("hg38|GRCh38",infor[1]))==1,"hg38",""))

# define input files
lifile=c(paste(sourceDir,"data/h_",spec,"_diffR.lst.gz ",sep=""),paste(sourceDir,"data/h38_",spec,"_diff.lst.gz ",sep=""));names(lifile)<-c("hg19","hg38");lifile=lifile[refg]
bedfile=c(paste(sourceDir,"data/h_",spec,"_diffR.pos.gz ",sep=""),paste(sourceDir,"data/h38_",spec,"_diff.pos.gz ",sep=""));names(bedfile)<-c("hg19","hg38");bedfile=bedfile[refg]
fafile=c(infor[1]);names(fafile)<-refg

# get the data through mpileup and overlap with diagnostic sites
cmnd<-paste("bcftools mpileup -d 2000 --ignore-RG -a FORMAT/DP,FORMAT/AD -T ",bedfile," --fasta-ref ",fafile," ",infor[3]," 2>/dev/null | zgrep -v '^#' | awk '{print $1,$2,$10,$5,$1}' | sed 's/:/\\t/g' | sed 's/ /:/' | sed 's/ /\\t/g' | awk -v OFS='\\t' -v FS='\\t' '{print $1,$3,$4,$5,$6}' | sort -k 1b,1 | join -j 1 - <(gunzip -c ",lifile,") ",sep="")
## Note that the tabs are here represented as '\\t' instead of '\t' (to run with R system() ). If you are debugging/testing step by step on the command line, this needs to be corrected.

ctab<-system(paste("/bin/bash -c ", shQuote(cmnd),sep=""),intern=T)

# check if bam file is sorted, otherwise generic error message
if (length(ctab)==0) { if(grep("SO:coordinate",system(paste("samtools view -T ",fafile," -H ",infor[3]," | grep @HD"),intern=T))==F) { print('Is your bam file sorted?');q() } else { print('It seems something went wrong with the bam file (or bcftools)!');q() } }

# processing
ctab<-do.call(rbind,strsplit(ctab,split=" "))
ctab[,4]<-gsub(",<X>|,<\\*>","",ctab[,4])
ctab[,5]<-gsub("chr","",ctab[,5])
ctab<-ctab[which(ctab[,2]>0),,drop=F]

# likely human contamination
csel<-which(ctab[,6]==ctab[,4]|ctab[,4]%in%c("<X>","<*>"))
ctab[,3]<-suppressWarnings(do.call(rbind,strsplit(as.character(ctab[,3]),split=","))[,2])
# warning message about number of columns can be ignored - only the number of reads supporting the second allele matters
csu<-sum(as.numeric(ctab[csel,2])-as.numeric(ctab[csel,3]))
nsu<-sum(as.numeric(ctab[,2]))
cn<-csu/nsu
cis<-c()
for (chr in c(1:22)) { cta<-ctab[which(ctab[,5]==chr),,drop=F];csel<-which(cta[,6]==cta[,4]|cta[,4]%in%c("<X>","<*>")); cis[chr]= sum(as.numeric(cta[csel,2,drop=F])-as.numeric(cta[csel,3,drop=F]))/sum(as.numeric(c(cta[csel,2,drop=F]))) }
sds<-sd(cis,na.rm=T)

# output: Percent contamination; 95% CI, number of sites considered
op<-c(round(as.numeric(c(cn,ifelse((cn-sds)<0,0,cn-sds),ifelse((cn+sds)>1,1,cn+sds)))*100,3))
cat(c(op,as.integer(nrow(ctab)),"\n"))

q()
