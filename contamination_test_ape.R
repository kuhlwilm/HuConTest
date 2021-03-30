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
lifile=c(paste(sourceDir,"data/h_",spec,"_diffR.lst.gz ",sep=""),paste(sourceDir,"data/h38_",spec,"_diff.lst ",sep=""));names(lifile)<-c("hg19","hg38");lifile=lifile[refg]
bedfile=c(paste(sourceDir,"data/h_",spec,"_diffR.bed.gz ",sep=""),paste(sourceDir,"data/h38_",spec,"_diff.bed ",sep=""));names(bedfile)<-c("hg19","hg38");bedfile=bedfile[refg]
fafile=c(infor[1]);names(fafile)<-refg

# get the data through mpileup and overlap with diagnostic sites
cmnd<-paste("samtools mpileup -u -v -R -t DP,DV -f ",fafile," -l <(gunzip -c ",bedfile,") ",infor[3]," 2>/dev/null | zgrep -v '^#' | awk '{print $1,$2,$10,$5,$1}' | sed 's/:/\\t/g' | sed 's/ /:/' | sed 's/ /\\t/g' | awk -v OFS='\\t' '{print $1,$3,$4,$5,$6}' | sort -k 1b,1 | join -j 1 - <(gunzip -c ",lifile,") | awk '$2!=0||$3!=0' ",sep="")
ctab<-system(paste("/bin/bash -c ", shQuote(cmnd),sep=""),intern=T)

if (length(ctab)==0) { print('It seems something went wrong with the bam file or samtools!');q() }

# processing
ctab<-do.call(rbind,strsplit(ctab,split=" "))
ctab[,4]<-gsub(",<X>|,<\\*>","",ctab[,4])
ctab[,5]<-gsub("chr","",ctab[,5])

# likely human contamination
csel<-which(ctab[,6]==ctab[,4]|ctab[,4]%in%c("<X>","<*>"))
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


