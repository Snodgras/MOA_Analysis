#!/bin/bash

# SNP selection for MOA and background

MOAfile=$1
BKGDfile=$2

# randomly pick 1 SNP from each peak
# Find AF and distance from closest gene for picked SNPs

cut -f 10 $MOAfile | sort | uniq > 01_peak.list

awk -v OFS='\t' 'NR < 1001 {print $0 >> sub_01_peak.list};
				 NR >= 1001 && NR < 2001 {print $0 >> sub_02_peak.list};
				 NR >= 2001 && NR < 3001 {print $0 >> sub_03_peak.list};
				 NR >= 3001 && NR < 4001 {print $0 >> sub_04_peak.list};
				 NR >= 4001 && NR < 5001 {print $0 >> sub_05_peak.list};
				 NR >= 1001 && NR < 2001 {print $0 >> sub_02_peak.list};' 01_peak.list

awk -v OFS='\t' 'NR < 1001 {print $0 >> sub_0000_peak.list};
				 %for loop 1..1000% NR >= iterator=1'

while read -r line ; do
	grep -w $line $MOAfile > 02_SNPsinpeak #gets all the SNPs for a given peak
	samplesize=$(cat 02_SNPsinpeak | wc -l) #what's the upper limit for a random number
	RNV=$(shuf -i 1-${samplesize} -n 1) #pick a random number, RNV=random number variable 
	awk -v OFS='\t' -v rnv=$(echo RNV) 'NR == rnv {print "snp."$1"."$3,$10,$5,$4}' 02_SNPsinpeak  >> 03_MOA_SNP_distrib.txt #cut the dist col and AF col for picked SNP and add to ref file
done < 01_peak.list #feed in the peak names

#If we need to know exactly which MOA SNP was picked for the distribution, we can also include that info into the MOA_SNP_distrib.txt 

# match those criteria to non-peak SNPs

while read -r line ; do
	dist=$(echo $line | awk -v OFS='\t' '{print $3}' - )
	AF=$(echo $line | awk -v OFS='\t' '{print $4}' -)
	awk -v OFS='\t' -v d=$dist -v af=$AF '$5 >= d-50 && $5 <= d+50 && $4 == af {print $0}' $BKGDfile >> 04_pool
	samplesize=$(cat 04_pool | wc -l) #what's the upper limit for a random number
	RNV=$(shuf -i 1-${samplesize} -n 1) #pick a random number, RNV=random number variable 
	awk -v OFS='\t' -v rnv=$(echo $RNV) 'NR == rnv {print $0}' 04_pool  >> 05_bkgd_SNP_distrib.txt #print picked background SNP record into a background distribution file
	rm pool
done < 03_MOA_SNP_distrib.txt


####
#R script

# load the module
# r-devtools
# r-tidyr
# r-dplyr
# module load r-devtools r-tidyr r-dplyr

## read table from featureCounts output
args <- commandArgs(T)
tag <- args[3]
#Btag <- tools::file_path_sans_ext(args[2])
MOA.SNPs <- read.table(args[1], sep="\t", stringsAsFactors=FALSE,
  header=FALSE)
BKGD.SNPs<-read.table(args[2], sep="\t", stringsAsFactors=FALSE,
  header=FALSE)
library(dplyr)
library(tidyr)

peaks<-group_by(MOA.SNPs, V10) #check to make sure this is what the column 10 will be called
MOA.distrib.SNPs<-sample_n(peaks,1) #test this to make sure the output is more than just 1 peak
#for the 10,000 SNP selection
MOA.distrib.SNPs<-ungroup(MOA.distrib.SNPs) %>% sample_n(., 10000)

AF_splitlist<-split(BKGD.snps, BKGD.snps$V4)
MOA_AF_splitlist<-MOA.distrib.SNPs %>% split(., .$V4)

BKGD.distrib.SNPs<-tibble(V1=NA, V2=NA, V3=NA,V4=NA,V5=NA)

for(AF in names(MOA_AF_splitlist)){ #for all the dataframes split by AF for MOA distribution SNPs
  if(AF %in% names(AF_splitlist)){ #make sure that AF is one of the possible dataframes for BKGD SNPs split by AF
    #print(c(AF, "matched"))
    s<-paste("^",AF,"$",sep = "") #writes exact match for AF
    t<-grep(s,names(MOA_AF_splitlist)) #finds the exact match for AF in dataframe list
    #print(nrow(MOA_AF_splitlist[[t]]))
  for(snp in 1:nrow(MOA_AF_splitlist[[t]])){ #for every snp in the MOA distribution for that AF
      #make a pool of BKGD snps that match AF and +/- 50 bp of the distance to closest gene
    b<-grep(s,names(AF_splitlist))  #finds the exact match for AF in background dataframe list
    pool<-AF_splitlist[[b]] %>% filter(., V5 >= (MOA_AF_splitlist[[t]][snp,5] - 500) & V5 <= (MOA_AF_splitlist[[t]][snp,5] + 500))
    #print(c(nrow(pool)))
      if(nrow(pool) != 0){ #so long as there's a snp in the pool
        bkgd_snp<-sample_n(pool,1) #go ahead and take a random SNP from it
        BKGD.distrib.SNPs<-add_row(BKGD.distrib.SNPs, V1=bkgd_snp$V1,V2=bkgd_snp$V2,V3=bkgd_snp$V3,V4=bkgd_snp$V4,V5=bkgd_snp$V5)
        #write that SNP to an empty tibble
        }
    }
  }
}
BKGD.distrib.SNPs<-na.omit(BKGD.distrib.SNPs)


write.table(MOA.distrib.SNPs, file=paste0(tag, "MOA_distrib.SNPs.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(BKGD.distrib.SNPs, file=paste0(tag, "BKGD_distrib.SNPs.txt"), sep="\t", row.names=FALSE, quote=FALSE)


