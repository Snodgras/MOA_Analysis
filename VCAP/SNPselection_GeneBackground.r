#!/usr/bin/env Rscript

# load the module
# r-devtools
# r-tidyr
# r-dplyr
# module load r-devtools r-tidyr r-dplyr

## uncomment these lines when running on the command line
args <- commandArgs(T)
tag <- tools::file_path_sans_ext(args[3])
MOA.snps <- read.table(args[1], sep="\t", stringsAsFactors=FALSE,header=FALSE)
geneBKGD.snps<-read.table(args[2], sep="\t", stringsAsFactors=FALSE,header=FALSE)

library(dplyr)
library(tidyr)

print("This is for selecting a genic background")

#Make the matching do AF bins of 0.05 and distances binned like <500, 500-1000, 1000-5000, >5000 
#make new (binning) variables in the full dataset and then split on that

#creates a function to return a split df based on AF steps
bin_by_AF.10<-function(df){
  df<-df %>% mutate(AFbin = case_when(V4 <= 0.10 ~ "AF.0-10",
                                      V4 > 0.10 & V4 <= 0.20 ~ "AF.10-20",
                                      V4 > 0.20 & V4 <= 0.30 ~ "AF.20-30",
                                      V4 > 0.30 & V4 <= 0.40 ~ "AF.30-40",
                                      V4 > 0.40 & V4 <= 0.50 ~ "AF.40-50",
                                      V4 > 0.50 & V4 <= 0.60 ~ "AF.50-60",
                                      V4 > 0.60 & V4 <= 0.70 ~ "AF.60-70",
                                      V4 > 0.70 & V4 <= 0.80 ~ "AF.70-80",
                                      V4 > 0.80 & V4 <= 0.90 ~ "AF.80-90",
                                      V4 > 0.90 & V4 <= 1 ~ "AF.90-100"),
                    Distbin = case_when(V5 <= 500 & V5 >= -500 ~ "dist.0-500",
                                        V5 <= 1000 & V5 > 500 ~ "dist.500-1000",
                                        V5 >= -1000 & V5 < -500 ~ "dist.500-1000",
                                        V5 <= 5000 & V5 > 1000 ~ "dist.1000-5000",
                                        V5 >= -5000 & V5 < -1000 ~ "dist.1000-5000",
                                        V5 > 5000 | V5 < -5000 ~ "dist.5000-max"))
  return(df)
}

bin_by_AF.05<-function(df){
  df<-df %>% mutate(AFbin = case_when(V4 >= 0.05 & V4 < 0.10 ~ "AF.05-10",
                                      V4 >= 0.10 & V4 < 0.15 ~ "AF.10-15",
                                      V4 >= 0.15 & V4 < 0.20 ~ "AF.15-20",
                                      V4 >= 0.20 & V4 < 0.25 ~ "AF.20-25",
                                      V4 >= 0.25 & V4 < 0.30 ~ "AF.25-30",
                                      V4 >= 0.30 & V4 < 0.35 ~ "AF.30-35",
                                      V4 >= 0.35 & V4 < 0.40 ~ "AF.35-40",
                                      V4 >= 0.40 & V4 < 0.45 ~ "AF.40-45",
                                      V4 >= 0.45 & V4 < 0.50 ~ "AF.45-50",
                                      V4 >= 0.50 & V4 < 0.55 ~ "AF.50-55",
                                      V4 >= 0.55 & V4 < 0.60 ~ "AF.55-60",
                                      V4 >= 0.60 & V4 < 0.65 ~ "AF.60-65",
                                      V4 >= 0.65 & V4 < 0.70 ~ "AF.65-70",
                                      V4 >= 0.70 & V4 < 0.75 ~ "AF.70-75",
                                      V4 >= 0.75 & V4 < 0.80 ~ "AF.75-80",
                                      V4 >= 0.80 & V4 < 0.85 ~ "AF.80-85",
                                      V4 >= 0.85 & V4 < 0.90 ~ "AF.85-90",
                                      V4 >= 0.90 & V4 < 0.95 ~ "AF.90-95",
                                      V4 >= 0.95 & V4 <= 1 ~ "AF.95-100"),
                    Distbin = case_when(V5 <= 500 & V5 >= -500 ~ "dist.0-500",
                                        V5 <= 1000 & V5 > 500 ~ "dist.500-1000",
                                        V5 >= -1000 & V5 < -500 ~ "dist.500-1000",
                                        V5 <= 5000 & V5 > 1000 ~ "dist.1000-5000",
                                        V5 >= -5000 & V5 < -1000 ~ "dist.1000-5000",
                                        V5 > 5000 | V5 < -5000 ~ "dist.5000-max"))
  return(df)
}

MOA.snps<-bin_by_AF.10(MOA.snps) #bin by whichever AF step you want to use
geneBKGD.snps<-bin_by_AF.10(geneBKGD.snps) %>% filter(V5 == 0) #this won't work because it's ANY genes, not genes NEAR MOA peaks, unless the input file is specifically just SNPs in genes near moa peaks


peaks<-group_by(MOA.snps, V10) #check to make sure this is what the column 6 will be called
MOA.distrib.SNPs<-sample_n(peaks,1) #test this to make sure the output is more than just 1 peak

geneBKGD_AF_splitlist<-split(geneBKGD.snps, geneBKGD.snps$AFbin) #split by the AF binning variable
MOA_AF_splitlist<-MOA.distrib.SNPs %>% split(., .$AFbin)

geneBKGD.distrib.SNPs<-tibble(V1=NA, V2=NA, V3=NA,V4=NA,V5=NA,AFbin=NA,Distbin=NA)

#don't need to match on distance to nearest gene
for(AF in names(MOA_AF_splitlist)){ #goes through each AF bin
  s<-paste("^",AF,"$",sep="") #gets name to exact match af bin
  t<-grep(s,names(MOA_AF_splitlist)) #finds exact match for AF in MOA dataframe list
  b<-grep(s,names(geneBKGD_AF_splitlist))  #finds the exact match for AF in background dataframe list
  n<-MOA_AF_splitlist[[t]] %>% nrow() #number of MOA snps that need a match
  if(nrow(geneBKGD_AF_splitlist[[b]] != 0 & n != 0)){ #as long as there's 1 background SNP and MOA snp in that set of bins
    geneBKGD_snp<-sample_n(geneBKGD_AF_splitlist[[b]], n , replace = FALSE) #sample the number of MOA snps in that AF bin
    geneBKGD.distrib.SNPs<-add_row(geneBKGD.distrib.SNPs, V1=geneBKGD_snp$V1,V2=geneBKGD_snp$V2,V3=geneBKGD_snp$V3,V4=geneBKGD_snp$V4,V5=geneBKGD_snp$V5,AFbin=geneBKGD_snp$AFbin,Distbin=geneBKGD_snp$Distbin)
    #write all that info to the background distribution table
  }
}

geneBKGD.distrib.SNPs<-na.omit(geneBKGD.distrib.SNPs) #remove initializing row
write.table(MOA.distrib.SNPs, file=paste0(tag, "MOA_distrib.SNPs.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(geneBKGD.distrib.SNPs, file=paste0(tag, "geneBKGD_distrib.SNPs.txt"), sep="\t", row.names=FALSE, quote=FALSE)
