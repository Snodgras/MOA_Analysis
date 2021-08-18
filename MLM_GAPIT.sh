#!/bin/bash
#MLM_GAPIT.sh
#This script will run the GAPIT pipeline for a MLM GWAS

module load r/3.5.0-py2-x335hrh # load R module
module load r-gplots
#module load r-tidyverse/1.2.1-py2-r3.5-zvowqw7
#module load r-devtools
#source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
#source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
#library(GAPIT3)

#ml r-gplots r-ape r-scatterplot3d r-multtest

#only first time
#module load r-rcpp
#versionyouwant <-"https://cran.r-project.org/src/contrib/Archive/LDheatmap/LDheatmap_0.99-8.tar.gz"
#install.packages(versionyouwant, repos=NULL, type="source")

R # open R

#! /usr/bin/env Rscript
# install packages (only the first time):
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("multtest")
#biocLite("snpStats")
#list.of.packages<-c(#"gplots","LDheatmap", #not available for R3.6.3"genetics",#"ape","EMMREML",#"scatterplot3d")
#new.packages<-list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

# when prompted: type y or yes to create a personal library and 60 to select appropriate mirror
library(gplots)
input.loc<-"/ptmp/LAS/snodgras/VCAP_test/" #MUST CHANGE!!!

#load needed libraries
#library(multtest)
#library(gplots)
#library(genetics)
#library(ape)
#library(EMMREML)
#library(compiler)
#library(scatterplot3d)
#library(LDheatmap)

# use paste0 to combine the input.locand file name into the file path
# functions can be downloaded from GitHub associated with paper
#source(paste0(input.loc, "GAPIT_functions.R"))
#source(paste0(input.loc, "EMMA_functions.R"))

source(paste0(input.loc, "GAPIT.library.R"))
source(paste0(input.loc, "gapit_functions.txt"))
# read in phenotype data
pheno<-read.table("NAM-ril-phenotypes/uniq_testsubset_NAM_phenos.gapit", head = TRUE)

# read in Kinship Matrix
myKI <-read.table("NAMRILs_projectedSNPs_CentIBS_wholegenomekinship.gapit.txt", head = FALSE)
#read in SNPs as hapmap format
myG<-read.table("NAM-ril-genotypes/NAM_rils_SNPs_allchr.hmp.txt", head=TRUE),
# Y = pheno, # set phenotype file

# Run GAPIT GWAS:
myGAPIT<-GAPIT(
Y = pheno # set phenotype file and select relevant columns (1stcolumn = Taxa, 8thcolumn = DF)
KI=myKI, #kinship user supplied option
PCA.total= 0, # set optimal number of PC to control population structure = 0 (see Table S3)
model = "CMLM", # choose GWAS model (other options include MLM, MLMM, etc)
# next 7 lines used to read in genotype files (GD and GM, from chromosome 1 to 11)
G=myG, #Genotypes
)