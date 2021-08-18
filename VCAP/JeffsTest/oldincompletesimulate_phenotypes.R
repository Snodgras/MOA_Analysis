library("tidyverse")
library()
setwd("Example_Kinships_1SNPperPeak/")
moa_K<-make_K("MOA.SNPs.bed_CentIBS_K") %>% 
  normalize_kinmat()
smoa_K<-svd(moa_K)

back_K<-make_K("BKGD.SNPs.bed_CentIBS_K") %>% 
  normalize_kinmat()
sback_K<-svd(back_K)

rest_K<-make_K("CentIBS_Rest_K") %>% 
  normalize_kinmat()
srest_K<-svd(rest_K)

# Three K matrices plus residuals
cmoa_K = sqrt(smoa_K$d) * t(smoa_K$u) 
cback_K = sqrt(sback_K$d) * t(sback_K$u) 
crest_K = sqrt(srest_K$d) * t(srest_K$u) 

# Create a grid of variance component proportions for the 4 variance components
h2s=data.frame(h21=c(rep(0.1,10),rep(0.5,10),rep(0.05,10)),h22=c(rep(0.5,10),rep(0.1,10),rep(0.05,10)),h23=c(rep(0.1,10),rep(0.1,10),rep(0.3,10))) %>%
  mutate(h24=1-h21-h22-h23)

# create 4 vectors of random variables with each kinship matrix (or uncorrelated residuals)
# column i of Y is a trait with variance component proportions h2s[i,]
Y=matrix(nrow=nrow(cmoa_K),ncol=dim(h2s)[1],NA)
for(i in 1:dim(h2s)[1]){
  U = matrix(0,nrow = nrow(cmoa_K),ncol = 4)
  U[,1] = crossprod(cmoa_K,rnorm(nrow(cmoa_K)))
  U[,2] = crossprod(cback_K,rnorm(nrow(cback_K)))
  U[,3] = crossprod(crest_K,rnorm(nrow(rest_K)))
  U[,4] = rnorm(nrow(cmoa_K))
  Y[,i] = U %*% t(sqrt(h2s[i,]))
}
apply(Y,2,var)

write.table(Y,file="phenos_round8.txt")


make_K<-function(rootname){
  main<-readGRM(rootname)
  main_K=matrix(0,max(main$grm$id1),max(main$grm$id2))
  main_K[upper.tri(main_K,diag=T)]=main$grm$grm
  main_K[lower.tri(main_K)]=t(main_K)[lower.tri(main_K)]
  return(main_K)
}

normalize_kinmat <- function(kinmat){
  #normalize kinship so that Kij \in [0,1]
  tmp=kinmat/mean(diag(kinmat))
  return(tmp)
}


readGRM <- function(rootname)
{
  #rootname="MOA.SNPs.bed_CentIBS_K"
  bin.file.name <- paste(rootname, ".grm.bin", sep="")
  n.file.name <- paste(rootname, ".grm.N.bin", sep="")
  id.file.name <- paste(rootname, ".grm.id", sep="")
  
  cat("Reading IDs\n")
  id <- read.table(id.file.name, colClasses="character")
  n <- dim(id)[1]
  cat("Reading GRM\n")
  bin.file <- file(bin.file.name, "rb")
  grm <- readBin(bin.file, n=n*(n+1)/2, what=numeric(0), size=4)
  close(bin.file)
  cat("Reading N\n")
  n.file <- file(n.file.name, "rb")
  N <- readBin(n.file, n=n*(n+1)/2, what=numeric(0), size=4)
  close(n.file)
  
  cat("Creating data frame\n")
  l <- list()
  for(i in 1:n)
  {
    l[[i]] <- 1:i
  }
  col1 <- rep(1:n, 1:n)
  col2 <- unlist(l)
  grm <- data.frame(id1=col1, id2=col2, N=N, grm=grm)	
  
  ret <- list()
  ret$grm <- grm
  ret$id <- id
  return(ret)
}
