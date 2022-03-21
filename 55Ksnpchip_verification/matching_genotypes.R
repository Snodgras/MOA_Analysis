#!/usr/bin/env Rscript
#args <- commandArgs(T)
#truth <- read.table(args[1], sep="\t", stringsAsFactors=FALSE,header=TRUE) #55Kv5 hmp
#query<-read.table(args[2], sep="\t", stringsAsFactors=FALSE,header=TRUE) #alignment or gatk hmp
#tag <- tools::file_path_sans_ext(args[3]) #output file name

truth <- read.table("SNP55K_maizeNAM_AGP5.hmp.txt", sep="\t", stringsAsFactors=FALSE,header=TRUE) #55Kv5 hmp
query<-read.table("filtered_gatkSNPs.tsv", sep="\t", stringsAsFactors=FALSE,header=TRUE) #alignment or gatk hmp
#align<-read.table("filtered_alignSNPs.tsv", sep="\t", stringsAsFactors=FALSE,header=TRUE) #alignment or gatk hmp

#Libraries
library(tidyverse)

#rename the columns for no issues with matching
colnames(truth)<-colnames(truth) %>% str_to_lower()
colnames(query)<-colnames(query) %>% str_to_lower()

#check to make sure chromosome columns are "chr[0-9]" instead of "[0-9]"
print("Checking chromosome column formatting...")
if(!any(str_detect(truth$chrom, "chr"))){ #if not all the chrom entries have "chr" 
  truth$chrom <- str_c("chr",truth$chrom, sep="")
  print("...changed truth chromosome column to have 'chr'...")
}
if(!any(str_detect(query$chrom, "chr"))){ #if not all the chrom entries have "chr" 
  query$chrom <- str_c("chr",query$chrom, sep="")
  print("...changed query chromosome column to have 'chr'...")
}

#remove unnecessary columns and put in same order
print("Selecting columns and reordering...")
truth<-select(truth, c("rs","alleles","chrom","pos","strand","b73","b97","cml103","cml228","cml247",
                       "cml277","cml322","cml333","cml52","cml69","hp301","il14h","ki11",
                       "ki3","ky21","m162w","m37w","mo18w","ms71","nc350","nc358","oh43",
                       "oh7b","p39","tx303"))
query<-select(query, c("rs","alleles","chrom","pos","strand","b73","b97","cml103","cml228","cml247",
                       "cml277","cml322","cml333","cml52","cml69","hp301","il14h","ki11",
                       "ki3","ky21","m162w","m37w","mo18w","ms71","nc350","nc358","oh43",
                       "oh7b","p39","tx303"))

#check to make sure that they are in the same order and named the same
if(all(colnames(truth) == colnames(query))){
  print("The column names are the same in both input files")
} else{ print("ERROR!!! COLUMN NAMES ARE NOT THE SAME IN BOTH INPUT FILES!!!")}

#turn the "NN" to NA####
print("Changing NN and empty to NA...")
truth<-truth %>% na_if("")
query<-query %>% na_if("")
truth<-truth %>% na_if("NN")
query<-query %>% na_if("NN")

#for testing
#truth.10<-filter(truth, chrom == "chr10")
#gatk.10<-filter(gatk, chrom == "chr10")
#align.10<-filter(align, chrom == "chr10")

#right_join(x=truth.10, y=gatk.10, by=c("chrom","pos"), 
#           suffix = c(".truth",".gatk")) %>% 
#  right_join(x = . , y = align.10, by = c("chrom", "pos")) %>% 
#  select(c(contains("alleles"),contains("rs"),"chrom","pos"))

matching_genotypes<-function(col.T, col.Q,strand.T,strand.Q){ #1
  if(!any(is.na(c(col.T,col.Q)))){ #checks for NA's 2
  if(length(c(col.T,col.Q)) == 2){ #checks to make sure we're just looking at one position 3
    #check if the strands are the same 4
    if(strand.T == "+"){ #5
      if(col.T == col.Q){ #if the strings are exactly the same 6
        m.value <- 2 #7
        } else{ #see if the strings are heterozygous and half match 8
          if(any(str_detect(col.Q, pattern = c(str_sub(col.T, 1,1),str_sub(col.T, 2,2))))){ #9
          m.value <- 1 #10
          } else{ #otherwise they don't match at all #11
            m.value <- 0 #12
          } #13
        } #14
    } else{ #if strand.T is negative, change it to positive genotype 15
        if(strand.T == "-"){ #16
          if(col.T == "AA"){col.T<-"TT"}else{ #17
            if(col.T == "TT"){col.T<-"AA"}else{ #18
              if(col.T == "GG"){col.T<-"CC"}else{ #19
                if(col.T == "CC"){col.T<-"GG"} #20
              } #21
            } #22
          } #23 #then do matching #24
    if(col.T == col.Q){ #if the strings are exactly the same 25
      m.value <- 2 #26
    } else{ #see if the strings are heterozygous and half match 27
      if(any(str_detect(col.Q, pattern = c(str_sub(col.T, 1,1),str_sub(col.T, 2,2))))){ #28
        m.value <- 1 #29
      } else{ #otherwise they don't match at all #30
        m.value <- 0 #31
      } #32
    }
    }#33
    } #34
  }
    else{#Length of matching elements is different than 1 35
    m.value <- NA} #36
  } #37
    else{ #39
    m.value<-NA #40
  } #41
  return(m.value) #42
  }
  
count_na<-function(x){sum(is.na(x))} 
#Matching loop
print("Starting matching loop...")
matching.output<-tibble(rs = NA,
                        alleles = NA,
                        chrom = NA,
                        pos = NA,
                        b73 = NA,
                        b97= NA,
                        cml103= NA,
                        cml228= NA,
                        cml247= NA,
                        cml277= NA,
                        cml322= NA,
                        cml333= NA,
                        cml52= NA,
                        cml69= NA,
                        hp301= NA,
                        il14h= NA,
                        ki11= NA,
                        ki3= NA,
                        ky21= NA,
                        m162w= NA,
                        m37w= NA,
                        mo18w=NA,
                        ms71= NA,
                        nc350= NA,
                        nc358= NA,
                        oh43= NA,
                        oh7b= NA,
                        p39= NA,
                        tx303= NA)
for(i in 1:nrow(truth)){
  matcher<- filter(query, chrom == truth$chrom[i] & pos == truth$pos[i])
  m.value <- c()
  for(j in 6:ncol(truth)){
   m.value<-c(m.value,matching_genotypes(truth[i,j],matcher[,j],
                                         truth[i,"strand"],matcher[,"strand"])) 
  }
  matching.output<-add_row(matching.output, #will have to edit for actual columns
                           rs = truth$rs[i],
                           alleles = truth$alleles[i],
                           chrom = truth$chrom[i],
                           pos = truth$pos[i],
                           b73 = m.value[1],
                           b97= m.value[2],
                           cml103= m.value[3],
                           cml228= m.value[4],
                           cml247= m.value[5],
                           cml277= m.value[6],
                           cml322= m.value[7],
                           cml333= m.value[8],
                           cml52= m.value[9],
                           cml69= m.value[10],
                           hp301= m.value[11],
                           il14h= m.value[12],
                           ki11= m.value[13],
                           ki3= m.value[14],
                           ky21= m.value[15],
                           m162w= m.value[16],
                           m37w= m.value[17],
                           mo18w= m.value[18],
                           ms71= m.value[19],
                           nc350= m.value[20],
                           nc358= m.value[21],
                           oh43= m.value[22],
                           oh7b= m.value[23],
                           p39= m.value[24],
                           tx303= m.value[25])
}
matching.output<-matching.output[-1,] #gets rid of initializing NA's
print("...finished matching loop...Calculating matching stats...")
matching.output<- mutate(matching.output, 
                         Total_match = rowSums(matching.output[,5:ncol(matching.output)], na.rm = TRUE),
                         Missing = apply(matching.output[,5:ncol(matching.output)], 1, count_na),
                         Prop_match = Total_match/((25-Missing)*2)) #edit 25 to be number of founder columns
write_csv(matching.output,path=paste0(tag))
paste("Output written to ",tag)