#Analyzing the Matching Output
library(tidyverse)
library(ggplot2)

gatk.results<-read_csv("gatk_matchingOutput.csv")
align.results<-read_csv("align_matchingOutput.csv")

#filter out rows where missing == 25
gatk.results<-filter(gatk.results, Missing != 25)
align.results<-filter(align.results, Missing != 25)

#distribution/summary statistics of missing-ness
select(gatk.results, Missing) %>% summary()
select(align.results, Missing) %>% summary()

ggplot(gatk.results, aes(x = Missing))+
  geom_histogram(binwidth = 1)+
  ggtitle("GATK Missing")+ theme_minimal()
ggsave("gatk.strandcorrected.Missing.png", device='png')
ggplot(align.results, aes(x = Missing))+
  geom_histogram(binwidth = 1, fill = "purple")+
  ggtitle("Alignment Missing")+theme_minimal()
ggsave("align.strandcorrected.Missing.png", device='png')

#distribution/summary statistics of matching
select(gatk.results, Prop_match) %>% summary()
select(align.results, Prop_match) %>% summary()

ggplot(gatk.results, aes(x = Prop_match))+
  geom_histogram(binwidth = 0.1)+
  ggtitle("GATK Proportion Matching 55SNP")+theme_minimal()
ggsave("gatk.strandcorrected.PropMatching.png", device='png')
ggplot(align.results, aes(x = Prop_match))+
  geom_histogram(binwidth = 0.1, fill="purple")+
  ggtitle("Alignment Proportion Matching 55SNP")+theme_minimal()
ggsave("align.strandcorrected.PropMatching.png", device='png')

#let's see how many overlap between the two and compare
align.results$rs
gatk.results$rs

joined.results<-full_join(x=select(gatk.results, c(rs, Total_match, Missing, Prop_match)),
          y=select(align.results, c(rs, Total_match, Missing, Prop_match)),
          by="rs") %>% na.omit()
ggplot(joined.results)+
  geom_histogram(aes(x=Missing.x), fill = "black", alpha = 0.5, binwidth = 1)+
  geom_histogram(aes(x=Missing.y), fill = "purple", alpha = 0.5, binwidth = 1)+
  ggtitle("Missing with SNPs shared bet. GATK and Align")+
  labs(x="Missing")+theme_minimal()
ggsave("shared.strandcorrected.Missing.png", device='png')

ggplot(joined.results)+
  geom_histogram(aes(x=Prop_match.x), fill = "black", alpha = 0.5, binwidth = 0.1)+
  geom_histogram(aes(x=Prop_match.y), fill = "purple", alpha = 0.5, binwidth = 0.1)+
  ggtitle("Prop Matched with SNPs shared bet. GATK and Align")+
  labs(x="Proportion Matched")+theme_minimal()
ggsave("shared.strandcorrected.PropMatching.png", device='png')

#compared to testing strand matching on chrom 10
gatk.results.10<-filter(gatk.results, chrom == "chr10")
matchedonstrand.10<-read_csv("testingStrandMatchingonChr10.csv")
matchedonstrand.10<-filter(matchedonstrand.10, Missing != 25)
ggplot()+
  geom_histogram(data = gatk.results.10, aes(x=Missing), 
                 fill = "black", binwidth = 1, alpha = 0.5)+
  geom_histogram(data = matchedonstrand.10, aes(x=Missing), 
                 fill = "red", binwidth = 1, alpha = 0.5)
ggplot()+
  geom_histogram(data = gatk.results.10, aes(x=Prop_match), 
                 fill = "black", binwidth = 0.1, alpha = 0.5)+
  geom_histogram(data = matchedonstrand.10, aes(x=Prop_match), 
                 fill = "red", binwidth = 0.1, alpha = 0.5)+
  theme_minimal()+ggtitle("Strand not matched vs matched for gatk, chr10")
ggsave("gatk.proportionmatched.strandvsnostrand.chr10.png",device = "png")

