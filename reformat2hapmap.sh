#!/bin/bash

#convert Summary_genotypes_AF.bed into hapmap format
#rs#    ["snp."$1"."$3]
#alleles    [$4"/"$5]
#chrom     $1
#pos     $3
#strand    NA
#assembly#   NA 
#center    NA
#protLSID    NA
#assayLSID    NA
#panelLSID    NA
#QCcode    NA
#B73_GT	[$4$4]
#[NAM]_GT	[if col_namGT == "n.a." ; then NA ; else if col_namGT == "1/1" ; then col5col5; else if col_namGT == "0/0" ; then col4col4 ; fi ; fi; fi ]

python 

with open('Alignment_SNPs.txt') as file:
    for line in file:
        entries = line.split('\t')
        num = len(entries)
        for i in range(5, num):
            if entries[i] == '1/1':
                entries[i] = entries[4]+entries[4]
            if entries[i] == '1/0':
                entries[i] = entries[4]+entries[3]
            if entries[i] == '0/0':
                entries[i] = entries[3]+entries[3]
            if entries[i] == '0/1':
                entries[i] = entries[3]+entries[4]
            if entries[i] == '1/n.a.':
                entries[i] = entries[4]+'N'
            if entries[i] == '0/n.a.':
                entries[i] = entries[3]+'N'
            if entries[i] == 'n.a.':
                entries[i] = 'NN'
        print("snp."+entries[0]+"."+entries[2]+'\t'+entries[3]+"/"+entries[4]+'\t'+entries[0]+'\t'+entries[2]+'\t'+"NA"+'\t'+"NA"+'\t'+"NA"+'\t'+"NA"+'\t'+"NA"+'\t'+"NA"+'\t'+"NA"+'\t'+entries[3]+entries[3]+'\t'+'\t'.join(entries[5:num]))
