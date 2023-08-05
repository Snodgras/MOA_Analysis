#!/bin/bash

outname=$1 #(Full_Phenotype_Results/LDAK_EG80WWall_full)
Klist=$2 #(kinship_list.txt)
pheno=$3 #(path to phenotype file)
module load ldak
ldak --reml $outname --mgrm $Klist --mpheno -1 --pheno $pheno --kinship-details NO --dentist YES --constrain YES
