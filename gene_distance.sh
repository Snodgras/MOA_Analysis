#!/bin/bash
GFF=$1
MOABED=$2

module load bedops
convert2bed --input=gff --output=bed < $GFF > gff.bed

#need to separate out only the gene coordinates (not exon, chromosome, mRNA, etc.)
awk -v OFS='\t' '$8 ~ /gene/ {print $0}' gff.bed | sort -k1,1 -k2,2n > sorted.gene.bed

#sort the peak bedfile
sort -k1,1 -k2,2n $MOABED > sorted.$MOABED

module load bedtools2
bedtools closest -D ref -t first -a sorted.$MOABED -b sorted.gene.bed > out.distance.${2}.${1%_*}.bed

#Slurm command
bash gene_disance ZmB73v5.gff EG80.WW.all_peaks.merged.inclCML333.bed				      
