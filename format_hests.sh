#!/bin/bash

d=$1 #directory where the output files are
outname=$2
phenotypefile=$3

cd $d/
for i in *.reml ; do
	n=$(echo $i | cut -f 2 -d \. )
	let s=$n+2
	pheno=$(head -n 1 $phenotypefile | cut -f $s)
	tail -n 5 $i | tr " " "\t" | datamash transpose | grep Her > 0.temp 
	awk -v OFS='\t' 'NR == 1 || NR == 2 {print $0}' 0.temp > 1.temp
	awk -v OFS='\t' 'NR == 1 || NR == 3 {print $0}' 0.temp > 2.temp
	paste 1.temp 2.temp | tr " " "\t" | awk -v OFS='\t' -v p=$pheno '{print $0,p}' - > her.$n.out
	echo "done with $i"
done

head -n 1 her.1.out > ../$outname.h_estimates.txt 
for i in her*.out ; do grep -v Component $i >> ../$outname.h_estimates.txt ; done

rm *.temp her.*.out

head -n 1 *.1.share > ../$outname.share.txt
for i in *.share ; do
	n=$(echo $i | cut -f 2 -d \. )
	let s=$n+2
	pheno=$(head -n 1 $phenotypefile | cut -f $s)
	cat $i | tr " " "\t" | awk -v OFS='\t' -v p=$pheno 'NR != 1 {print $0,p}' - >> ../$outname.share.txt
	echo "done with $i"
done

cd ..
pwd

#rm -r $d/