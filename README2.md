#Notes about MOA 

##Data from:
ftp://ftp.mpipz.mpg.de/

##Creating reference from 2 genome fastas
* need to change fasta id names so there's no conflicts
`bioawk -c fastx '{print ">Tzi8-"$name"\n"$seq}' <file name> | fold | less`
script: `namefastaids.sh`
```
sed -i -e 's/\$Ki11/Ki11/g' ref_Ki11.fasta
sed -i -e 's/MS71/Ms71/g' ref_Ms71.fasta
sed -i -e 's/\$Ki11/Ki11/g' ref_B73Ki11.fasta
sed -i -e 's/MS71/Ms71/g' ref_B73Ms71.fasta
```
* then cat the genomes together
`cat ref_B73.fasta ref_Oh43.fasta >> ref_B73Oh43.fasta`


for i in 00_reference_genomes/*; do 
if [[ ( ${i} != *B73* && ${i} != *Oh43* && ${i} != *CML69* ) && ( ${i} == *.fasta ) ]]
then	
	nam=${i#00*ref_}
	nam=${nam%.fasta}
	#cat ref_B73.fasta ${i} >> ref_B73${nam}.fasta
	echo ref_B73${nam}.fasta
fi ;
done


* then create an index file
script: `indexref.sh`
_DONE_

## Merge (NGMerge) paired-end reads to single-end reads
* Try merging the reads first (NGmerge) _DONE_
`merge_reads.sh` then `trim_merged_reads.sh`
	* QC the merged reads _DONE_
	* remove longer than 100bp reads _DONE_
	`aln_bwa_d5_NGmerge_100trimmed.fastq.sam	`
	`_sorted.bam`

## FastQC the raw reads
```
for i in A619 B97 CML277 CML322 CML333 CML69 IL14H Ki11 Ki3 Ky21 M162W Mo18W Ms71 NC358 Oh43 Tx303 Tzi8 ; do
>         cd 00_B73x${i}_data ;
>         for j in *.gz ; do
>                 echo fastqc -t 2 ${j} \> /work/LAS/mhufford-lab/snodgras/MOA/00_fastqc/fastqc_${j} >> /work/LAS/mhufford-lab/snodgras/MOA/commands.txt ;
>         done ;
>         cd .. ;
>         echo Done FastQC with ${i}
> done
python makeSLURMs.py 1 commands.txt 
#the script is modified to include the following lines:
w.write("module --ignore-cache load fastqc\n")
        w.write("module --ignore-cache load picard\n")
        w.write("PICARD_HOME=$(dirname $(which picard))\n")
        w.write("java -Xmx150g -Djava.io.tmpdir=$TMPDIR -jar ${PICARD_HOME}/picard.jar\n")

```
Trim Galore to remove extra adapter sequence
```
module load trimgalore/0.6.4-py3-36bynaj
trim_galore --paired *.R1.fq.gz *R2.fq.gz
```
FastQC of trimmed file
```
module --ignore-cache load fastqc
module --ignore-cache load picard
PICARD_HOME=$(dirname $(which picard))
java -Xmx150g -Djava.io.tmpdir=$TMPDIR -jar ${PICARD_HOME}/picard.jar
fastqc -t 2 A619d21_1_L4_val_1.fq.gz > fastqc_A619d21_1_trimmed.fq.gz

###Multiple files to QC
ml py-multiqc/1.5-py2-lqqx3ht
multiqc . #the . is for the directory of files to use
```

## Remove Adapters, merge, filter, map (merged read mapping)

Use `seqpurge_before_merge.sh` to remove adapter sequence
```
#different names for CML322 Ki11 Ki3 and Mo18W
for i in B97 CML277 CML333 IL14H Ky21 M162W Ms71 NC358 Tx303 Tzi8 ; do 
	for j in w d ; do 
		for k in {1..20} ; do 
			for l in 2 4 ; do 
				if [[ -f 00_B73x${i}_data/${i}${j}${k}_L${l}_1.fq.gz ]] ; then 
				echo bash scripts/seqpurge_before_merge.sh 00_B73x${i}_data ${i} ${j}${k} 01a_seqpurge_before_merge ${l} >> commands.txt ; 
				echo bash scripts/seqpurge_for_paired_end.sh 00_B73x${i}_data ${i} ${j}${k} 01b_seqpurge_for_paired_end ${l} >> commands.txt ; 
				fi ;
 			done ;
  		done ; 
  	done ; 
done

python makeSLURMs.py 1 commands.txt

for i in {20..29}; do sbatch --dependency=afterok:1028663 commands_${i}.sub ; done
for i in {30..39}; do sbatch --dependency=afterok:1028676 commands_${i}.sub ; done
for i in {40..49}; do sbatch --dependency=afterok:1028687 commands_${i}.sub ; done
for i in {50..59}; do sbatch --dependency=afterok:1028697 commands_${i}.sub ; done
for i in {60..69}; do sbatch --dependency=afterok:1028707 commands_${i}.sub ; done
for i in {70..79}; do sbatch --dependency=afterok:1028727 commands_${i}.sub ; done
for i in {80..89}; do sbatch --dependency=afterok:1028737 commands_${i}.sub ; done
for i in {90..99}; do sbatch --dependency=afterok:1028747 commands_${i}.sub ; done
for i in {100..119}; do sbatch --dependency=afterok:1033872	commands_${i}.sub ; done
```
jobs 1027149..1027267
Use `merge_reads.sh` to merge reads in `01a_seqpurge_before_merge`

```
for i in A619 ; do
        for j in {21..23} ; do
                for k in d w ; do
                        if [[ -f 01a_seqpurge_before_merge/${i}${k}${j}_1.trim.gz ]]  ; then
                        cd 01a_seqpurge_before_merge
                        bash /work/LAS/mhufford-lab/snodgras/MOA/scripts/merge_reads.sh ${i}${k}${j}_1.trim.gz ${i}${k}${j}_2.trim.gz ;
                        cd ..
                        fi ;
                done ;
        done ;
        mv 01a_seqpurge_before_merge/*NGmerge* .
done
```
To make it run in parallel
```
for i in CML322 CML69 Ki11 Ki3 Mo18W Oh43 B97 CML277 CML333 IL14H Ky21 M162W Ms71 NC358 Tx303 Tzi8 ; do
	for j in w d ; do
		for k in {1..20} ; do
			if [[ -f 01a_seqpurge_before_merge/${i}${j}${k}_1.trim.gz ]] ; then
				echo cd 01a_seqpurge_before_merge \; bash /work/LAS/mhufford-lab/snodgras/MOA/scripts/merge_reads.sh ${i}${j}${k}_1.trim.gz ${i}${j}${k}_2.trim.gz \; cd .. >> merge_commands.txt ;
			fi ; 
		done
	done
done

python makeSLURMs.py 1 merge_commands.txt

for i in {0..9}; do sbatch merge_commands_${i}.sub ; done [1037128:1037137]
for i in {10..19}; do sbatch --dependency=afterok:1037137 merge_commands_${i}.sub ; done [1037148:1037157]
for i in {20..29}; do sbatch --dependency=afterok:1037157 merge_commands_${i}.sub ; done [1037158:1037167]
for i in {30..39}; do sbatch --dependency=afterok:1037167 merge_commands_${i}.sub ; done [1037168:1037177]
for i in {40..49}; do sbatch merge_commands_${i}.sub ; done [1038822:1038831]
for i in {50..59}; do sbatch --dependency=afterok:1038831 merge_commands_${i}.sub ; done [1038834:1038843]
for i in {60..69}; do sbatch --dependency=afterok:1038843 merge_commands_${i}.sub ; done [1038844:1038853]
for i in {70..79}; do sbatch --dependency=afterok:1038853 merge_commands_${i}.sub ; done [1038854:1038863]
for i in {80..89}; do sbatch --dependency=afterok:1038863 merge_commands_${i}.sub ; done [1038864:1038873]
for i in {90..99}; do sbatch --dependency=afterok:1038873 merge_commands_${i}.sub ; done [1038874:1038883]

```
###Filtering out the 101 or longer merged reads
Testing on the A619 reads using `02.5a_filtered_merged_reads.sh`
```
mkdir 02.5a_filtered_merged_reads/
for i in 02_merged_reads/A619* ; do bash scripts/02.5a_filtered_merged_reads.sh $i ; done
```
The fastqc report looked good for length (though A619 d23 looked different from the rest)
Making the commands to allow it to run in parallel
```
for i in CML322 CML69 Ki11 Ki3 Mo18W Oh43 B97 CML277 CML333 IL14H Ky21 M162W Ms71 NC358 Tx303 Tzi8 ; do
	if ls 02_merged_reads/${i}*.gz &> /dev/null ; then 
		echo for j in 02_merged_reads/${i}* \; do bash scripts/02.5a_filtered_merged_reads.sh \$j \; done >> filter_commands.txt
	fi
done

python makeSLURMs.py 1 filter_commands.txt
for i in filter_commands*.sub ; do sbatch $i ;done [1055052:1055067]
```
### running the mapping on the merged  reads
```
for i in B97 CML277 CML322 CML333 CML69 IL14H Ki11 Ki3 Ky21 M162W Mo18W Ms71 NC358 Oh43 Tx303 Tzi8 ; do
	for j in {1..20} ; do
		if [ -f 02.5a_filtered_merged_reads/${i}d${j}_1.trim_${i}d${j}_2.trim_NGmerge.gz_100trimmed.fastq ] ; then
			echo bash bwa-map-and-process-diploid-3a.sh ref_B73${i}.fasta ${i}d${j}_1.trim_${i}d${j}_2.trim_NGmerge.gz_100trimmed.fastq >> commands.txt ; 
			echo bash bwa-map-and-process-diploid-3a.sh ref_B73${i}.fasta ${i}w${j}_1.trim_${i}w${j}_2.trim_NGmerge.gz_100trimmed.fastq >> commands.txt ;
			echo bash bwa-map-and-process-haploid-3b.sh ref_B73.fasta ref_${i}.fasta ${i}d${j}_1.trim_${i}d${j}_2.trim_NGmerge.gz_100trimmed.fastq >> commands.txt ; 
			echo bash bwa-map-and-process-haploid-3b.sh ref_B73.fasta ref_${i}.fasta ${i}w${j}_1.trim_${i}w${j}_2.trim_NGmerge.gz_100trimmed.fastq >> commands.txt ;
		fi ;
	done ; 
done

for i in 02.5a_filtered_merged_reads/* ; do if [[ ( ${i} != *A619* ) ]] ; then mv $i . ; fi ; done
for i in 00_reference_genomes/*.fasta ; do if [[ ( ${i} != *A619* ) ]] ; then mv $i . ; fi ; done

python makeSLURMs.py 1 commands.txt
```
Ran out of disk space
Moving the command scripts/files to `ptmp/LAS/snodgras`
```
for i in /work/LAS/mhufford-lab/snodgras/MOA/ref*.fasta;
	do ln -s ${i} ${i#/work/LAS/mhufford-lab/snodgras/MOA/} ;
	echo done linking ${i};
done
for i in /work/LAS/mhufford-lab/snodgras/MOA/*fastq ; 
	do ln -s ${i} ${i#/work/LAS/mhufford-lab/snodgras/MOA/} ;
	echo done linking ${i};
done
mv /work/LAS/mhufford-lab/snodgras/MOA/*.sh .
cp /work/LAS/mhufford-lab/snodgras/MOA/*.py . 
```
Some jobs timed out and left a bunch of temp files. Here's how to remove them:
```
ls *temp* | cut -f 1 -d . | uniq > removal_list.txt
while read -r line ; do rm ${line}*.bam ; rm ${line}*.txt ; done < removal_list.txt
```
Need to figure out which files "completed" but didn't write files. 
Should have `File write failed (wrong size)` in the .e files
```
grep -riIl "File write failed (wrong size)" commands*.e*
```
Testing Step 4 filtering with the NC358 reads
```
mv *ids.txt id-files
for j in B97 CML277 CML322 CML333 CML69 IL14H Ki11 Ki3 Ky21 M162W Mo18W Ms71 NC358 Oh43 Tx303 Tzi8 ; do 
	for f in $j*-diploid-mapped.bam ; 
		do g=$(basename $f | cut -f 1 -d "-") ; 
		cat id-files/*${g}* | sort --parallel 10 | uniq >> id-files/${g}-filtering.ids; 
		echo "../filter-bam-files-4.sh $f ${g}-filtering.ids"; 
	done >> filter.cmds ; 
done
```
There were some files that could not pass the filtering criteria and threw an error. 
These were B97d4, NC358w4, NC358d2, CML69w3, CML69w2. 
To get around this, I altered the `filter-bam-files-4.sh` to have lenient validation. 
This throws a warning error instead of halting the job altogether. 
```
picard FilterSamReads I=${bam} O=${bam%.*}-filtered.bam READ_LIST_FILE=${ids} FILTER=excludeReadList VALIDATION_STRINGENCY=LENIENT
```
Should not affect the way it handles these jobs compared to all the other ones that worked. 

## Remove Adapters, filter, map (paired read mapping)
Use `seqpurge_for_paired_end.sh`

## Checking the mapping stats for Tzi8
`picard_mappingstats.sh` run on `/ptmp/LAS/mhufford/NAM-MOA/Arun_notebook/scripts/id-files/Tzi8*filtered-RG.bam`

## Putting files on the gdrive for Thomas
```
rclone --drive-shared-with-me copy file gdrive:masked-NAM-genomes

rclone --drive-shared-with-me copy CML277w4_L2_2.fq.gz gdrive:masked-NAM-genomes/raw_reads/

for i in *NC358*.bam ; do rclone --drive-shared-with-me copy $i gdrive:masked-NAM-genomes/NC358-intermediate-files/ ; done
for i in *filtered.reads ; do rclone --drive-shared-with-me --progress copy $i gdrive:masked-NAM-genomes/filtered-reads-2.0/ ; done
```

# Using MOA-Seq to identify allele-specific TF binding sites (ASBs)

## For SNPs
Use your B73xNAM VCF; filter VCF (bcftools view) to include -i 'GT="AA"' (hom) -M 2 (biallelic) -type snsps (SNPs only) to produce ``$homo_M2_SNPsonly_vcf`
* In addition to filter reads we also need to remove variants which overlap with CNVs. E.g. A B73 region has a SNP to a region duplicated in NAM. The reads covering the duplication in NAM are removed, giving the B73 region a mapping bias.

## Filter $homo_M2_SNPsonly_vcf for SNPs which overlap to positions in $B73_multi_map_read_names_txt/$NAM_multi_map_read_names_txt to produce $homo_M2_SNPsonly_CNV_vcf

## For each SNP in $homo_M2_SNPsonly_CNV_vcf create second SNP entry with chromosome_NAM position instead of chromosome_B73 position; for each chromosome_NAM SNP reverse ref and alt allele (awk) to produce $diploid_SNPs_vcf 
This step is needed since VCF was for the haploid not diploid genome. 
To make ASEReadcounter or pileup work every SNP needs a B73 and NAM position.

## Change GT value (awk) for all SNPs in $diploid_SNPs_vcf from 1/1 to 0/1 (hom to het) for compatibilty with ASEReadCounter to produce $diploid_SNPs_hom_to_het_vcf

## Use $MAPQ13_reads_filtered_bam and $diploid_SNPs_hom_to_het_vcf as inputs for ASEReadCounter with --min-base-quality 20

## Calculate allele frequency ($AF) from ASEReadCounter: $AF = Ref_count chromosome_B73 SNP / (ref_count chromosome_B73 SNP + ref_count chromosome_NAM SNP) 
A $AF of 1 means 100% reads B73, a $AF of 0 means 100% reads NAM parent. 
Expected $AF in F1 = 0.5.

## Binomal test for $AF against expected 0.5 null hypothesis, adjust for multiple-testing.  
​
## For INDELS

Indels are first split into two groups: 
1) smaller than max MOA read length of 110 bp ($Small_INDELs) or 2) larger than 110 bp ($Large_INDELs)
For `$Small_INDELs` the genome position will be extended to 110 bp (incl. the INS) for `chromosome_B73` and `chromosome_NAM` position  
All INDELs are treated as insertions (DEL means B73 has insertion, INS means NAM has insertion). E.g. B73 has a 20 bp INS; Start and Stop of INS will be extended by 45 bp for B73 (20+45+45=110bp), 
however, for the NAM parent the region is only 90 bp (45+45).
The average per bp coverage between the `B73_region` (e.g. 110 bp) and `NAM_region` (e.g. 90bp) are compared.
​
For $Large_INDELs if any bp position in the INDEL region is above the cup-off (e.g. 150 reads) to be a MOA peak, it will be included as an ASB  
Since the INDEL is larger than MOA reads, any peak in this INDEL would be considered an allele-specific peak if it is above the peak threshold. 


# MOA Results
Key for the files
```
The Summary files contain values for all SNPs the were found in one or more of the 16 hybrids. SNPs were determined from the cactus hal file.
Columns in the Summary files:
Only once for all lines:
CHR: B73 chromosome
POS: B73 position
B73al: Base in B73 and NAM if 0/0
NAMal: Base in NAM if 1/1
For each NAM line:
ID: NAMchrosome_NAMposition
Genotype: 0/0 if B73, 1/1 if NAM has a SNP
B73_CPM: CPM in B73
NAM_CPM: CPM in NAM
B73_peak: 0 if below 0.3, 1 if above or equal 0.3
NAM_peak: 0 if below 0.3, 1 if above or equal 0.3
PF: "post frequency": B73_CPM/(B73_CPM+NAM_CPM), if both are 0, PF is denoted as n.r. (no reads), if both are below the peak threshold, the PF is calculated but n.p. is written in front of it.
Abbreviations:
n.r. no reads - no reads in either B73 or NAM
n.p. no peak - both allele CPM are below threshold for peak (0.3)
n.a. no liftover position at this SNP
```

# VCAP analysis
Using Tassel 5 pipeline`https://docs.google.com/document/d/1DvwgkFllI-WjEMrr2J69wCNTSMY7Zya42FrYE_28FNM/edit?pli=1#heading=h.7vm9br7kry8`
Practice using the files on Cyverse at `/iplant/home/shared/panzea/VCAP/`
```
ml singularity
singularity shell --bind $PWD /home/snodgras/irods.simg
iget -r /iplant/home/shared/panzea/VCAP
```
Creating the kinship matrices
```
module load tassel
run_pipeline.pl \
-Xmx14g \
-importGuess VCAP/genotypes/NAM/NAM_HM32_UnimpMAF01Cov05.hmp.txt.gz.lix \
-FilterSiteBuilderPlugin \
-bedFile VCAP/bedFiles/Zea_mays_AGPv3_CDS_chr10.bed \
-endPlugin \
-KinshipPlugin \
-method Centered_IBS \
-endPlugin \
-export NAM_HM32_UnimpMAF01Cov05_CentIBS_chr10kinshiptest.txt.gz \
-exportType SqrMatrix

run_pipeline.pl \
-Xmx14g \
-importGuess VCAP/genotypes/NAM/NAM_HM32_UnimpMAF01Cov05.hmp.txt.gz.lix \
-KinshipPlugin \
-method Centered_IBS \
-endPlugin \
-export NAM_HM32_UnimpMAF01Cov05_CentIBS_wholegenomekinshiptest.txt.gz \
-exportType SqrMatrix
```
We need to have kinship matrix for the stuff we want to partition.
Kinship remaining =  kinship of whole genome - kinship of partition 1 - kinship of partition 2
Given we only have the bed files for the CDS of each chromosome, we could do:
Whole genome kinship - CDS kinship = kinship of remaining genome parts
```
module load tassel
run_pipeline.pl \
-Xmx14g \
-importGuess VCAP/genotypes/NAM/NAM_HM32_UnimpMAF01Cov05.hmp.txt.gz.lix \
-FilterSiteBuilderPlugin \
-bedFile Zea_mays_AGPv3_CDS_allchr.bed \
-endPlugin \
-KinshipPlugin \
-method Centered_IBS \
-endPlugin \
-export NAM_HM32_UnimpMAF01Cov05_CentIBS_CDSkinshiptest.txt.gz \
-exportType SqrMatrix
```
_Note_: Need this argument to read into LDAK: `-exportType SqrMatrixBin`
Doing the rest of genome kinship subtraction
```
partialkin=$1
wholekin=$2
ml tassel
run_pipeline.pl -Xmx14g \
	-fork1 -importGuess $partialkin \
	-SubtractDistanceMatrixPlugin -wholeMatrix $wholekin \
	-endPlugin \
	-export mdp_genotype_rest \
	-exportType SqrMatrixBin
```
Rolling that all into 1:
```
module load tassel
run_pipeline.pl \
	-Xmx64g \
	-fork1 -importGuess VCAP/genotypes/NAM/NAM_HM32_UnimpMAF01Cov05.hmp.txt.gz.lix \
	-fork2 -FilterSiteBuilderPlugin -bedFile Zea_mays_AGPv3_CDS_allchr.bed -endPlugin \
	-input1 -KinshipPlugin -method Centered_IBS -endPlugin \
	-fork3 -export CDSallchr_Centered_IBS_kinship1 -exportType SqrMatrixBin -input2 \
	-fork4 -SubtractDistanceMatrixPlugin -wholeMatrix NAM_HM32_UnimpMAF01Cov05_CentIBS_wholegenomekinshiptest.txt.gz -endPlugin \
	-input2 -export CDSallchr_Centered_IBS_kinshipRest -exportType SqrMatrixBin

```

Worked!
Make kinship.txt
```
CDSallchr_Centered_IBS_kinship1
CDSallchr_Centered_IBS_kinshipRest
```
Running ldak manually
```
ldak --reml LDAK_CDSallchr --mgrm kinship_list.txt --pheno VCAP/phenotypes/NAM/NAM_DaysToAnthesis_multiblup.txt --kinship-details NO
```
The kinship details flag is if you don't have a `.details` file with your kinship matrices
How big a deal is this warning?
```
Reading kinship matrix with stem CDSallchr_Centered_IBS_kinship1
Warning, Kinship 1 has average diagonal value 1.481743 (usually this is 1)
Reading kinship matrix with stem CDSallchr_Centered_IBS_kinshipRest
Warning, Kinship 2 has average diagonal value 1.429918 (usually this is 1)
```

Let's just try this example code for VCAPScanPlugin from the powerpoint
```
run_pipeline.pl -importGuess *.hmp.txt.gz -VCAPScanPlugin -method Chromosome -phenotypeFile phenotype.txt -endPlugin
run_pipeline.pl -Xmx14g -importGuess VCAP/genotypes/NAM/NAM_HM32_UnimpMAF01Cov05.hmp.txt.gz -VCAPScanPlugin -method Chromosome -phenotypeFile VCAP/phenotypes/NAM/NAM_DaysToAnthesis_multiblup.txt -endPlugin
```

Will not feed into LDAK because the version calls are different (plugin calls for 4.9, HPC has 5.1), but the kinship worked 
So trying to run manually with the VCAPscan plugin files:
```
ldak --reml LDAK_test_list0 --mgrm kinship_list0.txt --pheno VCAP/phenotypes/NAM/NAM_DaysToAnthesis_multiblup.txt --kinship-details NO
```
output is:
```
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
LDAK - Software for obtaining Linkage Disequilibrium Adjusted Kinships and Loads More
Version 5.1 - Help pages at http://www.ldak.org
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

There are 4 pairs of arguments:
--reml LDAK_test_list0
--mgrm kinship_list0.txt
--pheno VCAP/phenotypes/NAM/NAM_DaysToAnthesis_multiblup.txt
--kinship-details NO

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

Performing generalized REML with 2 kinship matrices and 0 regions

Variance components can be negative (to prevent this, use "--constrain YES")

Consider using "--covar" to provide covariates
Consider using "--top-preds" to include (strongly-associated) predictors; these will be treated as fixed effects, except that the variance they explain will count towards total heritability

If memory is an issue, use "--memory-save YES" and kinship matrices will be read on-the-fly

If you would like to perform a permutation analysis, add "--permute YES"

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

There are 5258 samples common to the 2 kinship matrices
Checking responses for 4635 samples from VCAP/phenotypes/NAM/NAM_DaysToAnthesis_multiblup.txt
Due to missing phenotypic values, the number of samples is reduced from 5258 to 4418; to instead pad mising values, use "--dentist YES"
Will be using 4418 samples

Reading responses for 4635 samples from VCAP/phenotypes/NAM/NAM_DaysToAnthesis_multiblup.txt

Reading kinship matrix with stem ./Kinship_10_000000229563_000147291022
Warning, Kinship 1 has average diagonal value 1.432753 (usually this is 1)
Reading kinship matrix with stem ./Kinship_10_000000229563_000147291022Rest
Warning, Kinship 2 has average diagonal value 1.430988 (usually this is 1)

Iter	Her_K1	Her_K2	Her_All	Likelihood	Difference	Target	Num_Constrained
Start	0.2500	0.2500	0.5000	-9588.12	n/a		0.001000	0
1	0.2243	0.3157	0.5400	-9300.48	287.636808	0.001000	0
2	0.1949	0.4058	0.6007	-8979.70	320.779466	0.001000	0
3	0.1593	0.5305	0.6898	-8627.29	352.408799	0.001000	0
4	0.1135	0.7028	0.8163	-8302.09	325.200574	0.001000	0
5	0.0972	0.7581	0.8553	-8270.73	31.364272	0.001000	0
6	0.1011	0.7442	0.8453	-8268.34	2.389273	0.001000	0
7	0.1014	0.7433	0.8447	-8268.33	0.006478	0.001000	0
8	0.1014	0.7433	0.8447	-8268.33	0.000000	0.001000	0

Main results saved in LDAK_test_list0.reml

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
Mission completed. All your basepair are belong to us :)
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
```

####NAM RIL DATA
NAM SNP projections on the RIL population `/iplant/home/shared/NAM/Misc/NAM-SV-projected-V8`
```
ml singularity
singularity shell --bind $PWD /home/snodgras/irods.simg
iget /iplant/home/shared/NAM/Misc/NAM-SV-projected-V8/*SNPs-only*
```
Phenotype data from: Merritt Burch `all_NAM_phenos.csv`

cleaning up datasets to be converted to PLINK format
`http://dougspeed.com/phenotypes-and-covariates/`
```
head -n 1 Peiffer2014_cornelldownload.txt | cut -f 5,12,13 > Peiffer2014_heighttraits.txt
grep NAM Peiffer2014_cornelldownload.txt | cut -f 5,12,13 >> Peiffer2014_heighttraits.txt

head -n 1 Cook2012_SupplementalData1.txt | cut -f 13,15-17 > Cook2012_kernelcomptraits.txt
grep Z0 Cook2012_SupplementalData1.txt | cut -f  13,15-17 >> Cook2012_kernelcomptraits.txt

cut -f 1,4-10  Brown_INFLO7_BLUPs.txt > Brown2011_inflotraits.txt

awk -v FS=, -v OFS='\t' -v b=9 '{for (i=b; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}' all_NAM_phenos.csv > temp.txt
awk -v FS=, -v OFS='\t' '{print $1,$1}' all_NAM_phenos.csv | paste - temp.txt > all_NAM_phenos.plink
rm temp.txt
sed -i "s/\"//g" all_NAM_phenos.plink
```

Geting the MOA Data
```
#should really do this on the dtn...
ml singularity
singularity shell --bind $PWD /home/snodgras/irods.simg
ils /iplant/home/aseetharam/NAM_MOA/MOA_tables
iget -r /iplant/home/aseetharam/NAM_MOA/MOA_tables
```
Read me file
For each condition (WW- well watered and DS- drought stressed) there are two archives, _all and _uni. 
Each archive contains the count files for the respective condition, one file per chromosome. 
Counts were generated from STAR mapped, CPM normalised data. For _uni, only reads mapping uniquely to the diploid genome were kept. 
This means all reads that did not contain an SNP were not counted, which can lead to misleading peak shapes. 
On the other hand, everything that is counted for each allele does really belong there. 
So even for locations with 0/0 genotype, if there are counts, it indicates a SNP in the vicinity made it possible to get allele-specific data for this position. 
Thus the 0/0 counts (if existing) are very valuable in this case. 
For _all on the other hand, reads that mapped exactly twice were equally divided between the two positions. 
This gives more reliable peak shapes but the 0/0 values do not necessarily display the real situation. 
Since _all has more reads, CPM values are lower and we have two different cut-offs for peaks: *0.17 for _all and 0.3 for _uni.*
The file structure is the same as for the previous data. 
Columns 1-4 contain the general information for each SNP: 
Chromosome(B73) Position(B73) B73Allele NAMAllele; 
the remaining columns contain, for each NAM line, ID Genotype CPM(B73) CPM(NAM) Peak(B73) Peak(NAM) and PostFrequency. 
Peak is 0 if below the threshold and 1 if above. 
PostFrequency is (B73/(B73+NAM)), thus 1 for total bias towards B73 and 0 for total bias to NAM. 
If both alleles are below the peak threshold, n.p. is place before the PF value. 
If a SNP is absent in a NAM, all fields are filled with n.a. for this line.

*Check the new readme to see what is specifically different*