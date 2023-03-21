#Notes about MOA 
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

##NAM RIL DATA
NAM SNP projections on the RIL population `/iplant/home/shared/NAM/Misc/NAM-SV-projected-V8`
```
ml singularity
singularity shell --bind $PWD /home/snodgras/irods.simg
iget /iplant/home/shared/NAM/Misc/NAM-SV-projected-V8/*SNPs-only*
```
Phenotype data from: Merritt Burch `all_NAM_phenos.csv`

_cleaning up datasets to be converted to PLINK format_
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

_Geting the MOA Data_
```
#should really do this on the dtn and a screen 
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

### Testing VCAP with MOA files

```
#create the genotype files (see concatenate_files.sh)
#check to make sure the headers are all in the same order
for i in {1..10}; do zcat NAM_rils_projected-reseq-SNPs-only.all-RILs.chr-${i}.v8.hmp.txt.gz | head -n 1 > header-chr${i}.txt ; done
for i in {1..10}; do 
	for j in {1..10} ;do 
		echo Comparing chr${i} and chr${j} ; 
		cmp --silent header-chr${i}.txt header-chr${j}.txt && echo "Files are identical" || echo "Files are different" ; 
	done ; 
done
zcat NAM_rils_projected-reseq-SNPs-only.all-RILs.chr-1.v8.hmp.txt.gz > NAM_rils_SNPs_allchr.hmp.txt
for i in {2..10}; do zcat NAM_rils_projected-reseq-SNPs-only.all-RILs.chr-${i}.v8.hmp.txt.gz | tail -n +2 >> NAM_rils_SNPs_allchr.hmp.txt ; done
module load tassel
module load samtools
bgzip -c NAM_rils_SNPs_allchr.hmp.txt > NAM_rils_SNPs_allchr.hmp.txt.gz
#Create the LIX index:
run_pipeline.pl -Xmx64g -LIXPlugin -createIndex NAM_rils_SNPs_allchr.hmp.txt.gz


#create test file for phenotypes
cut -f 1-15 all_NAM_phenos.plink > testsubset_NAM_phenos.plink

#create test file for peaks
for i in *WW.peaks.names.sorted.bed ; do grep ^B73 $i >> B73_WW_peaks_allhybrids.bed ; done
# peak ID's are not shared between files
cut -f 1-3 B73_WW_peaks_allhybrids.bed | sort | uniq > B73_WW_uniqpeaks.bed
sed -i "s/B73-//g" B73_WW_uniqpeaks.bed 
```

Running the test VCAP on MOA peaks
```
module load tassel
run_pipeline.pl \
	-Xmx64g \
	-importGuess NAM-ril-genotypes/NAM_rils_SNPs_allchr.hmp.txt.gz.lix  \
	-KinshipPlugin \
	-method Centered_IBS \
	-endPlugin \
	-export NAMRILs_projectedSNPs_CentIBS_wholegenomekinship.txt.gz \
	-exportType SqrMatrix

run_pipeline.pl -Xmx64g \
        -fork1 -importGuess NAM-ril-genotypes/NAM_rils_SNPs_allchr.hmp.txt.gz.lix \
        -fork2 -FilterSiteBuilderPlugin -bedFile MOA_peak_calls/Peaks/EG80/WW/EG80.WW.all_peaks.merged.inclCML333.bed -endPlugin \
        -input1 -KinshipPlugin -method Centered_IBS -endPlugin \
        -fork3 -export EG80WWall_CentIBS_kinship1 -exportType SqrMatrixBin -input2 \
        -fork4 -SubtractDistanceMatrixPlugin -wholeMatrix NAMRILs_projectedSNPs_CentIBS_wholegenomekinship.txt.gz -endPlugin \
        -input2 -export EG80WWall_CentIBS_kinshipRest -exportType SqrMatrixBin

module load ldak
ldak --reml LDAK_EG80WWall --mgrm kinship_list.txt --mpheno -1 --pheno NAM-ril-phenotypes/testsubset_NAM_phenos.plink --kinship-details NO

```

Errors running LDAK

```
Error, Z003E0046___Z003E0046 appears twice in NAM-ril-phenotypes/testsubset_NAM_phenos.plink
#solution
uniq testsubset_NAM_phenos.plink > uniq_testsubset_NAM_phenos.plink
#problem
Error, to analyse multiple phenotypes, each sample must have either all phenotypes present or all missing; you should instead either analyse each phenotype separately, or use "--dentist YES" to pad missing values
```

`VCAP_test_MOAEG80WW_5.log` is the working LDAK command without `--constrain YES`
`VCAP_test_MOAEG80WW_6.log` is the final working log

Getting a summary from all the different phenotypes

```
for i *.reml ; do
	n=$(echo $i | cut -f 2 -d \. )
	let s=$n+2
	head -n 1 NAM-ril-phenotypes/testsubset_NAM_phenos.plink | cut -f $s >> Heritability_All.log
	tail -n 5 $i >> Heritability_All.log 
done
```
Creating a file that is easy to make an excel chart with

```
#name of phenotype | Heritability_K1 | Her_SD_K1 | Heritability_K2 | Her_SD_K2 | Heritability_All | Her_SD_All
grep -n Component FullEG80_heritability.log | cut -d : -f 1 >> linenumbers.txt
while read -r line ; do
	awk -v linenumber=$line -v IFS=" " -v OFS="\t" '{if(NR==linenumber-1 || NR==linenumber+1 || NR==linenumber+2 || NR==linenumber+4) print $1,$2,$3}' FullEG80_heritability.log >> FullEG80_Results.txt ; 
done < linenumbers.txt

```

To scrape the info out of the `*.reml` files, use `format_hests.sh`

## Creating the background distribution

1. Find the distance between MOA peaks and B73v5 gene models (``/work/LAS/mhufford-lab/snodgras/NAM_FT_genes/RNAseq_bam_files/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff`)

```
module load bedops
gff2bed gff_file > gene.bed

#need to separate out only the gene coordinates (not exon, chromosome, mRNA, etc.)
awk -v OFS='\t' '$8 ~ /gene/ {print $0}' sorted.gff.bed > sorted.gene.bed

#sort the peak bedfile
sort -k1,1 -k2,2n in.bed > in.sorted.bed

module load bedtools2
bedtools closest -D ref -t first -a in.sorted.bed -b gff

#slurm command

bash gene_distance.sh Peaks/EG80/WW/EG80.WW.all_peaks.merged.inclCML333.bed sorted.gene_Zm-B73v5REF.bed

```

Now binning the calculated distances

```
awk -v OFS='\t' '$1 <= -1000000  {print $1, "more_than_1000000_upstream"; }
	$1 <= -100000 && $1 > -1000000 {print $1, "more_than_100000_upstream"; }
	$1 <= -10000 && $1 > -100000 {print $1, "more_than_10000_upstream"; }
	$1 <= -1000 && $1 > -10000 {print $1, "more_than_1000_upstream"; }
	$1 <= -100 && $1 > -1000 {print $1, "more_than_100_upstream"; }
	$1 <= -10 && $1 > -100 {print $1, "more_than_10_upstream"; }
	$1 <= -1 && $1 > -10 {print $1, "more_than_1_upstream"; }
	$1 == 0 {print $1, "overlaps_with_gene"; }
	$1 >= 1 && $1 < 10 {print $1,"more_than_1_downstream" ; }
	$1 >= 10 && $1 < 100 {print $1,"more_than_10_downstream" ; }
	$1 >= 100 && $1 < 1000 {print $1,"more_than_100_downstream" ; }
	$1 >= 1000 && $1 < 10000 {print $1,"more_than_1000_downstream" ;}
	$1 >= 10000 && $1 < 100000 {print $1,"more_than_10000_downstream" ; }
	$1 >= 100000 && $1 < 1000000 {print $1,"more_than_100000_downstream" ; }
	$1 >= 1000000 {print $1,"more_than_1000000_downstream" ;}' peak_distances.txt > binned_peak_distances.txt
	
cut -f 2 binned_peak_distances.txt | sort | uniq -c
	 59 more_than_1000000_downstream
     85 more_than_1000000_upstream
  21527 more_than_100000_downstream
  20299 more_than_100000_upstream
 151673 more_than_10000_downstream
 151172 more_than_10000_upstream
 140928 more_than_1000_downstream
 140181 more_than_1000_upstream
  73448 more_than_100_downstream
  72762 more_than_100_upstream
   9025 more_than_10_downstream
   8872 more_than_10_upstream
    937 more_than_1_downstream
   5206 more_than_1_upstream
 341031 overlaps_with_gene
```

2. Calculate the length distribution of the peaks

```
awk -v OFS='\t' '{ print $3-$2 }' Peaks/EG80/WW/EG80.WW.all_peaks.merged.inclCML333.bed > peak_lengths.txt
awk -v OFS='\t' '$1 <10 {print $1, "less_than_10" ; }
	$1 >= 10 && $1 <25 {print $1,"between_10_24"; }
	$1 >= 25 && $1 <50 {print $1,"between_25_50"; }
	$1 >= 50 && $1 <100 {print $1,"between_50_100"; }
	$1 >= 100 && $1 <500 {print $1,"between_100_500"; }
	$1 >= 500 && $1 <1000 {print $1,"between_500_1000"; }
	$1 >= 1000 {print $1,"more_than_1000" ;}' peak_lengths.txt > binned_peak_lengths.txt

cut -f 2 binned_peak_lengths.txt | sort | uniq -c 
 
 102833 less_than_10
 249061 between_10_24
 258458 between_25_50
 328058 between_50_100
 182593 between_100_500 
  14540 between_500_1000
   1662 more_than_1000
#1137205 total peaks
```
Finding median length for each bin
```
grep less_than_10 MOA_peak_calls/binned_peak_lengths.txt | sort -k1,1n | awk -f ~/bin/median.awk
4
```
3. Create a script that subsets the genome by those bins and samples from them

```
#make chromosome bed file
grep chromosome Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.bed | awk -v OFS='\t' '{print $1,$3}' - > chromosome_Zm-B73v5REF.bed
#shuffle the peak windows
bedtools shuffle -i Peaks/EG80/WW/EG80.WW.all_peaks.merged.inclCML333.bed -g chromosome_Zm-B73v5REF.bed -noOverlapping -excl Peaks/EG80/WW/EG80.WW.all_peaks.merged.inclCML333.bed > shuffled_peak_windows.bed
#find the distances between shuffled windows and closest gene
sort -k1,1 -k2,2n shuffled_peak_windows.bed | bedtools closest -D ref -t first -a - -b sorted.gene_Zm-B73v5REF.bed > shuffled_distances.txt
cut -f 14 shuffled_distances.txt | awk -v OFS='\t' '$1 <= -1000000  {print $1, "more_than_1000000_upstream"; }
$1 <= -100000 && $1 > -1000000 {print $1, "more_than_100000_upstream"; }
$1 <= -10000 && $1 > -100000 {print $1, "more_than_10000_upstream"; }
$1 <= -1000 && $1 > -10000 {print $1, "more_than_1000_upstream"; }
$1 <= -100 && $1 > -1000 {print $1, "more_than_100_upstream"; }
$1 <= -10 && $1 > -100 {print $1, "more_than_10_upstream"; }
$1 <= -1 && $1 > -10 {print $1, "more_than_1_upstream"; }
$1 == 0 {print $1, "overlaps_with_gene"; }
$1 >= 1 && $1 < 10 {print $1,"more_than_1_downstream" ; }
$1 >= 10 && $1 < 100 {print $1,"more_than_10_downstream" ; }
$1 >= 100 && $1 < 1000 {print $1,"more_than_100_downstream" ; }
$1 >= 1000 && $1 < 10000 {print $1,"more_than_1000_downstream" ;}
$1 >= 10000 && $1 < 100000 {print $1,"more_than_10000_downstream" ; }
$1 >= 100000 && $1 < 1000000 {print $1,"more_than_100000_downstream" ; }
$1 >= 1000000 {print $1,"more_than_1000000_downstream" ;}' - > binned_shuffled_distances.txt
cut -f 2 binned_shuffled_distances.txt | sort | uniq -c

```
### Working with SNPs rather than peaks
1. Finding the distances between genome alignment SNPs and genes
Trying to understand why there are so many more peaks in the full peak list than what Thomas says
(~1 million vs 200/250K)

```
#EG80WW_peaklabeled_dist.txt was made later, but is almost the same as Peaks/EG80/WW/EG80.WW.all_peaks.merged.inclCML333.bed (same number entries and coordinates)
sort -k1,1 -k2,2n EG80WW_peaklabeled_dist.txt | bedtools merge -i - -c 5 -o distinct > testing_merge_theory.txt #gave same number of lines, no change after "merging"
sort -k1,1 -k2,2n EG80WW_peaklabeled_dist.txt | bedtools merge -i - > testing_merge_theory.txt# same number of lines, no change after merging (wc -l ==> 1137205 )
```

```
awk -v OFS='\t' '{print $1,$2,$3,$25}' Summary_genotypes_AF.bed > Summary_MOASNPs_AF.bed
sed -i 's/B73\-//g' Summary_MOASNPs_AF.bed
```
2. Intersect these alignment SNPs with the MOA peak regions to remove MOA SNPs 
```
awk -v OFS='\t' '{print $1,$2,$3, $4,$15}' out.distance.Summary_MOASNPs_AF.bed.B73v5.bed > alignSNP_AF_dist.txt
awk -v OFS='\t' '{print $1,$2,$3,$14}' out.distance.EG80.WW.all_peaks.merged.inclCML333.bed.B73v5.bed  > EG80WWpeak_dist.txt
module load bedtools2
bedtools intersect -v -a alignSNP_AF_dist.txt -b EG80WWpeak_dist.txt > nonEG80WWSNPs_AF_dist.txt
```
and also to make a MOA SNP AF Dist file
```
awk -v OFS='\t' '{print $0, "EG80WW_peak_"NR}' EG80WWpeak_dist.txt > EG80WW_peaklabeled_dist.txt
module load bedtools2
bedtools intersect -wb -a alignSNP_AF_dist.txt -b EG80WW_peaklabeled_dist.txt > MOA_EG80WWSNPs_AF_dist.txt
```

### How to sample a background distribution to match the AF of MOA alleles and distance from genes
The `nonEG80WWSNPs_AF_dist.txt` has the format:
chromosome \t start \t stop(actual SNP location) \t AF \t distance from nearest gene
```
bin the SNP coordinates by distance into files
figure out how many SNPs of each AF for each distance is needed
Randomly pick until those conditions are satisfied
ta-da background SNPs!
But how to define coordinates for the VCAP? Would want to keep the same bp in each component...

awk -v OFS='\t' '$5 <= -1000000 || $5 >= 1000000 {print $0 >> "1MB_distant_SNPs.txt" ;}
$5 > -1000000 && $5 <= -100000 {print $0 >> "100KB_distant_SNPs.txt";}
$5 < 1000000 && $5 >= 100000 {print $0 >> "100KB_distant_SNPs.txt";}
$5 > -100000 && $5 <= -10000 {print $0 >> "10KB_distant_SNPs.txt";}
$5 < 100000 && $5 >= 10000 {print $0 >> "10KB_distant_SNPs.txt";}
$5 > -10000 && $5 <= -1000 {print $0 >> "1KB_distant_SNPs.txt";}
$5 < 10000 && $5 >= 1000 {print $0 >> "1KB_distant_SNPs.txt";}
$5 > -1000 && $5 <= -100 {print $0 >> "100bp_distant_SNPs.txt";}
$5 < 1000 && $5 >= 100 {print $0 >> "100bp_distant_SNPs.txt";}
$5 > -100 && $5 <= -10 {print $0 >> "10bp_distant_SNPs.txt";}
$5 < 100 && $5 >= 10 {print $0 >> "10bp_distant_SNPs.txt";}
$5 > -10 && $5 <= -1 {print $0 >> "1bp_distant_SNPs.txt";}
$5 < 10 && $5 >= 1 {print $0 >> "1bp_distant_SNPs.txt";}
$5 == 0 {print $0 >> "genic_distant_SNPs.txt"}' nonEG80WWSNPs_AF_dist.txt 
```

3. Parse out the different AF between each distance file
`binbyAF.sh`
```
input=$1
for i in $input ; do 
	awk -v OFS='\t' '$4 <0.05 {print $0 >> "AF_lt05.txt";}
	$4 >= 0.05 && $4 < 0.1 {print $0 >> "AF_05-10.txt";}
	$4 >= 0.1 && $4 < 0.15 {print $0 >> "AF_10-15.txt";}
	$4 >= 0.05 && $4 < 0.1 {print $0 >> "AF_05-10.txt";}
	$4 >= 0.15 && $4 < 0.2 {print $0 >> "AF_15-20.txt";}
	$4 >= 0.2 && $4 < 0.25 {print $0 >> "AF_20-25.txt";}
	$4 >= 0.25 && $4 < 0.3 {print $0 >> "AF_25-30.txt";}
	$4 >= 0.3 && $4 < 0.35 {print $0 >> "AF_30-35.txt";}
	$4 >= 0.35 && $4 < 0.4 {print $0 >> "AF_35-40.txt";}
	$4 >= 0.4 && $4 < 0.45 {print $0 >> "AF_40-45.txt";}
	$4 >= 0.45 && $4 < 0.5 {print $0 >> "AF_45-50.txt";}
	$4 >= 0.5 && $4 < 0.55 {print $0 >> "AF_50-55.txt";}
	$4 >= 0.55 && $4 < 0.6 {print $0 >> "AF_55-60.txt";}
	$4 >= 0.6 && $4 < 0.65 {print $0 >> "AF_60-65.txt";}
	$4 >= 0.65 && $4 < 0.7 {print $0 >> "AF_65-70.txt";}
	$4 >= 0.7 && $4 < 0.75 {print $0 >> "AF_70-75.txt";}
	$4 >= 0.75 && $4 < 0.8 {print $0 >> "AF_75-80.txt";}
	$4 >= 0.8 && $4 < 0.85 {print $0 >> "AF_80-85.txt";}
	$4 >= 0.85 && $4 < 0.9 {print $0 >> "AF_85-90.txt";}
	$4 >= 0.9 && $4 < 0.95 {print $0 >> "AF_90-95.txt";}
	$4 >= 0.95 && $4 < 1 {print $0 >> "AF_95-100.txt";}' $i
	
	for j in AF*.txt ; do
		mv $j ${i%_SNPs.txt}_${j} ;
	done
done
```

4. Pick 1 SNP from each peak region and match to 1 non-peak SNP

a. randomly pick 1 SNP from each peak
b. Find AF and distance from closest gene for picked SNPs
```
cut -f 10 MOA_EG80WWSNPs_AF_dist.txt | sort | uniq > peak.list
while read -r line ; do
	grep $line MOA_EG80WWSNPs_AF_dist.txt > SNPsinpeak #gets all the SNPs for a given peak
	samplesize=$(wc -l SNPsinpeak) #what's the upper limit for a random number
	RNV=$(shuf -i 1-${samplesize} -n 1) #pick a random number, RNV=random number variable 
	awk -v OFS='\t' -v rnv=$(echo RNV) 'NR == rnv {print $5,$4}' MOA_EG80WWSNPs_AF_dist.txt  >> MOA_SNP_distrib.txt #cut the dist col and AF col for picked SNP and add to ref file
done < peak.list #feed in the peak names
```
c. match those criteria to non-peak SNPs
```
while read -r line ; do
	dist=$(cut -f 1 $line)
	AF=$(cut -f 2 $line)
	awk -v OFS='\t' -v d=$(echo $dist) -v af=$(echo $AF) '$5 == d && $4 == af {print $0}' nonEG80WWSNPs_AF_dist.txt >> pool
	samplesize=$(wc -l pool) #what's the upper limit for a random number
	RNV=$(shuf -i 1-${samplesize} -n 1) #pick a random number, RNV=random number variable 
	awk -v OFS='\t' -v rnv=$(echo $RNV) 'NR == rnv {print $0}' pool  >> bkgd_SNP_distrib.txt #print picked background SNP record into a background distribution file
done < MOA_SNP_distrib.txt
```
wrote up the above as a script and testing it as
```
head -n 1000 MOA_EG80WWSNPs_AF_dist.txt > test_MOAsnps.txt
head -n 5000 nonEG80WWSNPs_AF_dist.txt > test_BKGDsnps.txt
bash SNPselection.sh test_MOAsnps.txt test_BKGDsnps.txt
```

##Using genic background
Match only on AF
```
awk -v OFS="\t" '{print $4,$5,$6,$7}' out.distance.EG80.WW.all_peaks.merged.inclCML333.bed.B73v5.bed | sort -k4,4 | uniq > nearEG80WWMOApeak.B73v5gene.bed
#manually remove the first line
```
I can just use the nonMOA peak file and filter for gene distance == 0 --> this gets ANY genes not genes NEAR MOA peaks
To get genes NEAR MOA peaks, I'll intersect `nearEG80WWMOApeak.B73v5gene.bed` with `nonEG80WWSNPs_AF_dist.txt`

```
bedtools intersect -a nonEG80WWSNPs_AF_dist.txt -b nearEG80WWMOApeak.B73v5gene.bed > nonEG80WWSNPs_nearMOApeak.txt
```

```
mkdir geneKBKGD_Test
for i in {0..99}; do mkdir perm_$i ; done
#add in the R script for selection
chmod +x SNPselection_GeneBackground.R
for i in {0..99}; do 
	echo module load r r-devtools \; module load r-tidyr \; module load r-dplyr \; ./SNPselection_GeneBackground.R ../MOA_EG80WWSNPs_AF_dist.txt ../nonEG80WWSNPs_nearMOApeak.txt perm_${i}/coords_ >> commands.txt ; 
done
python ../makeSLURM.py 1 commands.txt
#run commands_0.sub as a test
for i in {1..99} ; do sbatch commands_${i}.sub ; done
```
Then run `slurm_runLDAK.sh` and then `format_hests.sh`

## Testing order of terms in LDAK
Need to change the order of the kinship list file
Using the 1 SNP/peak selection scheme, matched background on AF and distance to nearest gene

Alterations to `slurm_runLDAK.sh`

```
echo NAMRILs_projectedSNPs_CentIBS_wholegenomekinship.txt.gz_Rest_K > kinship.list
echo BKGD.SNPs.bed_CentIBS_K >> kinship.list
echo MOA.SNPs.bed_CentIBS_K >> kinship.list

bash /ptmp/LAS/snodgras/VCAP_test/VCAP_runLDAK.sh RBM_model kinship.list /ptmp/LAS/snodgras/VCAP_test/NAM-ril-phenotypes/uniq_height_NAM_phenos.plink
```

## Testing just matching gene distance
Re-write SNPselection to only match on distance `SNPselection_MatchDistOnly.R`
```
for i in {0..99}; do
	mkdir perm_$i ;  
	echo module load r r-devtools \; module load r-tidyr \; module load r-dplyr \; ./SNPselection_MatchDistOnly.R ../MOA_EG80WWSNPs_AF_dist.txt ../nonEG80WWSNPs_nearMOApeak.txt perm_${i}/coords_ >> commands.txt ; 
done

```

## Permutation of VCAP

### 0 : Create starting files
_Create allele frequency and distance containing SNP bed files for MOA peak and non-peak regions_
_If using new markers, make the full kinship matrix_

### 1 : Run background selection (permutation starts)
_Select a matched background_
	- run `SNPselection.R`
	- make sure the resulting files are being written to specific directories

commands to run `SNPselection.R`
```
module load r r-devtools r-tidyr r-dplyr
./SNPselection.R $MOA_AF_dist.txt $nonMOASNPs_AF_dist.txt $outputdirectory/filename
```
Use `makeSLURM.py` to be able to run it in parallel
```
#make a command file
for i in {0..99} ; do 
	mkdir EG80WW_perm_$i
	echo module load r r-devtools r-tidyr r-dplyr \; ./SNPselection.R MOA_EG80WWSNPs_AF_dist.txt  nonEG80WWSNPs_AF_dist.txt EG80WW_perm_$i/coords_$i_ >> 1_command.txt ;
done
python makeSLURM.py 1 1_command.txt
```
	
### 2 : Run VCAP model 
_Make the kinship matrices needed based on selected SNPs in 1_
_Requires a full kinship matrix of all markers (see 0)_
_make sure output writes to specific directory_

```
#add in slurm script to each directory
for i in {2..99} ; do 
	cp EG80WW_perm_1/slurm_VCAP_mkKinship.sh EG80WW_perm_$i/.
	cd EG80WW_perm_$i/ 
	sbatch slurm_VCAP_mkKinship.sh
	cd ..
done

for i in {1..99} ; do 
	rm EG80WW_perm_$i/MBrest_model* ; 
done
```
To find the specific directory of a failed job:
```
for i in EG80WW_perm* ; do if [[ -f $i/slurm-152619.out ]] ; then echo $i ; fi ; done
```
perm 88: had trouble writing the final kinship matrix (rest) (ERROR: bad address)
same error for perm 54
removing the created kinship files and rerunning the `slurm_VCAP_mkKinship.sh` didn't throw an error

perm 39: `ERROR net.maizegenetics.plugindef.AbstractPlugin - filterSitesByBedFile: problem reading: MOA.SNPs.bed line: ch9925176`
Doesn't appear to be a carry over from the `coords.txt` file so it's likely an error in making the bed file
Tried `tail -n +2 coords_MOA_distrib.SNPs.txt | cut -f 1-3 | grep ch992517 -` and didn't get anything; rerunning and seeing if its a repeat error.

Creating an extra phenotype test file having only the height phenotypes

```
awk -v OFS='/t' '{print $1,$2,$16,$17,$21,$22,$62,$63}' NAM-ril-phenotypes/uniq_all_NAM_phenos.plink > NAM-ril-phenotypes/uniq_height_NAM_phenos.plink
```
Run the LDAK part
```
#add in slurm script to each directory
for i in {2..99} ; do 
	cp EG80WW_perm_1/slurm_VCAP_runLDAKp.sh EG80WW_perm_$i/.
	cd EG80WW_perm_$i/ 
	sbatch slurm_VCAP_mkKinship.sh
	cd ..
done
```

### 3 : Create output summary
_Make sure the final file is written to the general directory (otherwise will be deleted)_
_Will delete all the rest of the files_
	- run `format_hests.sh`
	
### For testing biased peak

```
cut -f 1-3 CML277.TT.01.expr.bed > Test_biasedSNP.bed
sed -i '' -e 's/B73\-//g' Test_biasedSNP.bed

module load bedtools2
bedtools intersect -wb -a ../MOA_EG80WWSNPs_AF_dist.txt  -b Test_biasedSNP.bed | cut -f 10 | sort | uniq > biased_peaks.list
grep -w -f biased_peaks.list ../MOA_EG80WWSNPs_AF_dist.txt | cut -f 1-5,10 >> biased_EG80WW.bed
bedtools intersect -v -a ../MOA_EG80WWSNPs_AF_dist.txt -b biased_EG80WW.bed | cut -f 1-5,10 > unbiased_EG80WW.bed

```

Using the new file from Julia `Significant_biased_peaks_16NAM_limma_05.bed`

```
sed -i '' -e 's/B73\-//g' Significant_biased_peaks_16NAM_limma_05.bed
#manually delete last line "ID	-1"

bedtools intersect -a MOA_EG80WWSNPs_AF_dist.txt  -b Significant_biased_peaks_16NAM_limma_05.bed > biased_EG80WWsnps.txt
bedtools intersect -v -a MOA_EG80WWSNPs_AF_dist.txt -b Significant_biased_peaks_16NAM_limma_05.bed > unbiased_EG80WWsnps.txt
```

To do the SNP selection commands in parallel:

```
for i in {0..99}; do mkdir perm_$i ; done
for i in {0..99}; do 
	echo module load r r-devtools \; module load r-tidyr \; module load r-dplyr \; ./SNPselection_biasedpeaks.R ../biased_EG80WWsnps.txt ../unbiased_EG80WWsnps.txt perm_${i}/coords_ >> commands.txt ; 
done
python ../makeSLURM.py 1 commands.txt
#run commands_0.sub as a test
for i in {1..99} ; do sbatch commands_${i}.sub ; done

# checking for errors in the kinship runs
for i in {0..99} ; do echo $i ; grep ERROR perm_${i}/slurm*.out ; done
```

##Jeff's test
Simulated Phenotypes
Used Kinship matrices from `1SNPperPeak_Originalperms/EG80WW_perm_1/`
```
paste Example_Kinships_1SNPperPeak/BKGD.SNPs.bed_CentIBS_K.grm.id phenos.y.txt > JeffsPhenos.plink
sed -i 's/ /\t/g' JeffsPhenos.plink
```
For the second test, I removed the Rest kinship name from the kinship list so it should hopefully only run with the BKGD and MOA kinships

For the third test needed
```
cut -d " " -f 2- phenos_round8.txt | paste Example_Kinships_1SNPperPeak/BKGD.SNPs.bed_CentIBS_K.grm.id - > JeffsPhenos.3.plink
sed -i 's/ /\t/g' ../JeffsPhenos.3.plink 
```

Make sure the plink format is FID IID Trait names
```
#head -n 1 phenos.txt | echo "FID IID" - > phenos.plink
tail -n +2 Simulated_Phenotypes_highMOAh2.txt | paste Example_Kinships_1SNPperPeak/BKGD.SNPs.bed_CentIBS_K.grm.id -  >> simulatedphenotypes_highMOAh2.plink
sed -i 's/ /\t/g' phenos.plink
```
On the projected SNPs, I realized that the same kinship matrix was being loaded in for every matrix
So I fixed that so the different K matrices were loaded in correctly.
Made sure that the randomization was in so it wasn't the same trait value calculated over and over.
And created `Simulated_Phenotypes.Round2.txt`
To convert to plink format:
```
head -n 1 Simulated_Phenotypes.Round2.txt > Simulated_Phenotypes.Round2.plink
#manually add in the "FID \t IID"
tail -n +2 Simulated_Phenotypes.Round2.txt | paste BKGD.SNPs.bed_CentIBS_K.grm.id -  >> Simulated_Phenotypes.Round2.plink
sed -i '' -e 's/ /\t/g' Simulated_Phenotypes.Round2.plink
```

*CAN PROBABLY REMOVE EVERYTHING AFTER THIS LINE*
###Biased Peak Differential Expression Analysis

Arun reformatted the vcf and sent as `clean-file.tsv`
I renamed it to `NAM_SNPgenotypes.tsv`
Julia sent this biased peak file `CML277.TT.01.expr.bed` and I renamed to `Biased_MOAsnps.tsv`

1. filter to only get snps from biased file

```

sed -i '' -e 's/B73\-//g' Biased_MOAsnps.tsv
awk -v OFS='\t' '{print $1, $2-1, $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31}' NAM_SNPgenotypes.tsv > NAM_SNPgenotypes.bed
module load bedtools2
bedtools intersect -wa -wb -a NAM_SNPgenotypes.bed -b Biased_MOAsnps.tsv > Biased_SNPintersect.tsv

```

2. filter that to only get snps from +/- 3kb a gene

```
awk -v OFS='\t' '$NF > -3001 && $NF < 3001 {print $0}' Biased_SNPintersect.tsv > 3KB_biased_SNPintersect.tsv
awk -v OFS="\t" '{print $(NF-1)}' 3KB_biased_SNPintersect.tsv | sort | uniq > biased_gene.list
```

3. filter tpm to get only those gene ids/expression levels

```
ml singularity
singularity shell --bind $PWD /home/snodgras/irods.simg
iinit
icd /iplant/home/shared/NAM/Misc/to_Yinjie
iget NAM_pangene_expression_counts_per_tissue-TPM.tar.gz
tar xvzf NAM_pangene_expression_counts_per_tissue-TPM.tar.gz

grep -f biased_gene.list NAM_pangene_expression_counts_per_tissue-TPM/B73_full-tpm.tsv | cut -f 1-2 > biased_pangene.list
#should have 1842 lines if the grep worked
cut -f 2 biased_pangene.list|cut -d _ -f 1 | grep -v -f - biased_gene.list | cut -f 1-2 > missing_biased_pangene.list

#check that the number of pangene ids is == to what comes out of the grep 

for i in *full-tpm.tsv; do
	head -n 1 $i > biased_gene_${i%full-tpm.tsv}-tpm.tsv ;
	cut -f 1 ../biased_pangene.list | grep -w -f - $i >> biased_gene_${i%full-tpm.tsv}-tpm.tsv ; 
done

```

4. TPM ~ gene + MOAsnpGenotype (or separate tests for each gene?)

## Projecting the Alignment SNPs onto the NAM RILs

###97M vs GATK overlap
*How many SNPs in the GATK?*

```
module load bedtools2
bedtools intersect -wa -a /ptmp/LAS/snodgras/VCAP_test/MOA_peak_calls/Summary_genotypes_AF.bed -b /ptmp/LAS/snodgras/MOA_DE_Analysis/NAM_SNPgenotypes.bed > AlignGATKoverlap.out
```
`AlignGATKoverlap.out` is 27750427 lines long

*How many peaks in each?*

```
#number of peaks in overlap (747804)
bedtools intersect -wb -a AlignGATKoverlap.out -b /ptmp/LAS/snodgras/VCAP_test/MOA_peak_calls/EG80WW_peaklabeled_dist.txt > AlignGATKoverlap_EG80WWpeaklabeled.out
awk -v OFS="\t" '{print $NF}' AlignGATKoverlap_EG80WWpeaklabeled.out | sort | uniq | wc -l 

#number of peaks called (1137205)
awk -v OFS="\t" '{print $NF}' EG80WW_peaklabeled_dist.txt | sort | uniq | wc -l

#number of peaks in the 97M call (953628)
awk -v OFS="\t" '{print $NF}' MOA_EG80WWSNPs_AF_dist.txt | sort | uniq | wc -l

#number of 97M SNPs within peaks (9004469)
wc -l MOA_EG80WWSNPs_AF_dist.txt

#number of overlap SNPs within peaks (3556784)
wc -l AlignGATKoverlap_EG80WWpeaklabeled.out

#number of shared SNPs between 97M call and overlap set
echo "Number of values of overlap in 97M call"
bedtools intersect -a /ptmp/LAS/snodgras/VCAP_test/MOA_peak_calls/MOA_EG80WWSNPs_AF_dist.txt -b AlignGATKoverlap.out | wc -l 
echo "Number of values of 97M in overlap call"
bedtools intersect -a AlignGATKoverlap.out -b /ptmp/LAS/snodgras/VCAP_test/MOA_peak_calls/MOA_EG80WWSNPs_AF_dist.txt  | wc -l 
#Answer was the same both ways: 3,556,784

#number of shared peaks(747,804)
bedtools intersect -a /ptmp/LAS/snodgras/VCAP_test/MOA_peak_calls/MOA_EG80WWSNPs_AF_dist.txt -b AlignGATKoverlap.out > test.txt
awk -v OFS="\t" '{print $NF}' test.txt | sort | uniq | wc -l  

```

From `http://augustogarcia.me/statgen-esalq/Hapmap-and-VCF-formats-and-its-integration-with-onemap/`

```
sed -i 's/B73\-//g' Summary_genotypes_AF.bed
awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' Summary_genotypes_AF.bed > Alignment_SNPs.txt
python reformat2hapmap.py

awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' AlignGATKoverlap.out > Overlap_SNPindel.txt
python reformat2hapmap.py > Overlap_SNPindel.hmp.txt

```
Manually changed the header line, verify that it's still tab delimited and that the columns are in the right order

Split into crosses and chromosomes
```
for i in (list of NAM names or column numbers); do
	cut -f 1-generic hmp columns,$i > B73x${i}.overlap.hmp.txt
done

for i in {1..10} ; do
	c=$(echo chr${i})
	head -n 1 B73x...overlap.hmp.txt > B73x...overlap.${c}.hmp.txt
	grep -w $c B73x...overlap.hmp.txt >> B73x...overlap.${c}.hmp.txt
done
```
###63M vs GATK overlap
*How many SNPs in the GATK?*
`NAM_SNPgenotypes.bed` is 62224083 lines long
`All_SNPs_bialelic_chr1_10_new_04_2021_NAM_genotypes_cleaned_0423.bed` is 64787967 lines ling

Overlap with GATK that has indels
```
module load bedtools2
bedtools intersect -wa -a /ptmp/LAS/snodgras/All_SNPs_bialelic_chr1_10_new_04_2021_NAM_genotypes_cleaned_0423.bed -b /ptmp/LAS/snodgras/MOA_DE_Analysis/NAM_SNPgenotypes.bed > 0423_21_AlignGATKindeloverlap.out
```
`0423_21_AlignGATKindeloverlap.out` is 27849652 lines long

Overlap with GATK that doesn't have indels (`nova:/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/l-gatk-illumina/c-gatk/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.gz`)
```
module load bedtools2
bedtools intersect -wa -a /ptmp/LAS/snodgras/All_SNPs_bialelic_chr1_10_new_04_2021_NAM_genotypes_cleaned_0423.bed -b /ptmp/LAS/arnstrm/to-sam/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf > 0423_21_AlignGATK_NOindel_overlap.out
```
`0423_21_AlignGATK_NOindel_overlap.out` is 17655201 lines long

To get overlap with Peaks
```
##Calculate AF
## count number of columns 6:30 that == 0/0 (B73_allele), == 1/1 (NAM_allele), or == n.a.
## AF = # nam allele/(# nam allele + b73 allele)

gawk -v OFS='\t' '{if(NR==1){print $0, "Number_NAMalleles","Number_present","NAM_AF"};if(NR>1){for(j=0;j<=30;j=j+1){if($(6+j)=="1/1"){NAL++};if($(6+j)=="0/0" || $(6+j)=="1/1"){Pres++}};{if(NAL>0){print $0,NAL,Pres,NAL/Pres }};NAL=0;Pres=0}}' All_SNPs_bialelic_chr1_10_new_04_2021_NAM_genotypes_cleaned_0423.bed > 0423_alignSNPs_AF.bed
```
```
bash gene_distance.sh ZmB73v5.gff All_SNPs_bialelic_chr1_10_new_04_2021_NAM_genotypes_cleaned_0423.bed 

awk -v OFS='\t' '{print $1,$2,$3,$41}' out.distance.All_SNPs_bialelic_chr1_10_new_04_2021_NAM_genotypes_cleaned_0423.bed.ZmB73v5.gff.bed > 0423_alignSNPs_dist.txt
awk -v OFS='\t' '{print $1,$2,$3,$14}' out.distance.EG80.WW.all_peaks.merged.inclCML333.bed.B73v5.bed  > EG80WWpeak_dist.txt
module load bedtools2
#bedtools intersect -v -a 0423_alignSNPs_dist.txt -b EG80WWpeak_dist.txt > nonEG80WWSNPs_dist.txt
```
and also to make a MOA SNP AF Dist file
```
awk -v OFS='\t' '{print $0, "EG80WW_peak_"NR}' EG80WWpeak_dist.txt > EG80WW_peaklabeled_dist.txt
module load bedtools2
bedtools intersect -wb -a 0423_alignSNPs_dist.txt -b EG80WW_peaklabeled_dist.txt > MOA_EG80WW_0423SNPs_dist.txt
```

And now with the file that has AF:
```
bash gene_distance.sh ZmB73v5.gff 0423_alignSNPs_AF.bed 

awk -v OFS='\t' '{print $1,$2,$3,$33,$44}' out.distance.0423_alignSNPs_AF.bed.ZmB73v5.gff.bed > 0423_alignSNPs_AF_dist.txt
awk -v OFS='\t' '{print $1,$2,$3,$14}' out.distance.EG80.WW.all_peaks.merged.inclCML333.bed.ZmB73v5.gff.bed  > EG80WWpeak_dist.txt
awk -v OFS='\t' '{print $0, "EG80WW_peak_"NR}' EG80WWpeak_dist.txt > EG80WW_peaklabeled_dist.txt
module load bedtools2
bedtools intersect -v -a 0423_alignSNPs_AF_dist.txt -b EG80WWpeaklabeled_dist.txt > nonEG80WWSNPs_AF_dist.txt
bedtools intersect -wb -a 0423_alignSNPs_AF_dist.txt -b EG80WW_peaklabeled_dist.txt > MOA_EG80WW_0423SNPs_AF_dist.txt

```


*How many peaks in each?*

```
#number of peaks in overlap with indels  (724799)
bedtools intersect -wb -a 0423_21_AlignGATKindeloverlap.out -b EG80WW_peaklabeled_dist.txt | awk -v OFS="\t" '{print $NF}' - | sort | uniq | wc -l 

#number of peaks called (1137205)
#awk -v OFS="\t" '{print $NF}' EG80WW_peaklabeled_dist.txt | sort | uniq | wc -l

#number of peaks in the 64M call (870175)
#awk -v OFS="\t" '{print $NF}' MOA_EG80WW_0423SNPs_dist.txt | sort | uniq | wc -l

#number of 64M SNPs within peaks (6638244)
#wc -l MOA_EG80WW_0423SNPs_dist.txt

#number of overlap SNPs within peaks ()
#wc -l AlignGATKoverlap_EG80WWpeaklabeled.out

#number of shared SNPs between 97M call and overlap set
#echo "Number of values of overlap in 97M call"
#bedtools intersect -a /ptmp/LAS/snodgras/VCAP_test/MOA_peak_calls/MOA_EG80WWSNPs_AF_dist.txt -b AlignGATKoverlap.out | wc -l 
#echo "Number of values of 97M in overlap call"
#bedtools intersect -a AlignGATKoverlap.out -b /ptmp/LAS/snodgras/VCAP_test/MOA_peak_calls/MOA_EG80WWSNPs_AF_dist.txt  | wc -l 
#Answer was the same both ways: 3,556,784

#number of shared peaks(747,804)
#bedtools intersect -a /ptmp/LAS/snodgras/VCAP_test/MOA_peak_calls/MOA_EG80WWSNPs_AF_dist.txt -b AlignGATKoverlap.out > test.txt
#awk -v OFS="\t" '{print $NF}' test.txt | sort | uniq | wc -l  

```

Get the GBS SNP calls for each of the NAM families
`GBS-output/populations.snps.vcf`: file with SNP calls (GBS) for all NAM lines.

```
iinit
icd /iplant/home/shared/NAM/PANDA/SVs-impute
iget -K GBS-output.tar.gz
tar -xvzf GBS-output.tar.gz
```

Get the new file of alignment SNPs for all NAM families (not just those with MOA)

```
ml singularity
singularity shell --bind $PWD /home/snodgras/irods.simg
icd /iplant/home/aseetharam/NAM_MOA/MOA_tables
iget All_SNPs_bialelic_chr1_10_new_04_2021_NAM_genotypes_cleaned_0423.bed
```
and convert it to hapmap format 
```
python reformat2hapmap.py > All_NAM_genotypes_cleaned.hmp.txt
```
Need to split by chromosome and scaffolds
```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
for i in {1..10}; do
	head -n 1 All_NAM_genotypes_cleaned.hmp.txt > projections/NAM_founders_SNPs.chr${i}.hmp.txt
	awk -v OFS='\t' -v chr=$(echo chr${i}) '$3 == chr {print $0}' All_NAM_genotypes_cleaned.hmp.txt >> projections/NAM_founders_SNPs.chr${i}.hmp.txt
done
```
There are no scaffolds in the main alignment snp file.

Need to split the GBS VCF files by chromosome as well (Creating hapmap files for each NAM population in Rafa's pipeline)
This is the script `slurm_splitgbsvcf.sh`
```
# go to project folder
cd ~/projections

# create folder to save temporary files
mkdir GBS-output/tmp

# print header to each new chromosome vcf file
# awk will print only lines that start with # and quit after first mismatch
# (this is faster than using grep, which will go through the entire file and then quit)
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
  awk '{if(/^#/)print;else exit}' GBS-output/populations.snps.vcf > GBS-output/tmp/NAM_rils_SNPs.${chr}.vcf
done

# split vcf file by chromosome
module load parallel
for i in {1..10}; do
  sem -j +0 "grep -w '^chr$i' GBS-output/populations.snps.vcf >> GBS-output/tmp/NAM_rils_SNPs.chr${i}.vcf"
done
sem --wait
grep "^scaf_" GBS-output/populations.snps.vcf >> GBS-output/tmp/NAM_rils_SNPs.scaffs.vcf
```
From Rafa's pipeline:
The bottleneck of performance now is that I have 5000+ lines in the vcf file, and transforming this file into hapmap format with TASSEL will take a lot of time. Thus I used vcftools to create a vcf file for each NAM family (i.e. 25 vcf files with ~200 RILs each). But first, I had to create a file telling which RILs belong to which family based on the information on http://maizecoop.cropsci.uiuc.edu/nam-rils.php using scripts/create_file_with_nam_rils_info.py.

```
# go to project folder
cd /work/LAS/mhufford-lab/snodgras/projections/
module load python/3.6.5-fwk5uaj #needs to be python 3.6.6 according to the pipeline

# get list of all NAM RIL names (and parents) with GBS data
head -n 1 GBS-output/populations.sumstats.tsv | cut -f 2 > nam_ril_populations.txt
# rearrange information in a table
python scripts/create_file_with_nam_rils_info.py nam_ril_populations.txt

#while I don't need to filter out any SNPs, I do need to make separate files for each chromosome and cross
# so I'll do the "filter" step in Rafa's pipeline to get the correct splitting
#and remove the `exclude-bed` argument from the commands

sbatch slurm_splitvcfbycrossfamilychr.sh 

# once the above is done, i have to sort each vcf file and export to hapmap format
cd /work/LAS/mhufford-lab/snodgras/projections/GBS-output/tmp/

# commands for sorting and transforming to hapmap
cd ...../GBS-output/tmp/
for cross in $(ls -d B73x*); do
  echo bash /work/LAS/mhufford-lab/snodgras/projections/scripts/vcf2hmp_nam-rils_snps-not-in-svs.sh $cross >> commands.txt
done

python makeSLURM.py 1 commands.txt
```
Then collapse any duplicated SNPs
```
cd /work/LAS/mhufford-lab/snodgras/projections/GBS-output/tmp/

# collapse duplicated SNPs
for cross in $(ls -d B73x*); do
  for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 ; do
  echo "/work/LAS/mhufford-lab/snodgras/projections/scripts/collapse_GBS_markers.R $cross/NAM_rils_SNPs.${cross}.${chr}.hmp.txt $cross"
  done
done > /work/LAS/mhufford-lab/snodgras/projections/scripts/commands_for_collapse-GBS-SNPs.txt

python makeSLURM.py 1 commands_for_collapse-GBS-SNPs.txt #edit to include module load r r-data-table r-doparallel and cancel out all the other module loads
```
Merge all the collapsed files together for each line
`slurm_mergecollapsedhmpfiles.sh`
```
cd /work/LAS/mhufford-lab/snodgras/projections/GBS-output/tmp/

# merge hapmap files
for cross in $(ls -d B73x*); do
  cat ${cross}/NAM_rils_SNPs.${cross}.chr1.collapsed.hmp.txt > ${cross}/NAM_rils_SNPs.${cross}.not-imputed.hmp.txt
  for chr in chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 ; do
    sed 1d ${cross}/NAM_rils_SNPs.${cross}.${chr}.collapsed.hmp.txt >> ${cross}/NAM_rils_SNPs.${cross}.not-imputed.hmp.txt
  done
done

# check number of SNPs per pop
wc -l B73*/*.not-imputed.hmp.txt
# ~2.88M (with slightly different number per population)
```
QC to see why I'm getting all exactly the same wc -l by:
```
cd /work/LAS/mhufford-lab/snodgras/projections/GBS-output/tmp/

# merge hapmap files
for cross in $(ls -d B73x*); do
  wc -l ${cross}/NAM_rils_SNPs.${cross}.chr1* 
done
```
Hopefully if they start out as different line counts, I'll be able to see where there was a bad loop potentially.

Overlay alignment data into parental GBS data
`overlay_reseq-parental-SNPs_onto_GBS-data.R` from Rafa's pipeline
The differences in the headers for HMP of the NAM founder hapmap files for the alignment is causing a lot of issues
Manually changing:
```
From this: snp.Chr.Stop    B73_allele/NAM_allele   Chr     Stop    NA      NA      NA      NA      NA      NA      NA      B73_alleleB73_allele    B97     CML103  ...
To this: rs#	alleles	Chr	pos	strand	assembly#	center	protLSID	assayLSID	panelLSID	QCcode	B73...
```
Can use sed to change:
```
for i in {1..10} ; do
sed -i 's/snp.Chr.Stop\tB73_allele\/NAM_allele\tChr\tStop\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tB73_alleleB73_allele/rs\#\talleles\tChr\tpos\tstrand\tassembly\#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\tB73/g' NAM_founders_SNPs.chr${i}.hmp.txt
done
```
And then have to remove the "chr" from before the chromsome number
```
sed -i 's/chr//g' NAM_founders_SNPs.chr1.hmp
for i in {2..10} ; do sed -i 's/chr//g' NAM_founders_SNPs.chr${i}.hmp.txt ; echo done with $i ; done
```
```
cd /work/LAS/mhufford-lab/snodgras/projections/GBS-output/tmp/
for cross in $(ls -d B73x*); do
  echo "/work/LAS/mhufford-lab/snodgras/projections/scripts/overlay_reseq-parental-SNPs_onto_GBS-data.R ${cross}/NAM_rils_SNPs.${cross}.not-imputed.hmp.txt /work/LAS/mhufford-lab/snodgras/projections/NAM_founders_SNPs.chr1.hmp.txt ${cross} ${cross}/NAM_gbs-parents_SNPs.${cross}.align-overlay.hmp.txt"
done > /work/LAS/mhufford-lab/snodgras/projections/scripts/commands_for_overlay-parental-SNPs.txt

```

Finally we need to pick the best GBS markers
use the script `select_best_SNPs_per_pop.R` 

Will require the packages: r, r-data-table, r-ggplot2, r-foreach, r-doparallel 
So make sure to update the `makeSLURM.py` to module load each of those in this order
```
ml purge
ml r-devtools r-data-table r-foreach r-doparallel
ml r-ggplot2
```
```
#check to make sure you have this directory:
/work/LAS/mhufford-lab/snodgras/projections/qc/filter_best_SNPs

cd /work/LAS/mhufford-lab/snodgras/projections/GBS-output/tmp/

for cross in $(ls -d B73x*); do
	echo "/work/LAS/mhufford-lab/snodgras/projections/scripts/select_best_SNPs_per_pop.R ${cross} ${cross}/NAM_rils_SNPs.${cross}.not-imputed.hmp.txt ${cross}/NAM_gbs-parents_SNPs.${cross}.align-overlay.hmp.txt /work/LAS/mhufford-lab/snodgras/projections/qc/filter_best_SNPs --max_missing=0.3 --window_size=15 --window_step=1 --min_snps_per_window=5" >> /work/LAS/mhufford-lab/snodgras/projections/scripts/commands_for_selectbestgbssnp.txt
done

python makeSLURM.py 1 commands_for_selectbestgbssnp.txt

mv commands_for_selectbestgbssnp_*.sub /work/LAS/mhufford-lab/snodgras/projections/GBS-output/tmp/.
```

Then we need to correct SNP names
```
cd /work/LAS/mhufford-lab/snodgras/projections/GBS-output/tmp/

for cross in $(ls -d B73x*); do
	/work/LAS/mhufford-lab/snodgras/projections/scripts/correct_SNP-names_rils.R ${cross}/NAM_gbs-parents_SNPs.${cross}.align-overlay.hmp.txt ${cross}/NAM_rils_SNPs.${cross}.not-imputed.best-markers.hmp.txt
done
```
And make sure that there's the same number of markers in these kinds of files:
Not sure this is actually needed???
```
# gbs.parents.file <- "B73xTzi8/NAM_gbs-parents_SNPs.B73xTzi8.not-in-SVs.reseq-overlay.hmp.txt"
# gbs.rils.file <- "B73xTzi8/NAM_rils_SNPs.B73xTzi8.not-in-SVs.not-imputed.best-markers.correct-marker-names.hmp.txt"

for cross in $(ls -d B73x*) ; do
	wc -l ${cross}/NAM_gbs-parents_SNPs.${cross}.align-overlay.hmp.txt
	wc -l ${cross}/NAM_rils_SNPs.${cross}.not-imputed.best-markers.correct-marker-names.hmp.txt
done
```

Before the actual projection, we need to ensure that the files are sorted:
```
cd /work/LAS/mhufford-lab/snodgras/projections/GBS-output/tmp
ml tassel

# sort parents snps
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx64g -SortGenotypeFilePlugin -inputFile ${cross}/NAM_gbs-parents_SNPs.${cross}.align-overlay.hmp.txt -outputFile ${cross}/NAM_parents_SNPs.${cross}.sorted.hmp.txt -fileType Hapmap
done

# sort rils snps
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx64g -SortGenotypeFilePlugin -inputFile ${cross}/NAM_rils_SNPs.${cross}.not-imputed.best-markers.correct-marker-names.hmp.txt -outputFile ${cross}/NAM_rils_SNPs.${cross}.not-imputed.best-markers.correct-marker-names.sorted.hmp.txt -fileType Hapmap
  run_pipeline.pl -Xmx64g -importGuess ${cross}/NAM_rils_SNPs.${cross}.not-imputed.best-markers.correct-marker-names.sorted.hmp.txt -export ${cross}/NAM_rils_SNPs.${cross}.not-imputed.best-markers.correct-marker-names.sorted.hmp.txt -exportType HapmapDiploid
done

# make sure number of rows of parental and RIL data matches in each population
for cross in $(ls -d B73x*); do
  wc -l ${cross}/NAM_parents_SNPs.${cross}.sorted.hmp.txt
  wc -l ${cross}/NAM_rils_SNPs.${cross}.not-imputed.best-markers.correct-marker-names.sorted.hmp.txt
done
```
The projection Step
`project_SNPs.sh`
```
for cross in $(ls -d B73x*); do 
	echo bash /work/LAS/mhufford-lab/snodgras/projections/scripts/project_SNPs.sh ${cross} 3000 >> /work/LAS/mhufford-lab/snodgras/projections/scripts/commands_for_projections.txt 
done

#no modules needed
python makeSLURM.py 1 commands_for_projections.txt
```
Projection QC:
```
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx64g -importGuess /work/LAS/mhufford-lab/snodgras/projections/NAM_rils_SNPs.${cross}.best-markers.projected.hmp.txt -GenotypeSummaryPlugin -endPlugin -export /work/LAS/mhufford-lab/snodgras/projections/NAM_rils_SNPs_${cross}_OverallSummary
  (echo $cross && grep "Proportion Missing" /work/LAS/mhufford-lab/snodgras/projections/NAM_rils_SNPs_${cross}_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> /work/LAS/mhufford-lab/snodgras/projections/missing_data_best-markers_after_projection.txt
done
```
Merging the projections all into one file

Let's try the `mergeGenotypeTables` plugin with tassel5
```
ml tassel
run_pipeline.pl -fork1 -h NAM_rils_SNPs.B73xB97.best-markers.projected.hmp.txt -fork2 -h NAM_rils_SNPs.B73xCML333.best-markers.projected.hmp.txt -combine3 -input1 -input2 -mergeGenotypeTables -export B97_CML333_merge.hmp.txt -runfork1 -runfork2 
```
That seemed to work so writing a loop to do it
```
ml tassel

run_pipeline.pl -fork1 -h NAM_rils_SNPs.B73xB97.best-markers.projected.hmp.txt -fork2 -h NAM_rils_SNPs.B73xCML103.best-markers.projected.hmp.txt -combine3 -input1 -input2 -mergeGenotypeTables -export NAM_rils_SNPs.ALL.best-markers.projected.hmp.txt  -runfork1 -runfork2

for cross in B73xCML228 B73xCML247 B73xCML277 B73xCML322 B73xCML333 B73xCML52 B73xCML69 B73xHp301 B73xIl14H B73xKi11 B73xKi3 B73xKy21 B73xM162W B73xM37W B73xMo18W B73xMS71 B73xNC350 B73xNC358 B73xOh43 B73xOh7B B73xP39 B73xTx303 B73xTzi8 ; do
	run_pipeline.pl -fork1 -h NAM_rils_SNPs.ALL.best-markers.projected.hmp.txt -fork2 -h NAM_rils_SNPs.${cross}.best-markers.projected.hmp.txt -combine3 -input1 -input2 -mergeGenotypeTables -export NAM_rils_SNPs.ALL.best-markers.projected.hmp.txt -runfork1 -runfork2
done
```
The loop is running out of memory, if it fails on the huge node then
```
#Split by chromosome

for cross in B73xB97 B73xCML103 B73xCML228 B73xCML247 B73xCML277 B73xCML322 B73xCML333 B73xCML52 B73xCML69 B73xHp301 B73xIl14H B73xKi11 B73xKi3 B73xKy21 B73xM162W B73xM37W B73xMo18W B73xMS71 B73xNC350 B73xNC358 B73xOh43 B73xOh7B B73xP39 B73xTx303 B73xTzi8 ; do
	for chr in {1..10} ; do
		head -n 1 NAM_rils_SNPs.${cross}.best-markers.projected.hmp.txt > NAM_rils_SNPs.${cross}.best-markers.chr-${chr}.projected.hmp.txt
		awk -v OFS='\t' -v chr=$(echo $chr) '$3 == chr {print $0}' NAM_rils_SNPs.${cross}.best-markers.projected.hmp.txt >> NAM_rils_SNPs.${cross}.best-markers.chr-${chr}.projected.hmp.txt
	done
done

ml tassel
for chr in {1..10} ; do
	run_pipeline.pl -fork1 -h NAM_rils_SNPs.B73xB97.best-markers.chr-${chr}.projected.hmp.txt -fork2 -h NAM_rils_SNPs.B73xCML103.best-markers.chr-${chr}.projected.hmp.txt -combine3 -input1 -input2 -mergeGenotypeTables -export NAM_rils_SNPs.ALL.best-markers.chr-${chr}.projected.hmp.txt  -runfork1 -runfork2
	for cross in B73xCML228 B73xCML247 B73xCML277 B73xCML322 B73xCML333 B73xCML52 B73xCML69 B73xHp301 B73xIl14H B73xKi11 B73xKi3 B73xKy21 B73xM162W B73xM37W B73xMo18W B73xMS71 B73xNC350 B73xNC358 B73xOh43 B73xOh7B B73xP39 B73xTx303 B73xTzi8 ; do
		run_pipeline.pl -fork1 -h NAM_rils_SNPs.ALL.best-markers.chr-${chr}.projected.hmp.txt -fork2 -h NAM_rils_SNPs.${cross}.best-markers.chr-${chr}.projected.hmp.txt -combine3 -input1 -input2 -mergeGenotypeTables -export NAM_rils_SNPs.ALL.best-markers.chr-${chr}.projected.hmp.txt -runfork1 -runfork2
	done
done
```
It failed on the huge node but the above worked.
To check if the columns are in the same order:
```
for i in *ALL* ; do head -n 1 $i > header_${i} ; md5sum header_${i} ; done
```
Then we can simply concatenate them:
```
head -n 1 NAM_rils_SNPs.ALL.best-markers.chr-1.projected.hmp.txt > temp_NAM_rils_SNPs.ALL.best-markers.chr-ALL.projected.hmp.txt
for chr in {1..10} ; do
	tail -n +2 NAM_rils_SNPs.ALL.best-markers.chr-${chr}.projected.hmp.txt >> temp_NAM_rils_SNPs.ALL.best-markers.chr-ALL.projected.hmp.txt
done
cut -f 1-11,63- temp_NAM_rils_SNPs.ALL.best-markers.chr-ALL.projected.hmp.txt >> NAM_rils_SNPs.ALL.best-markers.chr-ALL.projected.hmp.txt

```

###Running VCAP on new projections
Step 0: making a full kinship matrix with the base files
_Note to self, I changed the argument in_ `VCAP_mkKinship.sh`_so that the WGKM would be directed to the projected file. The original had .lix so be on the lookout for an error._

Making phenotype files for simulation:
```
paste Simulated_Phenotypes.[0-9]*.txt > Simulated_Phenotypes.txt
echo FID IID $(head -n 1 Simulated_Phenotypes.txt ) > Simulated_Phenotypes.plink
tail -n +2 Simulated_Phenotypes.txt | paste BKGD.SNPs.bed_CentIBS_K.grm.id -  >> Simulated_Phenotypes.plink
sed -i 's/ /\t/g' Simulated_Phenotypes.plink
```

### Checking the accuracy of the different SNP calls
Using the 55K SNPs on v3 from the 282 panel (panzea, `SNP55K_maize282_AGP3_20190419.hmp.txt`)
1. filter out the columns for just the NAM lines
```
#make file with nam_names.txt
B73 B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 Hp301 Il14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7B P39 Tx303 Tzi8
#find out which of the columns in the original file you'll need to cut
head -n 1 SNP55K_maize282_AGP3_20190419.hmp.txt | tr '\t' '\n' | cat -n | grep -f nam_names.txt -
#cut the columns including the hmp first 11 columns of info, there's no column for Tzi8
cut -f 1-11,51,58,79,88,90,95,101,106,111,113,141,152,158,162,165,171,172,175,186,226,230,238,241,243,264 SNP55K_maize282_AGP3_20190419.hmp.txt > SNP55K_maizeNAM_AGP3.hmp.txt
```
2. Convert from v3 to v5 coordinates
```
#1.Take my current V3 SNP coordinates and add +/-50 bp to the start and stop coordinates
sort SNP55K_maizeNAM_AGP3.hmp.txt > sorted_SNP55K_maizeNAM_AGP3.hmp.txt
##This works if it's all on the same indexing 
awk -v OFS='\t' '{print $3, $4-50, $4+50, $1, $5}' sorted_SNP55K_maizeNAM_AGP3.hmp.txt > sorted_SNP55K_maizeNAM_AGP3.bed
##to convert from 1-based to 0-based indexing
awk -v OFS='\t' '{print $3, $4-51, $4+49, $1, $5}' sorted_SNP55K_maizeNAM_AGP3.hmp.txt > sorted_SNP55K_maizeNAM_AGP3.bed

#2. Get the fasta for these coodinates on V3
module load bedtools2
tail -n +2 sorted_SNP55K_maizeNAM_AGP3.bed | bedtools getfasta -fi Zea_mays.AGPv3.30.dna.genome.fa -bed - > SNP55K_maizeNAM_AGP3.fasta

#3. Use crossmap along with the chain files to lift over the coordinates from V3 to V4 and then V4 to V5.

module load py-crossmap/0.2.7-py2-f7ajtl2
#CrossMap.py bed chain_file input.bed output.bed
CrossMap.py bed AGPv3_to_B73_RefGen_v4.chain.gz sorted_SNP55K_maizeNAM_AGP3.bed sorted_SNP55K_maizeNAM_AGP4.bed
CrossMap.py bed B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain sorted_SNP55K_maizeNAM_AGP4.bed sorted_SNP55K_maizeNAM_AGP5.bed

#4. Meanwhile, align the fasta from step 2 to the B73v5 genome and get the coordinates for best mapping (minimap or blast)
module load minimap2
minimap2 -ax sr B73.PLATINUM.pseudomolecules-v1.fasta SNP55K_maizeNAM_AGP3.fasta > SNP55K_maizeNAM_AGP3_aln2B73v5.sam
samtools view -S -b SNP55K_maizeNAM_AGP3_aln2B73v5.sam > SNP55K_maizeNAM_AGP3_aln2B73v5.bam

In the sam file there's some that have flags for reverse (16), secondary (256), or reverse and secondary(272)
The number of alignments that are only primary (0 or 16 flags) is about 200 less than what I have for the IDs in the bed file.

#### Try using blast instead + e-value filter 
module load blast-plus 
makeblastdb -in B73.PLATINUM.pseudomolecules-v1.fasta -input_type fasta -dbtype nucl -parse_seqids
blastn -max_hsps 1 -subject B73.PLATINUM.pseudomolecules-v1.fasta -query SNP55K_maizeNAM_AGP3.fasta -out blast.B73v5_SNP55K.out -outfmt "6 qseqid sseqid pident qstart qend length sstart send qcovs evalue" 

#Output looks like (header added here)
qseqid 				sseqid  pident qstart,qend,length,sstart,send qcovs evalue
1:2934015-2934115	chr1	100.000	1	100	100	3027061	3027160	100	1.88e-45
1:2934015-2934115	chr9	90.526	6	100	95	160859401	160859307	95	1.16e-27
1:3726773-3726873	chr1	100.000	1	100	100	3836189	3836288	100	1.88e-45
1:4186864-4186964	chr1	100.000	1	100	100	4248409	4248508	100	1.88e-45

So we want $3 == 100.000 $6 == 100 $9 == 100
awk -v OFS='\t' '$3 == 100.000 && $6 == 100 && $9 == 100 {print $0}' blast.B73v5_SNP55K.out > filtered.blast.B73v5_SNP55K.out

#### Try using hisat2 instead of blast
module load hisat2
hisat2-build -p 12 B73.PLATINUM.pseudomolecules-v1.fasta B73v5
hisat2 -p 12 --mp 1,1 --no-softclip -f -x B73v5 -U SNP55K_maizeNAM_AGP3.fasta  > SNP55K_maizeNAM_AGP5.sam

module load bedops samtools
sam2bed < SNP55K_maizeNAM_AGP5.sam > SNP55K_maizeNAM_AGP5.bed

#5. compare the coordinates between the output of steps 3 and 4 and keep only those which overlap
module load bedtools2
awk -v OFS='\t' '{print "chr"$1, $2, $3, $4, $5}' sorted_SNP55K_maizeNAM_AGP5.bed > chr_SNP55K_maizeNAM_AGP5.bed
bedtools intersect -a chr_SNP55K_maizeNAM_AGP5.bed -b SNP55K_maizeNAM_AGP3_aln2B73v5.bam > liftoverv5_aln_intersect.bed

Worked, but I'm not sure I like it. Will need to think on the filtering of the bam file to remove secondary alignments
Maybe if I switch a and b around? Nope not quite. I lose the rs# names
Let's try

samtools view -f 0 -f 16 SNP55K_maizeNAM_AGP3_aln2B73v5.bam  > SNP55K_primaryAlnOnly.sam
```
head -n 1 `SNP55K_primaryAlnOnly.sam`
1:28422097-28422197	272	chr9	152096520	0	2S98M	*	0	0	*	*	tp:A:S	cm:i:4	s1:i:56	NM:i:6	ms:i:136	AS:i:136	nn:i:0
To get it to bed format would be:
chr9 152096520 152096520+100
```
samtools view -f 0 -f 16 SNP55K_maizeNAM_AGP3_aln2B73v5.bam  > SNP55K_primaryAlnOnly.sam
awk -v OFS='\t' '{print $3,$4,$4+100}' SNP55K_primaryAlnOnly.sam > SNP55K_primaryAlnOnly.bed
bedtools intersect -a chr_SNP55K_maizeNAM_AGP5.bed -b SNP55K_primaryAlnOnly.bed > liftoverv5_primaryAlnOnly_intersect.bed
```
Doing this step actually lead to different lengths in the v5 coordinates from 100, so it was impossible to figure out which position was the 
real SNP. Instead I'll use the coordinates from the filtered blast results:
```
awk -v OFS='\t' '{print $2,$7+51,$8-49,$1}' filtered.blast.B73v5_SNP55K.out > filtered.blast.B73v5_SNP55K.bed
awk -v OFS='\t' '{print $2,$7+50,$8-48,$1}' filtered.blast.B73v5_SNP55K.out > filtered.blast.B73v5_SNP55K.bed
```
Also difficult to get the one position 
```
#This seems to work from hisat2, but not sure if I now have messed up the getfasta part...
awk -v OFS='\t' '{print $1,$2+51,$3-49,$4}' SNP55K_maizeNAM_AGP5.bed > SNPpos_SNP55K_maizeNAM_AGP5.bed
```
To double check indexing issues:
```
#1 grab examples from our SNP list
##to convert from 1-based to 0-based indexing
##to convert from 1-based to 0-based indexing
awk -v OFS='\t' '{print $3, $4, $4, $1, $5}' sorted_SNP55K_maizeNAM_AGP3.hmp.txt | head > exampleTRUTH_SNP55K_maizeNAM_AGP3.bed
awk -v OFS='\t' '{print $3, $4-51, $4+49, $1, $5}' sorted_SNP55K_maizeNAM_AGP3.hmp.txt | head > example_SNP55K_maizeNAM_AGP3.bed

#2 get the fasta for those examples on b73v3
module load bedtools2
tail -n +2 example_SNP55K_maizeNAM_AGP3.bed | bedtools getfasta -fi Zea_mays.AGPv3.30.dna.genome.fa -bed - > example_SNP55K_maizeNAM_AGP3.fasta

#3 run the alignments on B73v3
module load hisat2
hisat2-build -p 12 Zea_mays.AGPv3.30.dna.genome.fa B73v3
hisat2 -p 12 --mp 1,1 --no-softclip -f -x B73v3 -U example_SNP55K_maizeNAM_AGP3.fasta  > example_SNP55K_maizeNAM_AGP5.sam

module load bedops samtools
sam2bed < example_SNP55K_maizeNAM_AGP5.sam > example_SNP55K_maizeNAM_AGP5.bed

#4 try the formatting to see if we get the same position
awk -v OFS='\t' '{print $1,$2+51,$3-49,$4}' example_SNP55K_maizeNAM_AGP5.bed #this gives the right answers
```
To connect the `SNPpos` to the genotypes in the hapmap file
```
awk -v OFS='\t' '{print $4,$1":"$2"-"$3}' sorted_SNP55K_maizeNAM_AGP3.bed > SNP55K.100bptoRS.key.txt 

sort -k4,4 SNPpos_SNP55K_maizeNAM_AGP5.bed | awk '!seen[$4]++' - > sorted_SNPpos_SNP55K_maizeNAM_AGP5.bed
sort -k2,2 SNP55K.100bptoRS.key.txt  > sorted_SNP55K.100bptoRS.key.txt 

join -1 2 -2 4 sorted_SNP55K.100bptoRS.key.txt sorted_SNPpos_SNP55K_maizeNAM_AGP5.bed | sort -k 2,2 | uniq -u - > RSjoined_SNPpos_SNP55K_maizeNAM_AGP5.bed
```
Now to join with the hapmap file and get the genotypes on v5 coodinates
```
head -n 1 sorted_SNP55K_maizeNAM_AGP3.hmp.txt > SNP55K_maizeNAM_AGP5.joined
sort -k2,2 RSjoined_SNPpos_SNP55K_maizeNAM_AGP5.bed > sorted_RSjoined_SNPpos_SNP55K_maizeNAM_AGP5.bed
tail -n +2 sorted_SNP55K_maizeNAM_AGP3.hmp.txt | sort -k1,1 | join -1 1 -2 2 - sorted_RSjoined_SNPpos_SNP55K_maizeNAM_AGP5.bed >> SNP55K_maizeNAM_AGP5.joined
#add col names for the joined columns manually

sed -i 's/ /\t/g' SNP55K_maizeNAM_AGP5.joined
#26 lines that have "wrong" in them and an extra field from the original hmp.txt; removing them while formatting (should have 50881 lines)
grep -v wrong SNP55K_maizeNAM_AGP5.joined | awk -v OFS='\t' '{print $1,$2,$38,$39,$5,"AGPv5",$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36}' - > SNP55K_maizeNAM_AGP5.hmp.txt
#manually change col name from "AGPv5" to "assembly#"
```

To get the founder genotypes for the reseq-snps
```
cp /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/l-gatk-illumina/c-gatk/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf .
#Use a huge node on nova
ml tassel
run_pipeline.pl -Xmx10g -importGuess B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf -export B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.hmp.txt -exportType HapmapDiploid
```

The matching script is taking too long on the matching loop, probably because the gatk/alignment file is so big. 
I'll filter it first with `bedtools intersect` so that the input files are smaller

```
#simplemerge.py from Arun
import pandas as pd
fileA = pd.read_csv("file1.txt", sep='\t', header=None, names=['FileA.1'])
fileB = pd.read_csv("file2.txt", sep='\t', header=0)
​
print(fileA.head())
print(fileB.head())
​
mergedFile = fileA.merge(fileB, left_on='FileA.1', right_on='rs')
print(mergedFile.head())
mergedFile.to_csv("merged_FileA-and-B.tsv", sep='\t', index=False)
```

```
module load bedtools2  py-pandas
awk -v OFS="\t" '{print $3,$4,$4+1,$1}' SNP55K_maizeNAM_AGP5.hmp.txt | tail -n +2 > SNP55K_maizeNAM_AGP5.coords.bed
awk -v OFS="\t" '{print "chr"$3,$4,$4+1,$1}' B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.hmp.txt | tail -n +2 > gatk.coords.bed
bedtools intersect -wa -a gatk.coords.bed -b SNP55K_maizeNAM_AGP5.coords.bed > filtered_gatk.coords.bed
cut -f 4 filtered_gatk.coords.bed > filtered_gatk_rsID.txt
python simplemerge.gatk.py #fileA = filtered_gatk_rsID.txt fileB = B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.hmp.txt

awk -v OFS="\t" '{print $3,$4,$4+1,$1}' ../All_NAM_genotypes_cleaned.hmp.txt |tail -n +2 > align.coords.bed
bedtools intersect -wa -a align.coords.bed -b SNP55K_maizeNAM_AGP5.coords.bed > filtered_align.coords.bed
cut -f 4 filtered_align.coords.bed > filtered_align_rsID.txt
python simplemerge.align.py #fileA = filtered_align_rsID.txt fileB = ../All_NAM_genotypes_cleaned.hmp.txt

#rm *coords.bed
```
Use the `matching_genotypes.R` script to match genotypes between files (make sure it's executable)
R won't read in the `#` in the hmp header names, so manually changing them
*change gatk name from ms37w to m37w*
```
module load r r-tidyverse
 ./matching_genotypes.R SNP55K_maizeNAM_AGP5.hmp.txt filtered_gatkSNPs.tsv gatk_matchingOutput.csv
./matching_genotypes.R SNP55K_maizeNAM_AGP5.hmp.txt filtered_alignSNPs.tsv align_matchingOutput.csv
```
