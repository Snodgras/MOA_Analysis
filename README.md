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

#VCAP of 21 lines with different peaks called for MOA
###Input files
From Julia:
- `EG80.WW.all_peaks.merged.21_NAMs.bed` which is a bed file with peak coordinates (unlabeled)
- `Genotypes_21_lines_AF_for_Sam.bed` which is the SNP coordinates and genotypes across all lines + allele frequency

##Step 0:
I need to create: *start at line 364*
_Create allele frequency and distance containing SNP bed files for MOA peak and non-peak regions_
_If using new markers, make the full kinship matrix_

Pulling out coordinates, AF 
```
#select the coordinate columns and allele frequency column
awk -v OFS='\t' '{print $1,$2,$3,$29}' Genotypes_21_lines_AF_for_Sam.bed > AllSNPs_AF.bed
```
Calculating distance from each peak to nearest gene (requires gff of B73v5)
Use script `gene_distance.sh` and `slurm_gene_distance.sh`
```
#!/bin/bash
GFF=$1
MOABED=$2
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load bedops
#convert2bed --input=gff --output=bed < $GFF > gff.bed

#need to separate out only the gene coordinates (not exon, chromosome, mRNA, etc.)
#awk -v OFS='\t' '$8 ~ /gene/ {print $0}' gff.bed | sort -k1,1 -k2,2n > sorted.gene.bed

#sort the peak bedfile
sort -k1,1 -k2,2n $MOABED > sorted.$MOABED

module load bedtools2
bedtools closest -D ref -t first -a sorted.$MOABED -b sorted.gene.bed > out.distance.${2}.${1%_*}.bed

#slurm command

bash gene_distance.sh ZmB73v5.gff AllSNPs_AF.bed #(manually remove the header line first)
#bash gene_distance.sh ZmB73v5.gff EG80.WW.all_peaks.merged.21_NAMs.bed #add distance to all the peaks, can probably get this from the next intersect step
```
Then to make our alignment and distance bed files needed for the SNP sampling script
```
awk -v OFS='\t' '{print $1,$2,$3, $4,$15}' out.distance.AllSNPs_AF.bed.ZmB73v5.gff.bed > alignSNP_AF_dist.txt
awk -v OFS='\t' '{print $1,$2,$3,$14}' out.distance.EG80.WW.all_peaks.merged.21_NAMs.bed.ZmB73v5.gff.bed  > EG80WWpeak_dist.txt
module load bedtools2
bedtools intersect -v -a alignSNP_AF_dist.txt -b EG80WWpeak_dist.txt > nonEG80WWSNPs_AF_dist.txt
```
and also to make a MOA SNP AF Dist file
```
awk -v OFS='\t' '{print $0, "EG80WW_peak_"NR}' EG80WWpeak_dist.txt > EG80WW_peaklabeled_dist.txt
bedtools intersect -wb -a alignSNP_AF_dist.txt -b EG80WW_peaklabeled_dist.txt > MOA_EG80WWSNPs_AF_dist.txt
```

The `nonEG80WWSNPs_AF_dist.txt` has the format:
chromosome \t start \t stop(actual SNP location) \t AF \t distance from nearest gene

### 1 : Run background selection (permutation starts)
_Select a matched background_
	- run `SNPselection.R`
	- make sure the resulting files are being written to specific directories

commands to run `SNPselection.R`
```
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load r r-devtools r-tidyr r-dplyr
./SNPselection.R $MOA_AF_dist.txt $nonMOASNPs_AF_dist.txt $outputdirectory/filename
```
Use `makeSLURM.py` to be able to run it in parallel

```
#make a command file
for i in {0..99} ; do 
	mkdir 1SNPperPeak/EG80WW_perm_$i
	echo module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/ \; module load r r-devtools r-tidyr r-dplyr \; ./SNPselection.R MOA_EG80WWSNPs_AF_dist.txt  nonEG80WWSNPs_AF_dist.txt 1SNPperPeak/EG80WW_perm_${i}/coords_${i}_ >> 1_command.txt ;
done
python makeSLURM.py 1 1_command.txt

```
then to submit the slurm scripts
```
sbatch 1_commands_0.sub #to make sure it works ok
for i in {1..99} ; do sbatch --dependency=afterok:[job-id] 1_command_${i}.sub ; done
```
There are 1220878 peak coordinates in the starting bed file.
1 SNP per Peak sampling results in 883684 SNP region coordinates
Is the discrepancy caused by no SNPs called in some peak regions? 

There are slightly more SNPs in the new 21 line MOA selection file (883684) than in the one original MOA SNP selection file I have from the first attempt (870175)
But there are 59,922,979 SNPs total... seems strange that some peaks wouldn't have SNPs?

```
cut -f 10 1SNPperPeak/EG80WW_perm_0/coords_0_MOA_distrib.SNPs.txt | sort - | uniq #to make sure that no peaks are repeated, see which might be missing
#try bedtools intersect -wa (instead of -wb) to maybe find those peaks without a SNP? 
```
###2 : Make the Kinship Matrices
_Change the whole genome kinship name and thus which genotypes are used in the slurm script for making the kinships_
Using the gatk resequenced SNPs (`NAM_rils_SNPs_allchr.hmp.txt`)

*For the first time making the whole genome kinship matrix*
specify exclusive use of a huge node and increase the heap size in the tassel argument to 
2/3 the total memory  of the huge node (2Tb) `-Xmx2000g` instead of `-Xmx64g`

_Add the slurm script to each permutation directory_
Before doing this, make sure that the name of the whole genome kinship matrix is correctly specified

Also change all of the files from `coords_*_BKGD_distrib.SNPs.txt` to `BKGD_distrib.SNPs.txt`
Same for MOA
```
for i in {1..99} ; do
	cd EG80WW_perm_${i}
	mv coords_${i}_BKGD_distrib.SNPs.txt BKGD_distrib.SNPs.txt
	mv coords_${i}_MOA_distrib.SNPs.txt MOA_distrib.SNPs.txt
	cd ..
done

for i in {1..99} ; do cp ../slurm_VCAP_mkKinship.sh EG80WW_perm_${i}/. ;done
```
_Submit the slurm scripts_
```
for i in {2..99} ; do
	cd EG80WW_perm_${i}
	sbatch --dependency=afterok:4435301 slurm_VCAP_mkKinship.sh
	cd ..
done
```

To see if there were any kinships that need re-running:
```
for i in {1..99} ; do echo ${i} ; ls EG80WW_perm_${i} | wc -l ; done #there should be 15 files
```
Needed to rerun perm 71

###3: running LDAK
_Get the phenotype data_

Phenotype data from: Merritt Burch `all_NAM_phenos.csv`

_cleaning up datasets to be converted to PLINK format_
`http://dougspeed.com/phenotypes-and-covariates/`
```
awk -v FS=, -v OFS='\t' -v b=9 '{for (i=b; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}' all_NAM_phenos.csv > temp.txt
awk -v FS=, -v OFS='\t' '{print $1,$1}' all_NAM_phenos.csv | paste - temp.txt > all_NAM_phenos.plink
rm temp.txt
sed -i "s/\"//g" all_NAM_phenos.plink
```

_test on 1 permutation_

```
#in the slurm script
bash /work/LAS/mhufford-lab/snodgras/MOA/VCAP_21_lines/VCAP_runLDAK.sh 1SNPperPeak_fullResults kinship.list /work/LAS/mhufford-lab/snodgras/MOA/VCAP_21_lines/all_NAM_phenos.plink
```
Got this error: 
```
Reading responses for 5187 samples from /work/LAS/mhufford-lab/snodgras/MOA/VCAP_21_lines/all_NAM_phenos.plink
Error, Z003E0046___Z003E0046 appears twice in /work/LAS/mhufford-lab/snodgras/MOA/VCAP_21_lines/all_NAM_phenos.plink
```
They appear to have identical values. So I'll manually remove one of them. 
Had the same issue for `Z006E0118___Z006E0118`
To systematically remove them, I'll use unique
```
uniq all_NAM_phenos.plink > all_NAM_phenos_uniq.plink
rm all_NAM_phenos.plink
mv all_NAM_phenos_uniq.plink all_NAM_phenos.plink
```

Making sure that the time request is appropriate before copying slurm script to the other perm directories

```
for i in {1..99} ; do cp EG80WW_perm_0/slurm_VCAP_runLDAK.sh EG80WW_perm_${i}/. ; done
for i in {1..99} ; do cd EG80WW_perm_${i} ; sbatch slurm_VCAP_runLDAK.sh ; cd .. ; done
```

###4: Formatting the results
_Make sure the final file is written to the general directory (otherwise will be deleted)_
_Will delete all the rest of the files_
	- run `format_hests.sh`
	
###5: QC
The results weren't exactly the same as when we ran it with 17 genotypes, so checking: 
1. Play around with the LDAK arguments
* Double check model convergences
Check number of traits that converged
Check which traits converged (or not) across perms

```
echo perm trait convergence_status > convergence.results
for i in {0..99} ; do 
	for j in {1..143} ; do
		if grep -q "Converged YES" EG80WW_perm_${i}/1SNPperPeak_fullResults.${j}.reml
		then
			echo $i $j "YES" >> convergence.results
		else
			echo $i $j "NO" >> convergence.results
		fi
	done
done
sed -i 's/ /\t/g' convergence.results
grep -c "NO" convergence.results
```
Out of 14300 traits x permutations, only 2409 failed to converge

```
grep "NO" convergence.results | cut -f 2 | sort | uniq -c
```
89 traits failed to converge at least once
(I took the output and put it into an excel spreadsheet for easier data sorting. I took the list of trait numbers and made a file `noConverge.trait.list` to make the next step easier)

```
while read -r line ; 
do
	let n=${line}+2
	head -n 1 ../all_NAM_phenos.plink | cut -f $n
done < noConverge.trait.list 
```
(I copied the output into the excel spreadsheet that was ordered in the same way for the trait number as the input list for while loop)
Looks like metabolite traits were the most likely to fail in the majority (or all) of the permutations
Other traits that failed < 10 times were things related to pentrometer, DTA, disease resistance, other metabolite traits, 
and a few plant architecture traits

```
grep "NO" convergence.results | cut -f 1 | sort | uniq -c
```

~15-30 traits in any given permutation would fail to converge, so it's not like one set of regions failed a lot more than others

* alter one of the arguments (constrain = TRUE) _was already used in the run_
2. Check if certain genotypes were identified as notable in the original GWAS which might explain differences with the new lines
3. Make side by side, more polished graphs for easier comparison of results from the 17 genotype run and the 21 genotype run.

##bQTL VCAPs
###0: format files

starting files: 
`WW-MOA_ALLmm_merged.FDR0_05.tsv`

####Plink clumping:
Want to clump the bQTL SNPs together by LD
1. Get the `.hmp.txt.gz` file into plink format (to use in `--file`)
First, the `NAM_rills_SNPs_allchr.hmp.txt.gz` was split into a hmp file for each chromosome
Going to try converting to VCF first (so we can avoid the non-biallelic errors)
```#slurm_makePlinkFile.sh
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=4:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue


module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load tassel

#run_pipeline.pl -Xmx64g -fork1 -importGuess NAM_rils_SNPs_chr10.hmp.txt -export -exportType VCF 

run_pipeline.pl -Xmx64g -fork1 -importGuess NAM_rils_SNPs_allchr.hmp.txt.gz -export -exportType VCF
```
If the above works for the full file, then we don't have to worry about the merge.

```#slurm_convertVCF.sh
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=4:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue


module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load plink 
plink --vcf NAM_rils_SNPs_allchr.vcf --make-bed --out NAM_rils_SNPs_allchr

```

Then use `plink --vcf reference.vcf --biallelic-only --out reference` to get rid of non-biallelic alleles
Then go ahead with the merge

Then the resulting out files (ped and map) needed to be merged together:
create a file `allfiles.txt` with the name of each ped/map file duo on a line (9 lines total)
```#slurm_mergePlinkFiles.sh
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=4:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load plink

plink --file NAM_rils_SNPs_chr10.plk --merge-list allfiles.txt --make-bed --out NAM_rils_SNPs_allchr
```
2. Turn the `.tsv` into a bed-like format??? I need to change the first column to split into the coordinates
```
echo "chr\tstart\tstop\tSNP\tP" |sed 's/\\t/\t/g' - > bQTL.padj_probgeno.bed
awk -v OFS='\t' '{split($1,a,"-"); split(a[2],b,"_"); print b[1],b[2],b[2]+1,$2,$13}' WW-MOA_ALLmm_merged.FDR0_05.tsv >> bQTL.padj_probgeno.bed
```
3. Run the `--clump` command on the plink data (check settings)
```
#slurm_clumpSNPs.sh
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue

module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load plink

plink --file NAM_rils_SNPs_allchr --clump bQTL.padj_probgeno.bed --clump-kb 0.150 --clump-r2 0.99
```

The above didn't work... `Warning: No significant --clump results.  Skipping.`

I'm going to try it again but add the variables to specify p-value column
`--clump-field P`

Same issue, will try again using 10% and 30% cutoffs for the "p-values" since they're really FDRs
`--clump-p1 0.1 --clump-p2 0.3`
and add a specific output file character string: `-out NAM_rils_SNPs_allchr.clumped`

Doesn't work because the "P" field has some values in scientific notation
First need to convert to completely numeric, then run clump:
`awk '{printf "%s\t%f\n", $0, $5}' bQTL.padj_probgeno.bed | cut -f 1-4,6 > bQTL.padj_probgeno.numeric.bed`
manually change the first row (header) from 0.00000 to P for the last column

Also making the SNP marker names the same: `sed -i -e 's/B73-chr/S/' bQTL.padj_probgeno.numeric.bed`

Retrying with all of the specific flags of the previous clump attempts (and it worked!)

4. Turn clumped output into a bed-ish file to run the intersect step with of the original pipeline to create bQTL and non-bQTL SNPs

Reformat the clumped index SNPs in bed file format with +/- 75 bp on either side:
```
tail -n +2 NAM_rils_SNPs.padj_probgeno.clumped.clumped | awk -v OFS='\t' '{print "chr"$1, $4-75, $4+75, $3}' - > NAM_rils_SNPs.padj_probgeno.clumped.bed
```
manually remove the two empty rows at the bottom of the file

Get allele frequency/distance bed files for bQTL SNPs and nonbQTL SNPs
Can start with the intersection step since the genotype file has already been formatted. 
```
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load bedtools2
bedtools intersect -v -a alignSNP_AF_dist.txt -b NAM_rils_SNPs.padj_probgeno.clumped.bed > nonbQTL_AF_dist.txt
bedtools intersect -wb -a alignSNP_AF_dist.txt -b NAM_rils_SNPs.padj_probgeno.clumped.bed > bQTL_AF_dist.txt
```
5. Run the SNP selection script
Modified to `bQTL.SNPselection.R` to not subsample the bQTL, it should only create a matched background

Need to create commands to run it in parallel
```
for i in {0..99} ; do mkdir bQTL/perm_${i} ; echo module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/ \; module load r r-devtools r-tidyr r-dplyr \; ./bQTL.SNPselection.R  bQTL_AF_dist.txt nonbQTL_AF_dist.txt bQTL/perm_${i}/coords_${i}_ >> bQTL.selection.commands ; done
python ../makeSLURM.py 1 bQTL.selection.commands
```

6. Making the Kinships
_Change the whole genome kinship name and thus which genotypes are used in the slurm script for making the kinships_
Using the gatk resequenced SNPs (`NAM_rils_SNPs_allchr.hmp.txt`)

*For the first time making the whole genome kinship matrix*
specify exclusive use of a huge node and increase the heap size in the tassel argument to 
2/3 the total memory  of the huge node (2Tb) `-Xmx2000g` instead of `-Xmx64g`

_Add the slurm script to each permutation directory_
Before doing this, make sure that the name of the whole genome kinship matrix is correctly specified

Also change all of the files from `coords_*_BKGD_distrib.SNPs.txt` to `BKGD_distrib.SNPs.txt`
Same for bQTL
```
for i in {0..99} ; do
	cd perm_${i}
	mv coords_${i}_BKGD_distrib.SNPs.txt BKGD_distrib.SNPs.txt
	mv coords_${i}_bQTL_distrib.SNPs.txt bQTL_distrib.SNPs.txt
	cd ..
done

#make sure to change the lines of "MOA" to "bQTL" before copying into other directories
for i in {1..99} ; do cp ../slurm_VCAP_mkKinship.sh perm_${i}/. ;done
```
_Submit the slurm scripts_

```
cd perm_0 
sbatch slurm_VCAP_mkKinship.sh
cd ..
for i in {1..99} ; do
	cd perm_${i}
	sbatch --dependency=afterok:4554670 slurm_VCAP_mkKinship.sh
	cd ..
done
```

To see if there were any kinships that need re-running:

```
for i in {0..99} ; do echo ${i} ; ls perm_${i} | wc -l ; done #there should be 15 files
```

None need to be rerun

7. Running LDAK

Use this slurm script for `slurm_VCAP_runLDAK.sh`
```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=10:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN #Can be edited out if submitting many jobs
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --exclusive
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
if [[ ! -f kinship.list ]] ; then 
	echo bQTL.SNPs.bed_CentIBS_K > kinship.list
	echo BKGD.SNPs.bed_CentIBS_K >> kinship.list
	echo CentIBS_Rest_K >> kinship.list
fi


bash /work/LAS/mhufford-lab/snodgras/MOA/VCAP_21_lines/VCAP_runLDAK.sh bQTL_Full kinship.list /work/LAS/mhufford-lab/snodgras/MOA/VCAP_21_lines/all_NAM_phenos.plink  
```

```
for i in {1..99} ; do cp perm_0/slurm_VCAP_runLDAK.sh perm_${i}/. ; done
for i in {1..99} ; do cd perm_${i} ; sbatch --dependency=afterok:4555071 slurm_VCAP_runLDAK.sh ; cd .. ; done
```

8. Getting the results formatted correctly
_Make sure the final file is written to the general directory (otherwise will be deleted)_
_Will delete all the rest of the files_
	- run `format_hests.sh`

9. Convergence QC:
Check which traits converged (or not) across perms

```
echo perm trait convergence_status > convergence.results
for i in {0..99} ; do 
	for j in {1..143} ; do
		if grep -q "Converged YES" perm_${i}/bQTL_Full.${j}.reml
		then
			echo $i $j "YES" >> convergence.results
		else
			echo $i $j "NO" >> convergence.results
		fi
	done
done
sed -i 's/ /\t/g' convergence.results
grep -c "NO" convergence.results
```
Out of 14300 traits x permutations, only 1488 failed to converge compared to 2409 using the 1 SNP per peak method.

```
grep "NO" convergence.results | cut -f 2 | sort | uniq -c
```
37 traits failed to converge at least once (compared to 89)
(I took the output and put it into an excel spreadsheet for easier data sorting. I took the list of trait numbers and made a file `noConverge.trait.list` to make the next step easier)

```
while read -r line ; 
do
	let n=${line}+2
	head -n 1 ../all_NAM_phenos.plink | cut -f $n
done < noConverge.trait.list 
```
(I copied the output into the excel spreadsheet that was ordered in the same way for the trait number as the input list for while loop)
Looks like vitamin E traits, tassel traits, and kernel traits like protein and totalk weight were the most likely to fail in the majority (or all) of the permutations
Other traits that failed < 10 times were things related to pentrometer, ASI, other metabolite traits, 
and a few plant architecture traits

```
grep "NO" convergence.results | cut -f 1 | sort | uniq -c
```

~10-20 traits in any given permutation would fail to converge, so it's not like one set of regions failed a lot more than others

### Re-doing the bQTL to be just the bQTL, not a 75+/- bp window around them 
Starting at step 4 from the clumping
4. Turn clumped output into a bed-ish file to run the intersect step with of the original pipeline to create bQTL and non-bQTL SNPs

```
tail -n +2 NAM_rils_SNPs.padj_probgeno.clumped.clumped | awk -v OFS='\t' '{print "chr"$1, $4-1, $4, $3}' - > NAM_rils_SNPs.padj_probgeno.clumped.strict.bed
```
manually remove the two empty rows at the bottom of the file

Get allele frequency/distance bed files for bQTL SNPs and nonbQTL SNPs
Can start with the intersection step since the genotype file has already been formatted. 
```
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load bedtools2
bedtools intersect -v -a alignSNP_AF_dist.txt -b NAM_rils_SNPs.padj_probgeno.clumped.strict.bed > nonbQTL_AF_dist.strict.txt
bedtools intersect -wb -a alignSNP_AF_dist.txt -b NAM_rils_SNPs.padj_probgeno.clumped.strict.bed > bQTL_AF_dist.strict.txt
```
Then double check it worked correctly:
```
wc -l nonbQTL_AF_dist.strict.txt #59822519
wc -l bQTL_AF_dist.strict.txt #100460
```

5. Run the SNP selection script
Modified to `bQTL.SNPselection.R` to not subsample the bQTL, it should only create a matched background

Need to create commands to run it in parallel
```
for i in {0..99} ; do mkdir bQTL_strict/perm_${i} ; echo module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/ \; module load r r-devtools r-tidyr r-dplyr \; ./bQTL.SNPselection.R  bQTL_AF_dist.strict.txt nonbQTL_AF_dist.strict.txt bQTL_strict/perm_${i}/coords_ >> bQTL_strict.selection.commands ; done
python ../makeSLURM.py 1 bQTL_strict.selection.commands
```
Double check that the SNP selection worked as intended

```
wc -l bQTL_strict/perm_0/coords_bQTL_distrib.SNPs.txt #100461
wc -l bQTL_strict/perm_0/coords_BKGD_distrib.SNPs.txt #100461
for i in {0..9} ; do md5sum bQTL_strict/perm_${i}/coords* ; done #should be the same for all bQTL, different for all BKGD
```

6. Making the Kinships

_Add the slurm script to each permutation directory_
Before doing this, make sure that the name of the whole genome kinship matrix is correctly specified

Also change all of the files from `coords_*_BKGD_distrib.SNPs.txt` to `BKGD_distrib.SNPs.txt`
Same for bQTL
```
for i in {0..99} ; do
	cd perm_${i}
	mv coords_BKGD_distrib.SNPs.txt BKGD_distrib.SNPs.txt
	mv coords_bQTL_distrib.SNPs.txt bQTL_distrib.SNPs.txt
	cd ..
done

#make sure to change the lines of "MOA" to "bQTL" before copying into other directories
for i in {1..99} ; do cp perm_0/slurm_VCAP_mkKinship.sh perm_${i}/. ;done
```
_Submit the slurm scripts_

```
cd perm_0 
sbatch slurm_VCAP_mkKinship.sh
cd ..
for i in {1..99} ; do
	cd perm_${i}
	sbatch --dependency=afterok:4564545 slurm_VCAP_mkKinship.sh
	cd ..
done
```

To see if there were any kinships that need re-running:

```
for i in {0..99} ; do echo ${i} ; ls perm_${i} | wc -l ; done #there should be 15 files
```
None needed to be rerun

7. Running LDAK

Use this slurm script for `slurm_VCAP_runLDAK.sh`
```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=10:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --exclusive
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
if [[ ! -f kinship.list ]] ; then 
	echo bQTL.SNPs.bed_CentIBS_K > kinship.list
	echo BKGD.SNPs.bed_CentIBS_K >> kinship.list
	echo CentIBS_Rest_K >> kinship.list
fi


bash /work/LAS/mhufford-lab/snodgras/MOA/VCAP_21_lines/VCAP_runLDAK.sh bQTL_strict kinship.list /work/LAS/mhufford-lab/snodgras/MOA/VCAP_21_lines/all_NAM_phenos.plink  
```

```
for i in {1..99} ; do cp perm_0/slurm_VCAP_runLDAK.sh perm_${i}/. ; done
for i in {1..99} ; do cd perm_${i} ; sbatch --dependency=afterok:4564690 slurm_VCAP_runLDAK.sh ; cd .. ; done
```

8. Getting the results formatted correctly
_Make sure the final file is written to the general directory (otherwise will be deleted)_
_Will delete all the rest of the files_
	- run `format_hests.sh`

9. Convergence QC:
Check which traits converged (or not) across perms


```
echo perm trait convergence_status > convergence.results
for i in {0..99} ; do 
	for j in {1..143} ; do
		if grep -q "Converged YES" perm_${i}/bQTL_strict.${j}.reml
		then
			echo $i $j "YES" >> convergence.results
		else
			echo $i $j "NO" >> convergence.results
		fi
	done
done
sed -i 's/ /\t/g' convergence.results
grep -c "NO" convergence.results
```
Out of 14300 traits x permutations, 866 failed to converge compared to 2409 using the 1 SNP per peak method.

```
grep "NO" convergence.results | cut -f 2 | sort | uniq -c
```
48 traits failed to converge at least once (compared to 89)
(I took the output and put it into an excel spreadsheet for easier data sorting. I took the list of trait numbers and made a file `noConverge.trait.list` to make the next step easier)

```
while read -r line ; 
do
	let n=${line}+2
	head -n 1 ../all_NAM_phenos.plink | cut -f $n
done < noConverge.trait.list 
```
(I copied the output into the excel spreadsheet that was ordered in the same way for the trait number as the input list for while loop)
The four traits that failed >50 times were ear row number, cob length, and vitamin E 
Most traits that failed, failed in only a few permutations (< 10). Things like vitamin E, kernel traits, brace roots, penetrometer, metabolites, leaf wax

```
grep "NO" convergence.results | cut -f 1 | sort | uniq -c
```

~5-10 traits in any given permutation would fail to converge, so it's not like one set of regions failed a lot more than others

###Zipping the directories of the VCAP runs for archival purposes
Using the `tar.gz` for zipping

```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=10:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#tar -zcvf 1SNPperPeak.tar.gz 1SNPperPeak/
tar -zcvf bQTL_strict.tar.gz bQTL_strict/
```
Check the size to decide how best to share the zipped directories
`20G Jul 19 15:27 1SNPperPeak.tar.gz`

##Simulated Phenotypes
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

##All figures made with `VCAP_Figures.Rmd`