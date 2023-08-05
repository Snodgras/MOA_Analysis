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

<<<<<<< HEAD
##All figures made with `VCAP_Figures.Rmd`
=======
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
>>>>>>> origin/master
