#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue

####arguments####
WGKM=$1 #name of the whole genome kinship matrix (NAMRILs_projectedSNPs_CentIBS_wholegenomekinship.txt.gz)
MOAK=$2 #path/to/myBedFile1.bed
BKGD_K=$3

module load tassel

if [[ ! -f "$WGKM" ]] ; then
	echo "$TIMESTAMP Beginning full kinship matrix estimation" > VCAPtest.log ;
	run_pipeline.pl \
		-Xmx64g \
		-importGuess NAM-ril-genotypes/NAM_rils_SNPs_allchr.hmp.txt.gz.lix  \
		-KinshipPlugin \
		-method Centered_IBS \
		-endPlugin \
		-export $WGKM \
		-exportType SqrMatrix
	echo "$TIMESTAMP Finished full kinship matrix." >> VCAPtest.log ;
fi

#This works for 1 kinship matrix
#run_pipeline.pl -Xmx64g \
#        -fork1 -importGuess NAM-ril-genotypes/NAM_rils_SNPs_allchr.hmp.txt.gz.lix \
#        -fork2 -FilterSiteBuilderPlugin -bedFile MOA_peak_calls/Peaks/EG80/WW/EG80.WW.all_peaks.merged.inclCML333.bed -endPlugin \
#        -input1 -KinshipPlugin -method Centered_IBS -endPlugin \
#        -fork3 -export EG80WWall_CentIBS_kinship1 -exportType SqrMatrixBin -input2 \
#        -fork4 -SubtractDistanceMatrixPlugin -wholeMatrix NAMRILs_projectedSNPs_CentIBS_wholegenomekinship.txt.gz -endPlugin \
#        -input2 -export EG80WWall_CentIBS_kinshipRest -exportType SqrMatrixBin
#echo "$TIMESTAMP Finished kinship filtering and subtraction. Beginning LDAK" >> VCAPtest.log

#This works for 2 kinship matrices
run_pipeline.pl -Xmx64g \
	-importGuess path/to/myGenotypes.hmp.txt.gz.lix \
	-FilterSiteBuilderPlugin -bedFile ${MOAK}  -endPlugin \
	-KinshipPlugin -method Centered_IBS -endPlugin \
	-export ${MOAK}_CentIBS_K -exportType SqrMatrix 
run_pipeline.pl -Xmx64g \
	-importGuess path/to/myGenotypes.hmp.txt.gz.lix \
	-FilterSiteBuilderPlugin -bedFile ${BKGDK}  -endPlugin \
	-KinshipPlugin -method Centered_IBS -endPlugin \
	-export ${BKGDK}_CentIBS_K -exportType SqrMatrix 
run_pipeline.pl -Xmx64g \
	-fork1 -importGuess ${MOAK}_CentIBS_K.txt \
	-fork2 -importGuess ${BKGDK}_CentIBS_K.txt \
	-combine3 -input1 -input2 \
	-SubtractDistanceMatrixPlugin -wholeMatrix $WGKM -endPlugin \
	-export ${WGKM%../../}_Rest_K -exportType SqrMatrixBin
###

run_pipeline.pl -Xmx64g \
	-fork1 -importGuess ../../NAM-ril-genotypes/NAM_rils_SNPs_allchr.hmp.txt.gz.lix \
	-fork2 -FilterSiteBuilderPlugin -bedFile ${MOAK} -endPlugin \
	-input1 -KinshipPlugin -method Centered_IBS -endPlugin \
	-fork3 -FilterSiteBuilderPlugin -bedFile ${BKGDK} -endPlugin \
	-input1 -KinshipPlugin -method Centered_IBS -endPlugin \
	-fork4 -export ${MOAK}_CentIBS_K  -exportType SqrMatrixBin -input2 \
	-fork5 -export ${BKGDK}_CentIBS_K -exportType SqrMatrixBin -input3 \
	-combine6 -input2 -input3 -SubtractDistanceMatrixPlugin -wholeMatrix $WGKM -endPlugin \
	-export CentIBS_Rest_K -exportType SqrMatrixBin

###

outname=$1 #(Full_Phenotype_Results/LDAK_EG80WWall_full)
Klist=$2 #(kinship_list.txt)
pheno=$3 #(path to phenotype file)

module load ldak
ldak --reml $outname --mgrm $Klist --mpheno -1 --pheno NAM-ril-phenotypes/uniq_all_NAM_phenos.plink --kinship-details NO --dentist YES --constrain YES
#echo "$TIMESTAMP Finished VCAP pipeline" >> VCAPtest.log