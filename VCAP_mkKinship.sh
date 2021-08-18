#!/bin/bash

WGKM=$1 #name of the whole genome kinship matrix (NAMRILs_projectedSNPs_CentIBS_wholegenomekinship.txt.gz)
MOAK=$2 #path/to/myBedFile1.bed
BKGDK=$3

module load tassel
module load samtools

if [[ ! -f "$WGKM" ]] ; then
	echo "$TIMESTAMP Beginning full kinship matrix estimation" > VCAPtest.log ;
	run_pipeline.pl \
		-Xmx64g \
		-importGuess /work/LAS/mhufford-lab/snodgras/projections/projected_hmpfiles/NAM_rils_SNPs.ALL.best-markers.chr-ALL.projected.hmp.txt  \
		-KinshipPlugin \
		-method Centered_IBS \
		-endPlugin \
		-export $WGKM \
		-exportType SqrMatrix
	bgzip -c /work/LAS/mhufford-lab/snodgras/projections/projected_hmpfiles/NAM_rils_SNPs.ALL.best-markers.chr-ALL.projected.hmp.txt > /work/LAS/mhufford-lab/snodgras/projections/projected_hmpfiles/NAM_rils_SNPs.ALL.best-markers.chr-ALL.projected.hmp.txt.gz
	run_pipeline.pl -Xmx64g -LIXPlugin -createIndex /work/LAS/mhufford-lab/snodgras/projections/projected_hmpfiles/NAM_rils_SNPs.ALL.best-markers.chr-ALL.projected.hmp.txt.gz
	echo "$TIMESTAMP Finished full kinship matrix." >> VCAPtest.log ;
fi

#This works for 2 kinship matrices

run_pipeline.pl -Xmx64g \
	-fork1 -importGuess /work/LAS/mhufford-lab/snodgras/projections/projected_hmpfiles/NAM_rils_SNPs.ALL.best-markers.chr-ALL.projected.hmp.txt.gz.lix \ 
	-fork2 -FilterSiteBuilderPlugin -bedFile ${MOAK} -endPlugin \
	-input1 -KinshipPlugin -method Centered_IBS -endPlugin \
	-fork3 -FilterSiteBuilderPlugin -bedFile ${BKGDK} -endPlugin \
	-input1 -KinshipPlugin -method Centered_IBS -endPlugin \
	-fork4 -export ${MOAK}_CentIBS_K  -exportType SqrMatrixBin -input2 \
	-fork5 -export ${BKGDK}_CentIBS_K -exportType SqrMatrixBin -input3 \
	-combine6 -input2 -input3 -SubtractDistanceMatrixPlugin -wholeMatrix $WGKM -endPlugin \
	-export CentIBS_Rest_K -exportType SqrMatrixBin

