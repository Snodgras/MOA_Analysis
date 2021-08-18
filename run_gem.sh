#!/bin/bash
module load jdk

#Experiment=HybMoa_0819_WWvsDS
#ECO=Oh43
#Samples=B73vs"$ECO"
#treat=ds
#ctrl=ww
#AGP=AGPv4
#MAPQ=13 
pvalue=2
fold=3
nd=1
smooth=3

input=$1 #path to aligned read file (aln_bwa_d5_NGmerge_100trimmed.fastq_ref_B73.sam)
output=${input#aln_bwa_}
out=${output%.fastq} 
output=${out%.sam} #where to put the output
length=$2 #path to length file (ref_B73-length.txt)



#Output_DIR=$netscratch_DIR/Thomas/Results/$Experiment/GEM/BWA/"$Samples".T_"$treat".C_"$ctrl".q"$pvalue"."$fold"fc.nd"$nd".sm"$smooth"
#Input_DIR=$netscratch_DIR/Thomas/Results/$Experiment/BWA/bias_corrected/merged_replica
#Genome_DIR=$netscratch_DIR/Michael_Thomas/Genomes/Zea_mays/B73/$AGP/chr1_10
​
if [ ! -d $Output_DIR ]
then
mkdir -p $Output_DIR 
fi
​
cd $Output_DIR
echo "Selecting best mapping reads of $Samples mapped against B73 and $ECO!"

java -Xmx128G -XX:+UseParallelGC -XX:ParallelGCThreads=10 \
	-jar gem/gem.jar \
	--d gem/'Read_Distribution_default.txt' \ 
	--g ${length} \
	--min 5000 \
	--mrc 5000 \
	--nd ${nd} \
	--smooth ${smooth} \
	--t 16 \
	--expt ${input} \
	--f SAM \
	--sl \
	--out ${output} \
	--q "$pvalue" --outBED  --fold "$fold" --nrf
​
#--d path to the read distribution model file
#--g path to genome information file (chromosome lengths)
#--min Min number of events called by GPS/GEM to get read distribution and cont. to next round
#--mrc Max read count on a base position, helps remove duplicate reads
#--nd noise distribution model 1=uniform noise model
#--smooth width (bp) to smooth read distribution
#--t number of threads to run GMEM in parallel, =< physical CPU
#--expt1 Path to aligned read files for experiment
#--f read file format
#--sl sort GEM output by location
#--out Output folder and file name prefix
#--q significance level for q-value
#--outBED Output binding events in bed format
#--fold fold cutoff to filter predicted events
#--nrf do not filter duplicate reads