#!/bin/bash
#For merged reads pipeline (has no quality cut off)
Fastq_DIR=$1
Samples=$2
ID=$3
SeqPurge_DIR=$4
lane=$5

module load miniconda3
source activate ngsbitsenv

SeqPurge \
 -min_len 20 \
 -threads 16 \
 -qcut 0 \
 -prefetch 5000 \
 -in1 $Fastq_DIR/${Samples}${ID}_1_L${lane}.fq.gz \ 
 -in2 $Fastq_DIR/${Samples}${ID}_2_L${lane}.fq.gz \
 -out1 $SeqPurge_DIR/${Samples}${ID}_1.trim.gz \
 -out2 $SeqPurge_DIR/${Samples}${ID}_2.trim.gz \
 -summary $SeqPurge_DIR/${Samples}_1.trim.stats \
 -progress 120000 \
 -qc $SeqPurge_DIR/${Samples}_1.trim.qc
 
source deactivate

#For paired-end reads pipeline (default quality cut off)

SeqPurge \
 -min_len 20 \
 -threads 16 \
 -prefetch 5000 \
 -in1 $Fastq_DIR/${Samples}${ID}_1_L${lane}.fq.gz \ 
 -in2 $Fastq_DIR/${Samples}${ID}_2_L${lane}.fq.gz \
 -out1 $SeqPurge_DIR/${Samples}${ID}_1.trim.gz \
 -out2 $SeqPurge_DIR/${Samples}${ID}_2.trim.gz \
 -summary $SeqPurge_DIR/${Samples}_1.trim.stats \
 -progress 120000 \
 -qc $SeqPurge_DIR/${Samples}_1.trim.qc


###SLURM Command
now=$(date +"%T")
for i in A619 ; do
    for j in w d ; do
        for k in 21 22 23 ; do 
        	echo Starting seqpurge with ${i}...${now}
        	bash scripts/seqpurge_before_merge.sh 00_B73x${i}_data ${i} ${j}${k} 01a_seqpurge_before_merge 4
        	echo Finished with before merge\, moving onto for paired end
        	bash scripts/seqpurge_for_paired_end.sh 00_B73x${i}_data ${i} ${j}${k} 01b_seqpurge_for_paired_end 4
       	 echo Finished with seqpurge for paired end 
       	 done
    done
    echo Finished with ${i} ...${now} 
done 00_B73xCML69_data/CML69w4_2.trim.gz

###For the CML69 and Oh43 files that have different naming conventions
#For the merged reads
Fastq_DIR=$1
Samples=$2
ID=$3
SeqPurge_DIR=$4

if [[ -f ${Fastq_DIR}/${Samples}${ID}_1.trim.gz ]] ; then
SeqPurge \
  -min_len 20 \
  -threads 16 \
  -qcut 0 \
  -prefetch 5000 \
  -in1 ${Fastq_DIR}/${Samples}${ID}_1.trim.gz \ 
  -in2 ${Fastq_DIR}/${Samples}${ID}_2.trim.gz \
  -out1 ${SeqPurge_DIR}/${Samples}${ID}_1.trim.gz \
  -out2 ${SeqPurge_DIR}/${Samples}${ID}_2.trim.gz \
  -summary ${SeqPurge_DIR}/${Samples}${ID}_1.trim.stats \
  -progress 120000 \
  -qc ${SeqPurge_DIR}/${Samples}${ID}_1.trim.qc

#For paired end reads
SeqPurge \
  -min_len 20 \
  -threads 16 \
  -prefetch 5000 \
  -in1 ${Fastq_DIR}/${Samples}${ID}_1.trim.gz \ 
  -in2 ${Fastq_DIR}/${Samples}${ID}_2.trim.gz \
  -out1 ${SeqPurge_DIR}/${Samples}${ID}_1.trim.gz \
  -out2 ${SeqPurge_DIR}/${Samples}${ID}_2.trim.gz \
  -summary ${SeqPurge_DIR}/${Samples}${ID}_1.trim.stats \
  -progress 120000 \
  -qc ${SeqPurge_DIR}/${Samples}${ID}_1.trim.qc

fi