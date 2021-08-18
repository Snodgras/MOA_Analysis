#!/bin/bash
# load modules
module load samtools
module load bamtools
#module load bedtools2
module load hisat2
#module load star
# variables
index="$1" #the basename of the index (not including the .#.ht2)
R1="$2" #comma separated list of files with unpaired reads to be aligned (1.fq,2.fq,3.fq...)
#nam="$4"
#R2="$3"
nam="$3"

hisat2 -p 16 --mp 1,1  \
    --no-spliced-alignment \
    --end-to-end \
    --rdg 10000,10000 \
    --rfg 10000,10000 \
    -q -x ${index} -U ${R1} 1> ${nam}_mapped_paired.sam 2> ${nam}_mapping_stats.txt
    #-q -x ${index} -1 ${R1} -2 ${R2}  1> ${nam}_mapped_paired.sam 2> ${nam}_mapping_stats.txt
#--no-softclip
#STAR --runMode alignReads --runThreadN 16 --genomeDir ${index} --outFileNamePrefix ${nam}-star --readFilesIn ${markers}
#--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3
# sort and covert to bam
#samtools view -b -o ${nam}_mapped.bam ${nam}_mapped.sam
#${nam}-starAligned.out.sam
#samtools sort -o ${nam}_sorted.bam ${nam}_mapped.bam
#samtools view -q 30 -b ${nam}_sorted.bam > ${nam}_sorted-q30.bam
#bamtools filter -tag XM:0 -in ${nam}_sorted-q30.bam -out ${nam}_sorted_noMismatch.bam

#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=3   # number of nodes
#SBATCH --ntasks-per-node=16   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue

bash scripts/alignment.sh ref_B73Oh43 d5_1.trim d5_1