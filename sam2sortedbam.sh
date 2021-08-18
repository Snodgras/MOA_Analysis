#!/bin/bash
module load samtools

SAM="$1"
samtools view --threads 36 -b -o ${SAM%.*}.bam ${SAM}
samtools sort -o ${SAM%.*}_sorted.bam -T ${SAM%.*}_temp --threads 36 ${SAM%.*}.bam

#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=8:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=3   # number of nodes
#SBATCH --ntasks-per-node=16   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue

bash scripts/sam2sortedbam.sh aln_bwa_d5_1.trim.100bp_5prime_d5_2.trim.100bp_3prime.sam
bash scripts/sam2sortedbam.sh aln_bwa_d5_1.trim_d5_2.trim.sam
bash scripts/sam2sortedbam.sh d5_1_mapped_rnd2.sam #untrimmed
bash scripts/sam2sortedbam.sh d5_mapped_rnd2.sam #trimmed to 100bp