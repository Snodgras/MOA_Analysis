#!/bin/bash
ml bwa
ref=$1
R1=$2
#R2=$3
nam=${R1%.fq}_${ref%.fasta}  #_${R2%.fq}

bwa mem -t 16 ${ref} ${R1} ${R2} > aln_bwa_${nam}.sam

#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=8:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue

bash scripts/bwaalign.sh ref_B73Oh43.fasta d5_1.trim.fq d5_2.trim.fq
bash scripts/bwaalign.sh ref_B73Oh43.fasta d5_1.trim.100bp_5prime.fq d5_2.trim.100bp_3prime.fq