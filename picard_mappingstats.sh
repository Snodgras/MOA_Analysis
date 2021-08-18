#!/bin/bash
module load picard
ref=$1 #reference fasta file
input=$2 #sam or bam file
out=mapping_stats_${input%.bam} #what file it should write to

picard CollectAlignmentSummaryMetrics \
          R=$1 \
          I=$2 \
          O=${out}.txt

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

bash scripts/picard_mappingstats.sh ref_B73Oh43.fasta aln_bwa_d5_1.trim.100bp_5prime_d5_2.trim.100bp_3prime_sorted.bam
bash scripts/picard_mappingstats.sh ref_B73Oh43.fasta aln_bwa_d5_1.trim_d5_2.trim_sorted.bam
bash scripts/picard_mappingstats.sh ref_B73Oh43.fasta d5_1_mapped_rnd2_sorted.bam #untrimmed
bash scripts/picard_mappingstats.sh ref_B73Oh43.fasta  d5_mapped_rnd2_sorted.bam #trimmed to 100bp