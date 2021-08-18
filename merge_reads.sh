#!/bin/bash

#/work/LAS/mhufford-lab/snodgras/MOA/NGmerge-0.2
#Removing residual adapter pieces and merging forward and reverse reads
#Software NGmerge https://github.com/jsh58/NGmerge
R1=$1
R2=$2

/work/LAS/mhufford-lab/snodgras/MOA/NGmerge-0.2/NGmerge \
	-1 ${R1} \
	-2 ${R2} \
	-o ${R1}_${R2}_NGmerge \
	-p 0.2 \
	-m 15 \
	-d \
	-e 30 \
	-z \
	-n 48 \
	-v \
	-c ${R1}_${R2}_adapter_out

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

bash scripts/merge_reads.sh B73xOh43_data/d5_1.trim.gz B73xOh43_data/d5_2.trim.gz
bash scripts/merge_reads.sh B73xOh43_data/d6_1.trim.gz B73xOh43_data/d6_2.trim.gz
bash scripts/merge_reads.sh B73xOh43_data/d7_1.trim.gz B73xOh43_data/d7_2.trim.gz
bash scripts/merge_reads.sh B73xOh43_data/d8_1.trim.gz B73xOh43_data/d8_2.trim.gz

bash scripts/merge_reads.sh B73xOh43_data/w5_1.trim.gz B73xOh43_data/w5_2.trim.gz
bash scripts/merge_reads.sh B73xOh43_data/w6_1.trim.gz B73xOh43_data/w6_2.trim.gz
bash scripts/merge_reads.sh B73xOh43_data/w7_1.trim.gz B73xOh43_data/w7_2.trim.gz
bash scripts/merge_reads.sh B73xOh43_data/w8_1.trim.gz B73xOh43_data/w8_2.trim.gz
echo Finished merging B73xOh43 reads

bash scripts/merge_reads.sh B73xCML69_data/d1_1.trim.gz B73xCML69_data/d1_2.trim.gz
bash scripts/merge_reads.sh B73xCML69_data/d2_1.trim.gz B73xCML69_data/d2_2.trim.gz
bash scripts/merge_reads.sh B73xCML69_data/d3_1.trim.gz B73xCML69_data/d3_2.trim.gz
bash scripts/merge_reads.sh B73xCML69_data/d4_1.trim.gz B73xCML69_data/d4_2.trim.gz

bash scripts/merge_reads.sh B73xCML69_data/w1_1.trim.gz B73xCML69_data/w1_2.trim.gz
bash scripts/merge_reads.sh B73xCML69_data/w2_1.trim.gz B73xCML69_data/w2_2.trim.gz
bash scripts/merge_reads.sh B73xCML69_data/w3_1.trim.gz B73xCML69_data/w3_2.trim.gz
bash scripts/merge_reads.sh B73xCML69_data/w4_1.trim.gz B73xCML69_data/w4_2.trim.gz
echo Finished merging B73xCML69 reads