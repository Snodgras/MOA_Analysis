#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue

module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load tassel

#run_pipeline.pl -Xmx64g -fork1 -importGuess NAM_rils_SNPs_chr10.hmp.txt -export -exportType VCF 

run_pipeline.pl -Xmx64g -fork1 -importGuess NAM_rils_SNPs_allchr.hmp.txt.gz -export -exportType VCF
