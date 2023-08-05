#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=3:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue

module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load plink

plink --bfile NAM_rils_SNPs_allchr --clump bQTL.padj_probgeno.numeric.bed --clump-kb 0.150 --clump-r2 0.99 --clump-field P --clump-p1 0.1 --clump-p2 0.3 -out NAM_rils_SNPs.padj_probgeno.clumped
#plink --bfile NAM_rils_SNPs_allchr --clump bQTL.padj_probgeno.bed --out defaulttest.clump --clump-p1 0.1 --clump-p2 0.3
