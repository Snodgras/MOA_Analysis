#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load plink

plink --file NAM_rils_SNPs_chr10.plk --merge-list allfiles.txt --make-bed --out NAM_rils_SNPs_allchr
