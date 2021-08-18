#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=6:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue

if [[ ! -f MOA.SNPs.bed ]] ; then 
	tail -n +2 BKGD_distrib.SNPs.txt | cut -f 1-3 > BKGD.SNPs.bed
	tail -n +2 MOA_distrib.SNPs.txt | cut -f 1-3 > MOA.SNPs.bed
fi


bash /work/LAS/mhufford-lab/snodgras/VCAP_mkKinship.sh /work/LAS/mhufford-lab/snodgras/NAM_RILs_0423projectedSNPs_CentIBS_wholegenomekinship.txt.gz MOA.SNPs.bed BKGD.SNPs.bed