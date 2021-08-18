#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN #Can be edited out if submitting many jobs
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --exclusive

if [[ ! -f kinship.list ]] ; then 
	echo MOA.SNPs.bed_CentIBS_K > kinship.list
	echo BKGD.SNPs.bed_CentIBS_K >> kinship.list
	echo CentIBS_Rest_K >> kinship.list
fi


bash /work/LAS/mhufford-lab/snodgras/VCAP_runLDAK.sh MBRest_model kinship.list /work/LAS/mhufford-lab/snodgras/VCAP_Runs/JeffsSimulations/simulatedphenotypes_highMOAh2.plink