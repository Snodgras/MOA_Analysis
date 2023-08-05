#!/bin/bash
#Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=1:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue

#bash SNPselection.sh MOA_EG80WWSNPs_AF_dist.txt nonEG80WWSNPs_AF_dist.txt

module load r r-devtools 
module load r-dplyr
module load r-tidyr

./SNPselection.R MOA_EG80WW_0423SNPs_AF_dist.txt nonEG80WWSNPs_AF_dist.txt VCAP_Runs/1SNPperPeak/
