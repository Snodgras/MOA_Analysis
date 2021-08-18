#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=6:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue

module load bioawk

bioawk -c fastx '{print ">B73_AB10-"$name"\n"$seq}' B73_AB10.pseudomolecules-v1.fasta > ref_B73_AB10.fasta
bioawk -c fastx '{print ">B73-"$name"\n"$seq}' B73.PLATINUM.pseudomolecules-v1.fasta > ref_B73.fasta
bioawk -c fastx '{print ">B97-"$name"\n"$seq}' B97.pseudomolecules-v1.fasta > ref_B97.fasta
bioawk -c fastx '{print ">CML103-"$name"\n"$seq}' CML103.pseudomolecules-v1.fasta > ref_CML103.fasta
bioawk -c fastx '{print ">CML228-"$name"\n"$seq}' CML228.pseudomolecules-v1.fasta > ref_CML228.fasta
bioawk -c fastx '{print ">CML277-"$name"\n"$seq}' CML277.pseudomolecules-v1.fasta > ref_CML277.fasta
bioawk -c fastx '{print ">CML247-"$name"\n"$seq}' CML247.pseudomolecules-v1.fasta > ref_CML247.fasta
bioawk -c fastx '{print ">CML322-"$name"\n"$seq}' CML322.pseudomolecules-v2.fasta > ref_CML322.fasta
bioawk -c fastx '{print ">CML333-"$name"\n"$seq}' CML333.pseudomolecules-v1.fasta > ref_CML333.fasta
bioawk -c fastx '{print ">CML52-"$name"\n"$seq}' CML52.pseudomolecules-v2.1.fasta > ref_CML52.fasta
bioawk -c fastx '{print ">CML69-"$name"\n"$seq}' CML69.pseudomolecules-v1.fasta > ref_CML69.fasta
bioawk -c fastx '{print ">HP301-"$name"\n"$seq}' HP301.pseudomolecules-v1.fasta > ref_HP301.fasta
bioawk -c fastx '{print ">IL14H-"$name"\n"$seq}' IL14H.pseudomolecules-v2.fasta > ref_IL14H.fasta
bioawk -c fastx '{print ">Ki11-"$name"\n"$seq}' Ki11.pseudomolecules-v1.fasta > ref_Ki11.fasta
bioawk -c fastx '{print ">Ki3-"$name"\n"$seq}' Ki3.pseudomolecules-v1.fasta	 > ref_Ki3.fasta
bioawk -c fastx '{print ">Ky21-"$name"\n"$seq}' Ky21.pseudomolecules-v1.fasta > ref_Ky21.fasta
bioawk -c fastx '{print ">M162W-"$name"\n"$seq}' M162W.pseudomolecules-v1.fasta > ref_M162W.fasta
bioawk -c fastx '{print ">M37W-"$name"\n"$seq}' M37W.pseudomolecules-v1.fasta > ref_M37.fasta
bioawk -c fastx '{print ">Mo18W-"$name"\n"$seq}' Mo18W.pseudomolecules-v1.fasta > ref_Mo18W.fasta
bioawk -c fastx '{print ">MS71-"$name"\n"$seq}' MS71.pseudomolecules-v1.fasta > ref_MS71.fasta
bioawk -c fastx '{print ">NC350-"$name"\n"$seq}' NC350.pseudomolecules-v1.fasta  > ref_NC350.fasta
bioawk -c fastx '{print ">NC358-"$name"\n"$seq}' NC358.pseudomolecules-v1.fasta > ref_NC358.fasta
bioawk -c fastx '{print ">Oh43-"$name"\n"$seq}' Oh43.pseudomolecules-v1.fasta > ref_Oh43.fasta
bioawk -c fastx '{print ">Oh7b-"$name"\n"$seq}' Oh7b.pseudomolecules-v2.fasta > ref_Oh7b.fasta
bioawk -c fastx '{print ">P39-"$name"\n"$seq}' P39.pseudomolecules-v2.fasta > ref_P39.fasta
bioawk -c fastx '{print ">Tx303-"$name"\n"$seq}' Tx303.pseudomolecules-v1.fasta > ref_Tx303.fasta
bioawk -c fastx '{print ">Tzi8-"$name"\n"$seq}' Tzi8.pseudomolecules-v1.fasta > ref_Tzi8.fasta
