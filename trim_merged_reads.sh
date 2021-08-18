#!/bin/bash

#removing reads longer 100 bp from merged fastq

input=$1

module load bioawk
zcat ${input} | bioawk -c fastx 'length($seq)<=100 {print "@"$name" "$comment"\n"$seq"\n+\n"$qual}' > ${input}_100trimmed.fastq

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

for i in *NGmerge.gz; do bash scripts/trim_merged_reads.sh $i ; done

###File = 02.5a_filtered_merged_reads.sh
#!/bin/bash

#removing reads longer 100 bp from merged fastq

input=$1
nam=${input#02_merged_reads/}

module load bioawk
zcat ${input} | bioawk -c fastx 'length($seq)<=100 {print "@"$name" "$comment"\n"$seq"\n+\n"$qual}' > 02.5a_filtered_merged_reads/${nam}_100trimmed.fastq
