#!/bin/bash

module load bamtools
input=$1

bamtools filter -tag H0:i:2 -in ${input} -out ${input}_H0_2.bam #H0 = perfect hits
bamtools filter -tag NH:i:2 -in ${input} -out ${input}_NH_2.bam #NH = num reported alignments of query
bamtools filter -tag X0:i:2 -in ${input} -out ${input}_X0i_2.bam #X0 = num of best hits

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

bash scripts/filterbybamtags.sh aln_bwa_d5_NGmerge_100trimmed.fastq_sorted.bam
bash scripts/filterbybamtags.sh d5_NGmerged_trimmed100_mapped_rnd2_sorted.bam

echo ${input} has $(grep -c X0:i:1 ${input}) lines have X0:i:1 tag (number of best hits) >> sam_tag_info.txt
echo ${input} has $(grep -c XT:A:U ${input}) lines have XT:A:U tag (number of unique hits, no quality check) >> sam_tag_info.txt
echo ${input} has $(grep -c H0:i:2 ${input}) lines have H0:i:2 tag (number of perfect hits =?= 2) >> sam_tag_info.txt
echo ${input} has $(grep -c H0:i:1 ${input}) lines have H0:i:1 tag (number of perfect hits =?= 1) >> sam_tag_info.txt
echo ${input} has $(grep -c NH:i:2 ${input}) lines have NH:i:2 tag (number of reported alignments of query =?= 2) >> sam_tag_info.txt
echo ${input} has $(grep -c NH:i:1 ${input}) lines have NH:i:1 tag (number of reported alignments of query =?= 1) >> sam_tag_info.txt
echo ${input} has $(grep -c XT:A:M ${input}) lines have XT:A:M tag (mapped in multiple places along ref genome) >> sam_tag_info.txt

aln_bwa_d5_NGmerge_100trimmed.fastq.sam
d5_NGmerged_trimmed100_mapped_rnd2.sam