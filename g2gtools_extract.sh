#!/bin/bash
#For g2gtools need to create a condo environment

#module load miniconda2/4.3.30-gy2db2b
#conda create --prefix /work/LAS/mhufford-lab/snodgras/bin/conda/g2gtools
#source activate /work/LAS/mhufford-lab/snodgras/bin/conda/g2gtools
#conda install -c kbchoi g2gtools
#g2gtools
ref=$1
ref_nam=${ref%.fasta}
peaks=$2 #d5_NGmerge_100trimmed.fastq_ref_B73.GPS_events.bed
peak_nam=${peaks%.fastq*}


module load miniconda2/4.3.30-gy2db2b
source activate /work/LAS/mhufford-lab/snodgras/bin/conda/g2gtools
g2gtools extract -i ${ref} -b ${peaks} > ${peak_nam}_${ref_nam}.fasta
echo Finished with ${ref_nam} and ${peak_nam}
source deactivate /work/LAS/mhufford-lab/snodgras/bin/conda/g2gtools