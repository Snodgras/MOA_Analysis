#Notes about MOA 

##Data from:
ftp://ftp.mpipz.mpg.de/

##Creating reference from 2 genome fastas
* need to change fasta id names so there's no conflicts
bioawk `___pattern____ {print statement}`
bioawk -c fastx '{print ">Tzi8-"$name"\n"$seq}' <file name> | fold | less
script: `namefastaids.sh`
_DONE_
* then cat the genomes together
`cat ref_B73.fasta ref_Oh43.fasta >> ref_B73Oh43.fasta`
_DONE_
* then create an index file
script: `indexref.sh`
_DONE_

##Alignment
* optimize it to 50-70%
* hisat2
output is a bam file that can be used for alntools

script: `alignment.sh`
_DONE_

* figure out the summary statistics for bwa-mem alignments
	* if using a 3rd party, run on both bwa and hisat2 sam files _DONE_
	* could use picard tools (alignment summary metrics, try default first)
* trim the read files to see if it improves the hisat2 _DONE_
* try reducing the stringency of mapping parameters _DONE_
* do default options for bwa-mem _DONE_

```trim_galore --hardclip5 100 <R1 file>
trim_galore --hardclip3 100 <R2 file>```

* Try merging the reads first (NGmerge) _DONE_
`merge_reads.sh` then `trim_merged_reads.sh`
	* QC the merged reads _DONE_
	* remove longer than 100bp reads _DONE_
	* re-run both hisat and bwa-mem --> summary statistics
	`aln_bwa_d5_NGmerge_100trimmed.fastq.sam	`
	`d5_NGmerged_trimmed100_mapped_rnd2.sam`
	`_sorted.bam`
To make the comparable mapping results
head -n 7 mapping stats > tmp.txt
tail -n +8 mapping states | awk -v OFS='\t' -v n=merged_hisat2 '{print n,$0}' >> tmp.txt 

##What Thomas suggests
1. Map reads to the diploid genome with BWA-mem. 
	* To cat'ed genomes _DONE_
	* To individual parental genomes _DONE_
2. use unfiltered bam files and run GEM (ChIP-Seq peak calling software) in GPS mode to determine all MOA binding peaks. 
	* GEM installation? _DONE_
	* make chr length files with `bioawk -c fastx '{print $name"\t"length($seq)}' ref_Oh43.fasta > ref_Oh43-length.txt`
	Testing: _Submitted batch job 797938_
3. Use e.g g2gtools or deep tools to extract all reads that overlap with GEM MOA peaks 
4. Convert overlapping bam files with alntools to an incidence matrix so they can be used as emase input. 
5. Use bet file from GEM to generate the length and position files of MOA peaks required as inputs for emase.
6. Run emase, use output with region quantification of peaks to determine allele-specific peaks.

##alntools 

##emase

