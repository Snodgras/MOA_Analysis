# MOA_Analysis
Mapping MOA reads to diploid genomes

## Data from:
ftp://ftp.mpipz.mpg.de/

## Creating reference from 2 genome fastas
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

## Alignment
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

## What Thomas suggests
1. Map reads to the diploid genome with BWA-mem. 
	* To cat'ed genomes _DONE_
	* To individual parental genomes _DONE_
2. use unfiltered bam files and run GEM (ChIP-Seq peak calling software) in GPS mode to determine all MOA binding peaks. 
	* GEM installation? _DONE_
	* make chr length files with `bioawk -c fastx '{print $name"\t"length($seq)}' ref_Oh43.fasta > ref_Oh43-length.txt`
	Testing: _DONE_
3. Use e.g g2gtools or deep tools to extract all reads that overlap with GEM MOA peaks 
	* Have to make a database from the annotation files first
	Use cufflinks to transform the gff to gtf
	gffread -E Zm-Oh43-REFERENCE-NAM-1.0_Zm00039a.1.gff -T -o Oh43.gtf
	
	Tried making the database in an interactive session, got this error: 
	(/work/LAS/mhufford-lab/snodgras/bin/conda/g2gtools) [snodgras@condofree019 MOA]$ g2gtools gtf2db -i Oh43.gtf -o Oh43.db
GTF FILE: /work/LAS/mhufford-lab/snodgras/MOA/Oh43.gtf
DB File: /work/LAS/mhufford-lab/snodgras/MOA/Oh43.db
Parsing GTF file...
Traceback (most recent call last):
  File "/work/LAS/mhufford-lab/snodgras/bin/conda/g2gtools/bin/g2gtools", line 4, in <module>
    __import__('pkg_resources').run_script('g2gtools==0.1.31', 'g2gtools')
  File "/work/LAS/mhufford-lab/snodgras/bin/conda/g2gtools/lib/python2.7/site-packages/pkg_resources/__init__.py", line 666, in run_script
    self.require(requires)[0].run_script(script_name, ns)
  File "/work/LAS/mhufford-lab/snodgras/bin/conda/g2gtools/lib/python2.7/site-packages/pkg_resources/__init__.py", line 1462, in run_script
    exec(code, namespace, namespace)
  File "/work/LAS/mhufford-lab/snodgras/bin/conda/g2gtools/lib/python2.7/site-packages/g2gtools-0.1.31-py2.7.egg-info/scripts/g2gtools", line 117, in <module>
    G2GToolsApp()
  File "/work/LAS/mhufford-lab/snodgras/bin/conda/g2gtools/lib/python2.7/site-packages/g2gtools-0.1.31-py2.7.egg-info/scripts/g2gtools", line 75, in __init__
    getattr(self, args.command)()
  File "/work/LAS/mhufford-lab/snodgras/bin/conda/g2gtools/lib/python2.7/site-packages/g2gtools-0.1.31-py2.7.egg-info/scripts/g2gtools", line 99, in gtf2db
    g2gtools.g2g_commands.command_gtf2db(sys.argv[2:], self.script_name + ' gtf2db')
  File "/work/LAS/mhufford-lab/snodgras/bin/conda/g2gtools/lib/python2.7/site-packages/g2gtools/g2g_commands.py", line 151, in command_gtf2db
    gtf2db(args.input, args.output)
  File "/work/LAS/mhufford-lab/snodgras/bin/conda/g2gtools/lib/python2.7/site-packages/g2gtools/gtf_db.py", line 249, in gtf2db
    ensembl_id = record.attributes['exon_id']
KeyError: 'exon_id'
	
4. Convert overlapping bam files with alntools to an incidence matrix so they can be used as emase input. 
5. Use bet file from GEM to generate the length and position files of MOA peaks required as inputs for emase.
6. Run emase, use output with region quantification of peaks to determine allele-specific peaks.

## alntools 

## emase

