#!/bin/bash

#Make the chain files
#Using flo from: https://github.com/wurmlab/flo
#needs a configuration file (.yaml)
#(.yaml) includes: 
	#location of source, target assembly in FASTA format,
	#Location of GFF3 file(s) of annotations on source assembly (omission will make flo stop after creating the chain file)
	#BLAT parameters (opt) for if you're chaining two different species
	#Number of CPU cores to use, cannot be greater than the number of scaffolds in the target assembly

#Here, it's important to note that flo can only work with transcripts and their child exons 
#and CDS. Transcripts can be annotated as: mRNA, transcript, or gene. 
#However, if you have a 'gene' annotation for each transcript, you will need to remove that
#For genes with multiple transcripts, select the largest transcript for each gene
################################
module load flo

#need some way to specify the path to the source/target_fa in the .yaml

#modify the gffs
/path/to/flo/gff_remove_feats.rb gene xx_genes.gff > xx_transcripts.gff #deals with extra features
/path/to/flo/gff_longest_transcripts.rb xx_genes.gff > xx_longest_transcripts.gff #selects longest transcript

#runs flo to get chain files
rake -f /path/to/flo/Rakefile 

#Make VCF using CrossMap
#From: http://crossmap.sourceforge.net/

CrossMap.py vcf chain_file input_VCF_file target_genome_file output_file