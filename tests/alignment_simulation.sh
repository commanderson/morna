#!/bin/bash

#PREREQ:
#Set up ERR922713_1.fastq.gz and ERR922713_2.fastq.gz
#and downsample_fastqs.py in the current directory

#Prepare files and directories
mkdir alignments

#Get hisat2 index
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg38.tar.gz
tar -zxvf hg38.tar.gz

#Establish ground truth from full-file alignment
#First pass of full-file alignment
hisat2 -x hg38/genome -1 ERR922713_1.fastq.gz -2 ERR922713_2.fastq.gz --novel-splicesite-outfile alignments/pass_1_novel_splicesites | samtools view -bS > alignments/pass_1_alignment.bam

#Second pass of full-file alignment
hisat2 -x hg38/genome -1 ERR922713_1.fastq.gz -2 ERR922713_2.fastq.gz --novel-splicesite-infile alignments/pass_1_novel_splicesites | samtools view -bS > alignments/pass_2_alignment.bam

#Downsampling operation
python downsample_fastqs.py -1 ERR922713_1.fastq.gz -2 ERR922713_2.fastq.gz -d 40631462 -o 406314621

mv ERR922713_1.fastq.gz.downsamp.gz ERR922713_1_downsamp.fastq.gz
mv ERR922713_2.fastq.gz.downsamp ERR922713_2_downsamp.fastq.gz

#ON DOWNSAMPLED FILES:
mkdir alignments/downsampled
#hisat 2 2pass 
mkdir alignments/downsampled/2pass
#First pass of downsampled file
hisat2 -x hg38/genome -1 ERR922713_1_downsamp.fastq.gz -2 ERR922713_1_downsamp.fastq.gz --novel-splicesite-outfile alignments/downsampled/2pass/ds_pass_1_novel_splicesites | samtools view -bS > alignments/downsampled/2pass/ds_pass_1_alignment.bam

#Second pass of downsampled alignment
hisat2 -x hg38/genome -1 ERR922713_1_downsamp.fastq.gz -2 ERR922713_1_downsamp.fastq.gz --novel-splicesite-infile alignments/downsampled/2pass/ds_pass_1_novel_splicesites | samtools view -bS > alignments/downsampled/2pass/ds_pass_2_alignment.bam

#hisat 2 morna 
mkdir alignments/downsampled/morna
python ../morna.py align -x ../v2master/v2master -v -d -m --aligner hisat2 -i hg38/genome -1 ERR922713_1_downsamp.fastq.gz -2 EERR922713_1_downsamp.fastq.gz -p1 alignments/downsampled/morna/morna_p1.sam -p2 alignments/downsampled/morna/morna_p2.sam --junction-file ../pooled_combined_jns.tsv.gz

#hisat 2 annotated
#download gene annotation
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz
gunzip gencode.v26.annotation.gtf.gz

#create splice sites file
berk=$(which hisat2_extract_splice_sites.py)
python $berk gencode.v26.annotation.gtf > splicesites.txt

#run alignment with known splicesites
mkdir alignments/downsampled/annotated
hisat2 -x hg38/genome -S alignments/downsampled/annotated/ds_ann_alignment -1 ERR922713_1.fastq.gz.downsamp -2 ERR922713_2.fastq.gz.downsamp --known-splicesite-infile splicesites.txt
