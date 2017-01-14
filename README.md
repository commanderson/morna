# rna-sim-search
Tools for searching for existing RNA-seq samples with expression patterns similar to new sample
Intropolis data was downloaded from figshare as part of "intropolis: exon-exon junctions across publicly available RNA-seq samples on the Sequence Read Archive" https://figshare.com/articles/intropolis_exon-exon_junctions_across_publicly_available_RNA-seq_samples_on_the_Sequence_Read_Archive/3811680/1
Exon-exon junction data:
intropolis.v1.hg19.tsv.gz
https://ndownloader.figshare.com/files/5935809
Junction to SRA mapping data:
intropolis.idmap.v1.hg19.tsv
https://ndownloader.figshare.com/files/5935797

Created filtered files using: 
gzip -cd intropolis.v1.hg19.tsv.gz | awk -F',' 'NF >= 199'>intropolis_junctions_100samples.tsv
-To restrict to only junctions present in at least 100 samples
gzip -cd intropolis.v1.hg19.tsv.gz | awk -F',' 'NF >= 1999'>intropolis_junctions_1000samples.tsv
-To restrict to only junctions present in at least 1000 samples.

