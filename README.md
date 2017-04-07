# morna
Tools for searching for existing RNA-seq samples with expression patterns similar to new sample
Intropolis data was downloaded from figshare as part of "intropolis: exon-exon junctions across publicly available RNA-seq samples on the Sequence Read Archive"
https://figshare.com/articles/intropolis_exon-exon_junctions_across_publicly_available_RNA-seq_samples_on_the_Sequence_Read_Archive/3811680/1
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


# Pooling Junction Lists
In order to create a master index using morna index, we had to create a pooled junction list combining the junction contents of tcga, sra, and gtex. Additionally, in order to include metadata information with this index, we had to combine the metadata files which mapped the sample ids in each source to a set of metadata

This was achieved by downloading the files using the following commands:
intropolis.v1.hg19.tsv.gz:
wget https://ndownloader.figshare.com/files/5935809
first_pass_gtex_junctions.tsv.gz and  first_pass_tcga_junctions.tsv.gz:
available at intropolis.rail.bio.

metadata files
gtex_metadata.tsv, sra_metadata.tsv, tcga_metadata.tsv:
available at intropolis.rail.bio

Those files were renamed and then concatenated with the following commands:

mv intropolis.v1.hg19.tsv.gz sra_junctions.tsv.gz
mv first_pass_gtex_junctions.tsv.gz gtex_junctions.tsv.gz
mv first_pass_tcga_junctions.tsv.gz tcga_junctions.tsv.gz

for f in gtex sra tcga; 
do gzip -cd ${f}_junctions.tsv.gz|awk -v src="$f" '{print $0 "\t" src}'; 
done > combined_junctions.tsv

gzip combined_junctions.tsv
gzip -cd combined_junctions.tsv.gz | sort > sorted_cjs.tsv

gzip sorted_cjs.tsv
python pool_junction_lists.py -i sorted_cjs.tsv.gz -s sra tcga gtex -m -n 25000 -o pooled_combined_jns.tsv.gz

# Creating Master Index:
mkdir mastershards
python morna.py index --intropolis pooled_combined_jns.tsv.gz -x mastershards/master_index --n-trees 10 -v -b 32768 -m sra_tcga_gtex_metadata.tsv 


utilities:
all python scripts run with Python 2.7.12
[GCC 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2336.11.00)]

