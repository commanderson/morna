morna
=====
Morna is a collection of tools for indexing existing RNA-seq data, searching for samples with similar expression patterns, and using this information to informa alignment.
Morna search can be used to test if a sample’s global expression patterns are similar to those of samples from the same tissue or cell type and may uncover significant contamination. Morna can also aid alignment by improving detection of low-coverage splice junctions without the same sacrifices as employing a full annotation, and may be useful to impute junctions that might have been found had a sample been sequenced more deeply.

Components
--------
The three main tools morna offers are index, search, and align.

**morna index**

Using a gzipped file of junctions with associated sample ids in intropolis-like format, Morna Index can create an index of samples with junction information represented in a dimensionality-reduced space.  

**Usage**

    python morna.py index --intropolis sra_junctions.tsv.gz -x sra/sra_index --n-trees 10 -v -b 51200 -m  sra_metadata.tsv

**Features**

Creates each of the following index files for a data set of junctions:
-An index using the [Annoy](https://github.com/spotify/annoy) module index storing feature-hashed junctions by internal sample id
-A frequency index storing the frequency of each sample in the junctions file
-A map from sample ids (in the input file) to sequential internal ids
-A stats file including number of samples and size of hashed feature space
-A 100-sharded sqlite3 database storing the junction indexes found in each sample
-(Optionally) A sqlite3 database associating metadata with each sample id.

**Arguments**

'--intropolis' (required)
A file path to a gzipped file recording junctions across samples in "intropolis-like"  format

'-x' or '--basename' (Not required, default is "morna")
A base file path for the index files to create; they will be called {basename}.annoy.mor, {basename}.freq.mor, {basename}.map.mor, {basename}.stats.mor, {basename}.sh00.junc.mor through {basename}.sh99.junc.mor, and {basename}.meta.mor with respect to the order mentioned above.

'--features' (not required, default is 3000)
The integer dimensional size of the hashed feature space, used for the size of the annoy index; lower numbers will have faster indexing and searching while higher numbers will potentially give more accurate approximate nearest neighbor search results.

'--n-trees' (not required, default is 200)
An integer argument passed to the annoy module in building the annoy index, this controls the number of trees in the forest that annoy uses to efficiently search for approximate nearest neighbors. A larger value will give more accurate results, but larger indexes.

'-s' or '--sample-count' (not required, calculated if missing)
The integer number of distinct sample ids in the intropolis file. If this argument is not provided, morna index will parse through the file line by line and count distinct sample ids. As this can be cumbersome for very large files, it can save considerable time to provide this option, but you must be correct or the index will be built incorrectly.

'-t' or '--sample-threshold' (not required, default is 100)
The integer minimum number of samples in which a junction must be present to be added to the index. This will not keep it from being added to the junctions-by-sample database index.

'-b' or '--buffer-size' (not required, default 1024)
When writing the junctions-by-sample database, each sample has its own "buffer" which is a list of strings containing junction information in run-length-encoded format. When the size of theses string lists in bytes exceeds this argument, they are written to the table and cleared.

'-v' or '--verbose' (option flag)
If this option flag is included, morna index will print various status and progress messages during index operation to stdout.

**morna search**
Using an index produced by morna index (or downloaded) perform approximate nearest neighbor search to find samples which resemble a given sample in the reduced-dimensional space.

**Usage**

    samtools view 201_UHR.bam|python morna.py search -x sra/sra_index -v -d -m -f sam
    
or
    python morna.py search -x v2master/v2master -v -d -m -q 134

**Features**

Harvest junctions from incoming stream of sam-, bed-, or raw-formatted sample file, calculate its representation in the reduced-dimensional space of a specific morna index, and return the approximate nearest neighbors. 
-Can optionally include distance information and metadata for each result
-Alternately, can search for nearest neighbors to a specific sample id within the index
-Can also perform exact nearest neighbor search within the reduced-dimensional space

**Arguments**

'-x' or '--basename' (required)
The path to the basename of the junction index; it expects to find an Annoy index (basename.annoy.mor), a dictionary that maps each junction to the number of samples in which it's found (basename.freq.mor), a stats file including the total number of samples (basename.stats.mor) in order to run

'-v' or '--verbose' (option flag)
If this option flag is included, morna search will print various status and progress messages during query building and searching to stderr.

'--search-k' (not required, default is 100)
The integer number of nodes to inspect when searching the Annoy index. Larger numbers give more accurate results (better approximation of nearest neighbors) but increase search time.

'-f' or '--format' (not required, default is "sam")
A string that is one of {sam, bed, raw} which defines the format of the query file. Raw format is tab separated lines containing chromosome, start position, end position, and coverage, like this:
chr1    10253   180889  0

'-d' or '--distances' (option flag)
If this option flag is enabled, include cosine distances to each result along with the result.

'-m' or '--metadata' (option flag)
If this option flag is enabled, include associated metadata with results; metadata file must be found at (basename.meta.mor) and must include entries for all sample ids included in results.

'-c' or '--converge' (not required, default is False)
Attempt to converge on a solution with concordance equal to the given argument as a percentage (truncated to 0-100). Not functional, tested, or recommended, will be replaced with exponential back-off algorithm soon(tm).

'-q','--query-id' (not required, default is None)
Integer sample id that is found in the index; if provided, the search is for nearest neighbors to the sample (format argument will be ignored in this case).

'-e' or '--exact' (option flag)
If this option flag is enabled, search for exact nearest neighbor to query within morna's annoy index rather than using annoy's hyperplane division algorithm.

'-r' or '--results' (not required, default is 20)
The number of nearest neighbor results to report

**morna align**

**Usage**
 
**Features**



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

Morna Index



Commands:


'-m' or '--metafile'
A file path to a metadata file with sample index in first column and other keywords in other columns, whitespace delimited.

An example, typical call:
python morna.py index --intropolis ../../Downloads/intropolis_data_dl/sra/intropolis.v1.hg19.tsv.gz -x intropv1/intropv1 --n-trees 10 -v -b 51200 -m  ../../Downloads/intropolis_data_dl/sra/sra_metadata.tsv

Explanation:
index 
    specifies we are using morna index
--intropolis ../../Downloads/intropolis_data_dl/sra/intropolis.v1.hg19.tsv.gz
    specifies that the relative path to the intropolis-like file is ../../Downloads/intropolis_data_dl/sra/intropolis.v1.hg19.tsv.gz
 -x intropv1/intropv1 
    specifies that the relative base file path for each of the output index files will be intropv1/intropv1
--n-trees 10 
    specifies that the annoy index will be built with 10 trees
-v
    specifies verbose mode
-b 51200
    specifies that the buffer size for each sample's junction buffer at 51200 bytes
-m  ../../Downloads/intropolis_data_dl/sra/sra_metadata.tsv
    specifies the relative path to the metadata file is ../../Downloads/intropolis_data_dl/sra/sra_metadata.tsv

Morna Search

Morna Align


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

