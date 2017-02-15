TESTS

Find all pancreas samples in gtex and get their 20 nearest neighbors:
First, a morna index with accompanying metadata should be constructed for gtex. 
The command would look something like this:
python morna.py index --intropolis ../../Downloads/intropolis_data_dl/gtex/first_pass_gtex_junctions.tsv.gz -x gtex_index -s 9662 -m ../../Downloads/intropolis_data_dl/gtex/gt_metadata.ssv 

With this index, as well as the metadata files, we can extract the sample ids of 
pancreas samples and query our morna index for their nearest neighbors.
The process is outline in all_gtex_pancreas_nns.bash, and is as follows,
noting that the bash script assumes it is being run from the morna directory:

To get a file listing the sample ids:

#extract accession number from metadata file
for an in $(cat ../../Downloads/intropolis_data_dl/gtex/SRP012682.tsv |grep male_pancreas|cut -f 4); 
do 
#extract corresponding sample id from mapping
    id=$(grep $an gtex_sample_index_to_run_accession_num.tsv); 
    echo "$id" ; 
done|cut -f 1 > pancreas_sample_ids.tsv #store sample ids in file

To run a morna search for each sample id and store the result:

#Assuming index has already been created for gtex with basename gtex_index:
# run morna search on each sample id number 
#and store results in files labeled by sample id in a subdirectory
for line in $(cat ../../Downloads/intropolis_data_dl/gtex/pancreas_sample_ids.tsv); 
do 
    echo $line; python morna.py search -x gtex_index -v -m -q $line -f raw > pancreas_neighbors_results/$line.nns; 
done

Lastly, to rename the files to the accession numer instead of the sample id:

for line in $(ls -1 pancreas_neighbors_results/); 
do 
    id=${line%????}; 
    nm=$(grep -w $id ../../Downloads/intropolis_data_dl/gtex/gtex_sample_index_to_run_accession_num.tsv|cut -f 2);
    cp pancreas_neighbors_results/$line pancreas_neighbors_final_results/${nm%?}.nns;
done