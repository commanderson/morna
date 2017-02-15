#!/bin/bash
#Requirements: metadata file with tissue identifier and uuid/accession number,
#as well as accession number to sample id mapping file
#To get sample ids of all pancreas files in gtex as pancreas_sample_ids.tsv

#extract accession number from metadata file
for an in $(cat SRP012682.tsv |grep male_pancreas|cut -f 4); 
do 
#extract corresponding sample id from mapping
    id=$(grep $an gtex_sample_index_to_run_accession_num.tsv); 
    echo "$id" ; 
done|cut -f 1 > pancreas_sample_ids.tsv #store sample ids in file

#Assuming index has already been created for gtex with basename gtex_index:
# run morna search on each sample id number 
#and store results in files labeled by sample id in a subdirectory
for line in $(cat ../../Downloads/intropolis_data_dl/gtex/pancreas_sample_ids.tsv); 
do 
    echo $line; python morna.py search -x gtex_index -v -m -q $line -f raw > pancreas_neighbors_results/$line.nns; 
done

#Then, to rename the files to the accession numer instead of the sample id
for line in $(ls -1 pancreas_neighbors_results/); 
do 
    id=${line%????}; 
    nm=$(grep -w $id ../../Downloads/intropolis_data_dl/gtex/gtex_sample_index_to_run_accession_num.tsv|cut -f 2);
    cp pancreas_neighbors_results/$line pancreas_neighbors_final_results/${nm%?}.nns;
done