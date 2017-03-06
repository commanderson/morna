TESTS

1.Find all pancreas samples in gtex and get their 20 nearest neighbors:

First, a morna index with accompanying metadata should be constructed for gtex. 
The command would look something like this:
python morna.py index --intropolis ../../Downloads/intropolis_data_dl/gtex/first_pass_gtex_junctions.tsv.gz -x gtex_index -s 9662 -m ../../Downloads/intropolis_data_dl/gtex/gt_metadata.ssv 

With this index, as well as the metadata files, we can extract the sample ids of 
pancreas samples and query our morna index for their nearest neighbors.
The process is outlined in all_gtex_pancreas_nns.bash, and is as follows,
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



2.Test the accuracy and precision of the norm estimator e2 produced by test_norm_estimator.py:

run the following bash script (located in test_estimator_params.bash):

echo -e -n "Num. Blocks (s)\tNum. Loops (l)\tReal value\t" > estimator_test_results.tsv
echo -e -n "Mean e2\te2 variance\tMean difference\t" >> estimator_test_results.tsv
echo -e "Mean e2/rv ratio\tX-norm\tMean y-norm\tY-norm variance" >> estimator_test_results.tsv

#outer loop is on a list of desired s values,
#inner loop is on a list of desired random redraw loops, or l values
for s in {1,2,5,10};
do	
	echo s of $s...
	for l in {10,100,500};
	do
		echo ...l of $l;
		python test_norm_estimator.py -s $s -l $l>>estimator_test_results.tsv;
	done;
done

3. Test the degree to which the approximate nearest neighbor search carried out by annoy replicates the results of the exact nearest neighbor search of all vectors in the index (carried out by morna search with the --exact parameter) at various levels of num_trees (an indexing parameter) and search_k (a search parameter.The process is outlined in compare_approx_exact_search.bash, and is as follows,
noting that the bash script assumes it is being run from the morna directory:

First, a morna index with accompanying metadata should be constructed for gtex, for each different value of n_trees you want to test. 
The command would look something like this:

for nt in {10,20,50,100,200};
do 
    echo $nt trees;
    python morna.py index --intropolis ../../Downloads/intropolis_data_dl/gtex/first_pass_gtex_junctions.tsv.gz -x gtex_index_nt$nt -s 9662 --n-trees $nt -m ../../Downloads/intropolis_data_dl/gtex/gt_metadata.ssv;
done

Then, use a sample query to test the similarity of results between exact and approximate nearest neighbor search for each such index at each desired search_k. The command would look like this:
for nt in {10,20,50,100,200};
do
    echo $nt trees;
    for sk in {10,20,50,100,200};
    do
        echo $sk search_k;
        cat ../../Downloads/intropolis_data_dl/gtex/query_beds/junctions.SRR1076268_SRS524526_SRX408563_male_stomach.bed | python morna.py search -x gtex_index_nt$nt -f bed --search-k $sk --exact > parameter_test_results/$nt-$sk-ex.res;
        cat ../../Downloads/intropolis_data_dl/gtex/query_beds/junctions.SRR1076268_SRS524526_SRX408563_male_stomach.bed | python morna.py search -x gtex_index_nt$nt -f bed --search-k $sk > parameter_test_results/$nt-$sk-ap.res;
    done
done



