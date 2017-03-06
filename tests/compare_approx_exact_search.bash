#Assumes run from morna folder containing morna.py
#Make morna index for each n_trees you want to test
for nt in {10,20,50,100,200};
do 
    echo $nt trees;
    python morna.py index --intropolis ../../Downloads/intropolis_data_dl/gtex/first_pass_gtex_junctions.tsv.gz -x gtex_index_nt$nt -s 9662 -v --n-trees $nt -m ../../Downloads/intropolis_data_dl/gtex/gt_metadata.ssv;
done

#Then, use a sample query to test the similarity of results between exact and approximate nearest neighbor search for each such index at each desired search_k.
for nt in {10,20,50,100,200};
do
    echo $nt trees;
    for sk in {10,20,50,100,200};
    do
        echo $sk search_k;
        cat ../../Downloads/intropolis_data_dl/gtex/query_beds/junctions.SRR1076268_SRS524526_SRX408563_male_stomach.bed | python morna.py search -x gtex_index_nt$nt -f bed --search-k $sk -m --exact > parameter_test_results/$nt-$sk-ex.res;
        cat ../../Downloads/intropolis_data_dl/gtex/query_beds/junctions.SRR1076268_SRS524526_SRX408563_male_stomach.bed | python morna.py search -x gtex_index_nt$nt -f bed --search-k $sk -m > parameter_test_results/$nt-$sk-ap.res;
    done;
done
