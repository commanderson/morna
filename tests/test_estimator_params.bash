#!/bin/bash

#First, construct header as though on a single run of estimator_test_results.tsv with -i.
#Assumes it's being run from tests directory, where test_norm_estimator.py is located.

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
