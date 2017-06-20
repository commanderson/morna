#!/usr/bin/env python
"""
alignment_comparison

compare a list of .bam alignment files with a master .bam alignment
to determine what percentage of reads in each listed file have an identical
highest-scoring mapping (coordinate + cigar string) in the master file

Note the concept of fractional mappings - if master file has 4 equally high scoring mappings, and 3 of these are present in a compared file, then .75 of an identical mapping is present for that read in that file.

"""



def go(master, subfiles, detailed_stats = True, verbose = False):
    """ Runs the actual comparison; assumes that all subfiles and the masterfile 
        have been sorted by read name and FLAG in alphanueric order - that is,
        if any subfile has e.g. reads 
        ERR922713.160930574_HS3_322:7:2206:10273:60668  147
        ERR922713.160930574_HS3_322:7:2206:10273:60668  99
        they must be in that order, and they must also be in that 
        order in the master file
        use the same sorting operation for both types of files!
        
    """
    
    process1 = subprocess.Popen(["samtools","view","-F","4",
                                 master], stdout=subprocess.PIPE)
    #some high-level variables needed before looping                             
    sub_iterators = []
    current_read_groups = []
    aligned_reads_in_master = 0
    
    aligned_reads_in_subfiles = [0 for _ in subfiles]
    precision_scores = [0.0 for _ in subfiles]
    recall_scores = [0.0 for _ in subfiles]
    sub_finished = False #flag to shortcut if all subalignments are finished
    
    for filename in subfiles:
        process_sub = subprocess.Popen(["samtools","view","-F","4",
                                     filename], stdout=subprocess.PIPE)
        sub_iterators.append(groupby(iter(process_sub.stdout.readline, ''),
                                lambda x: x.split()[0] + "_" + x.split()[1]))

    for i,iterator in enumerate(sub_iterators):
        try:
            current_read_groups.append(next(iterator))
            aligned_reads_in_subfiles[i] += 1
        except StopIteration:
            sys.stdout.write("Can't read any aligned reads from " + subfiles[i]
                                                + "\n Exiting...")
            exit()
    
    
    #We must loop through every read name group in master file
    for name, group in groupby(iter(process1.stdout.readline, ''),
                                lambda x: "_".join(x.split()[:2])):
        aligned_reads_in_master += 1
                
        if verbose:
            sys.stdout.write("Considering master read " + name + "\n")
    
        #roll the current read group for each subfile forward 
        #until it isn't behind the masterfile
        for i, subgroup in enumerate(current_read_groups):
            while ((current_read_groups[i][0] < name) 
                    and current_read_groups[i][0] != None): 
                #print("in file " +str(i) + " subgroup name " 
                #+ current_read_groups[i][0] + " sorts before " + name 
                #+ " so continuing.")
                try:
                    current_read_groups[i] = next(sub_iterators[i])
                    aligned_reads_in_subfiles[i] += 1
                except StopIteration:
                    #If our iterators have reached their end, note this with
                    #a tuple of None values
                    current_read_groups[i] = (None,None)
                    #If all iterators are finished, flag to shortcut 
                    if [(None,None) for pair 
                            in current_read_groups] ==  current_read_groups:
                        sub_finished = True
                
        
        #print("For master group " + name + ", stabilised with subgroup names ")
        #for i, subgroup in enumerate(current_read_groups):
        #    print(subgroup[0])
            
        if sub_finished:
            continue
            
        current_subfile_names = [pair[0] for pair in current_read_groups]
        #if the group is not present in any subfile we can skip it
        if name in current_subfile_names:
            #once we confirm that one or more subfile read names matches the 
            #current master file read name, we get to work finding the highest
            #scoring alignments for the master file
            max_alignment_score=-99999999999
            alignments_list = []
            for alignment in group:
                fields = alignment.split()
                score = int(fields[11].strip("AS:i"))
                #Is it a max scoring alignment?
                #If it's the new max, make a new list with just it
                if score > max_alignment_score:
                    max_alignment_score = score
                    alignments_list = [(fields[0],fields[3],fields[5],)]
                #If it ties the max, add it to the list of max scorers
                elif score == max_alignment_score:
                    alignments_list.append((fields[0],fields[3],fields[5],))
                #else do nothing
        
            #if verbose:
                #sys.stdout.write("----------------------mainalign-------"
                #                                            + "------------\n")
                #sys.stdout.write("Alignments in this group:\n")
                #for entry in alignments_list:
                #    sys.stdout.write(str(entry) + "\n")
                #sys.stdout.write("Max score: " + str(max_alignment_score) 
                #                                                    + "\n")
            #Then, we must find the max scored alignments for the subfile(s)
            #with matching read name
            for i,file in enumerate(subfiles):
                #if the current group from that file has the same read name
                if current_subfile_names[i] == name:
                    #aligned_reads_in_subfiles[i] += 1
                    #then we need to look at the alignments in that group
                    #to see how many from our list are matched

                    matching = 0.0
                    max_sub_alignment_score=-99999999999
                    sub_alignments_list = []
                    for sub_alignment in current_read_groups[i][1]:
                        sub_fields = sub_alignment.split()
                        sub_score = int(sub_fields[11].strip("AS:i"))
                        #Is it a max scoring sub-alignment?
                        #If it's the new max, make a new list with just it
                        if sub_score > max_sub_alignment_score:
                            max_sub_alignment_score = sub_score
                            sub_alignments_list = [(sub_fields[0],
                                                sub_fields[3],sub_fields[5],)]
                        #If it ties the max, add it to the list of max scorers
                        elif sub_score == max_sub_alignment_score:
                            sub_alignments_list.append((sub_fields[0],
                                                sub_fields[3],sub_fields[5],))
                        #else do nothing
                
                
                 #   if verbose:
                 #       sys.stdout.write("^^^^^^^^^^^^^subalign " 
                 #                           + "^^^^^^^^^^^^^^^^^^^^^^^^^^^\n")
                 #       
                 #       sys.stdout.write("Alignments in this group:\n")
                 #       for entry in alignments_list:
                 #           sys.stdout.write(str(entry) + "\n")
                 #       sys.stdout.write("Max score: " + 
                 #                       str(max_alignment_score) + "\n")
                    
                    #Now, we take the size of the intersection of the master 
                    #file's list of highest-scoring alignments and 
                    #this sub alignment's list
                    for sub_alignment in sub_alignments_list:
                        if sub_alignment in alignments_list:
                            if verbose:
                                sys.stdout.write("In file " + file + ",\n")
                                sys.stdout.write("Subfile alignment " + str(sub_alignment) + " matches master alignment.\n")
                            matching += 1.0
                    #precision is what proportion of the highest-scoring 
                    #alignments in the sub file are in the master file
                    precision_scores[i] += matching / len(sub_alignments_list)
                    
                    #recall is what proportion of highest-scoring alignments in 
                    #the master file are in the sub file
                    recall_scores[i] += matching / len(alignments_list)
    
    #unless we spun out all of our sub files' iterators already, we need to 
    #finish them to get accurate count of aligned reads in subfiles  
    if not sub_finished:
        for i, subgroup in enumerate(current_read_groups):
            while (current_read_groups[i][0] != None): 
                #print("in file " +str(i) + " subgroup name " 
                #+ current_read_groups[i][0] + " sorts before " + name 
                #+ " so continuing.")
                try:
                    current_read_groups[i] = next(sub_iterators[i])
                    aligned_reads_in_subfiles[i] += 1
                except StopIteration:
                    #If our iterators have reached their end, note this with
                    #a tuple of None values
                    current_read_groups[i] = (None,None)
                    #If all iterators are finished, flag to shortcut 
                    if [(None,None) for pair 
                            in current_read_groups] ==  current_read_groups:
                        sub_finished = True
                
    
    sys.stdout.write("Of " + str(aligned_reads_in_master) + " aligned reads in " 
            + "the master file:\n")

    #for i,score in enumerate(precision_scores):
    #    precision_scores[i] /= aligned_reads_in_subfiles[i]
    #    recall_scores[i] /= aligned_reads_in_master
    

            
    for i,file in enumerate(subfiles):
        sys.stdout.write("In subfile " + file + ":\n")
        if detailed_stats:
            sys.stdout.write("*" + str(aligned_reads_in_subfiles[i]) 
                                                        + " aligned reads\n")
            sys.stdout.write("*Precision score: " + str(precision_scores[i]) 
                                                                        + "\n")
        sys.stdout.write("Precision of: " 
            + str(precision_scores[i] / aligned_reads_in_subfiles[i]) + "\n")
            
        if verbose:
            sys.stdout.write("*Recall score: " + str(recall_scores[i]) + "\n")
        sys.stdout.write("Recall of: " 
                    + str(recall_scores[i] / aligned_reads_in_master) + "\n")
                    
    

if __name__ == '__main__':
    import argparse
    import subprocess
    import sys
    import time
    from itertools import groupby

    _help_intro = \
    """ alignment_comparison determines what portion of reads in each of a list 
        of bam files have an identical mapping (coordinate + cigar string) 
        present in a master bam file. 
        NOTE: Requires that reads be sorted by read name (aka the QNAME field)
        and FLAG in alphanumeric order
        use sort -n for this effect
        to get files sorted this way (MUST USE SAME SORT ON MASTER + ALL SUBS)
    """

    parser = argparse.ArgumentParser(description=_help_intro)
    
    parser.add_argument('-m', '--master', metavar='<file>', type=str,
            required=True,
            help=('path to .bam alignment file which serves as master')
        )
    parser.add_argument('-s', '--subfiles', metavar='<m1>',
            type=str, nargs='+',
            required=False,
            help='paths to .bam alignment files to evaluate'
        )
    parser.add_argument('-d', '--detailed-stats', action='store_const',
            const=True,
            default=False,
            help='include more detailed output statistics'
        )
    parser.add_argument('-v', '--verbose', action='store_const',
            const=True,
            default=False,
            help='be talkative'
        )
    parser.add_argument('-t', '--test', action='store_const',
            const=True,
            default=False,
            help='run tests'
        )
        
    args = parser.parse_args()
    if args.test:
        #run tests here
        #test for aligned_read counts : 
        #for i in $(ls of test sams); 
        #do samtools view -F 4 $i|grep -v "^@" |cut -f 1,2|uniq|wc -l; 
        #done
        print()
        
    start_time = time.time()
    go(master = args.master, subfiles = args.subfiles, 
        detailed_stats = args.detailed_stats, verbose = args.verbose)
    sys.stderr.write("DONE with alignment_comparison.py;\ntime=%0.3f s" 
                                        % (time.time() - start_time)+ "\n")
    