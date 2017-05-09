from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
from BIOINFO_M260B.helpers import read_reads


def read_assembly_reads(read_fn):
    reads = read_reads(read_fn)
    output_reads = [_[0] for _ in reads]
    # Only taking one end of the read works okay, but
    # this is an obvious area for improvement.
    return output_reads


def simple_de_bruijn(sequence_reads, k):
    """
    Creates A simple DeBruijn Graph with nodes that correspond to k-mers of size k.
    :param sequence_reads: A list of reads from the genome
    :param k: The length of the k-mers that are used as nodes of the DeBruijn graph
    :return: A DeBruijn graph where the keys are k-mers and the values are the set
                of k-mers that
    :return: A reversed DeBruijn graph where the keys are k-mers next, and the values
    are of kmers pre

    """
    de_bruijn_counter = defaultdict(Counter)    # default dict of type counter
                                                # to filter out bad data from our thing
    de_bruijn_counter_rev = defaultdict(Counter)
    # You may also want to check the in-degree and out-degree of each node
    # to help you find the beginnning and end of the sequence.
    for read in sequence_reads:
        # Cut the read into k-mers
        kmers = [read[i: i + k] for i in range(len(read) - k)]
        # associate
        # each kmer is a counter
        # store the next kmer as its key

        for i in range(len(kmers) - 1):
            pvs_kmer = kmers[i]
            next_kmer = kmers[i + 1]
            de_bruijn_counter[pvs_kmer].update([next_kmer])
            de_bruijn_counter_rev[next_kmer].update([pvs_kmer])

            # counter is a module served like a dictionary, only designed to increment by 1
            # every time you see it
            # it's going to set counter[next_kmer] += 1
            # e.g. dbj_counter_['CATG'] = {'ATGG': 20, 'ATCG' : 1}
            # if in the next read is CATGG
            # then dbj_counter_['CATG'] = {'ATGG': 21 'ATCG': 1}


    # This line removes the nodes from the DeBruijn Graph that we have not seen enough.
    de_bruijn_graph = {key: {val for val in de_bruijn_counter[key] if de_bruijn_counter[key][val] > 10}
                       for key in de_bruijn_counter}
    # de_bruijn_graph is a dictionary of previous_kmers
    # previous_kmers is a Counter of next_kmers

    # This line removes the empty nodes from the DeBruijn graph
    de_bruijn_graph = {key: de_bruijn_graph[key] for key in de_bruijn_graph if de_bruijn_graph[key]}




    de_bruijn_graph_rev = {key: {val for val in de_bruijn_counter_rev[key] if de_bruijn_counter_rev[key][val] > 10}
                       for key in de_bruijn_counter_rev}
    de_bruijn_graph_rev = {key: de_bruijn_graph_rev[key] for key in de_bruijn_graph_rev if de_bruijn_graph_rev[key]}



    return de_bruijn_graph, de_bruijn_graph_rev
    #
    # pre_kmer: {next_kmer if the deB counter [previous kmer][next_kmer] > 1}}
    # if we something twcie, we are goint to include in out dB graph
    # improvement: if reads is 30, set to > 20?
    # also, to have a idea of how oftern these sequecne show up
    # this only show the single edge bettwen the provious kmer and the next kmer
    # we can also keep this in a dictionalry of previous kmer to the next kmer , and the number
    # of times that next_kmer shows up
    # => count the expdcted number of times you would see the next kmer
    # for repetitive sequecne
    #
    #
    #
    # '''
    # also we might want to see the sequence that we chopped out at the very end, and recover some
    # of them as well
    #
    # '''





def de_bruijn_reassemble(de_bruijn_graph, de_bruijn_graph_rev):
    # takes in a dict of sets, dict of dict or dict of counter
    """
    Traverses the DeBruijn Graph created by simple_de_bruijn and
    returns contigs that come from it.
    :param de_bruijn_graph: A De Bruijn Graph
    :return: a list of the
    """

    # find a better place to start the deB graph : has more output than inputs? by simple search

    '''
    implement backtracking!!!
    '''
    assembled_strings = []
    while True:
        n_values = sum([len(de_bruijn_graph[k]) for k in de_bruijn_graph])
        if n_values == 0:
            break

        good_starts = [k for k in de_bruijn_graph if de_bruijn_graph[k]]
        print type(good_starts)





        # print 'good start = ', good_starts
        # You may want to find a better start
        # position by looking at in and out-degrees,
        # but this will work okay.
        current_point = good_starts[0]
        assembled_string = current_point
        while True:
            try:
                next_values = de_bruijn_graph[current_point] # set next kmer to previous kmer
                next_edge = next_values.pop()# pop is useful to set in py
                # if pop, the size of the counter will decrease by 1
                # what ever value popped out get stored and returned

                # could have fail because the still an outgoing edge, but already being pop out

                assembled_string += next_edge[-1]
                # add the last character of the next edge
                de_bruijn_graph[current_point] = next_values
                # editing dbj graph as we go through = next_value
                current_point = next_edge
                #
            except KeyError:
                # if run into a keyerror, yiel the assemble string, and start over from
                # another random key from out dbj graph and start over
                assembled_strings.append(assembled_string)
                break
    return assembled_strings


if __name__ == "__main__":
    chr_name = 'hw3all_A_3_chr_1'
    input_folder = '../data/{}'.format(chr_name)
    reads_fn = join(input_folder, 'reads_{}.txt'.format(chr_name))
    reads = read_assembly_reads(reads_fn)
    db_graph, db_graph_rev = simple_de_bruijn(reads, 25)

    # for k in db_graph.keys():
    #     print k, db_graph[k]
    print len(db_graph.keys())

    output = de_bruijn_reassemble(db_graph, db_graph_rev)
    output_fn_end = 'assembled_{}.txt'.format(chr_name)
    output_fn = join(input_folder, output_fn_end)
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + chr_name + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(output))
