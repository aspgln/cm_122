import numpy as np
import pandas as pd
import itertools
from collections import defaultdict



# Given 4 viral strains, 10k-length each;
# Find the relative frequencies of each different viral genome

# 1) Do some alignment to each of 4 genomes
# 2) Identify which genome each read comes from
#      (this may be impossible)  -- most of the viral genomes
#      are identical to each other.
# 3) Once we quantify the number of reads for each genome
#    Compute the set of frequencies that most likely generates
#    them.




# Reduce the problem:

## Just consider the SNPs ==> 5 SNPs
## We observe the aggregate frequency of each SNP
## Compute the frequencies of each individual strain.




## Generate some test data

# 4 strains x 5 SNPs

def generate_test_SNPs(coverage=100, set_seed=True, strain_freqs=[.4, .2, .15, .1, .1, .05], n_snps=10):
    ## We're going to assume there is overlap between
    ## Strains for each SNP.

    if set_seed: np.random.seed(96000)

    n_strains = len(strain_freqs)
    assert sum(strain_freqs) == 1

    powerset_iterable = itertools.chain.from_iterable(itertools.combinations(range(n_strains), r)
                                                      for r in range(2, n_strains))

    frozen_powerset = [index_tuple for index_tuple in powerset_iterable]

    while True:
        these_indices = np.random.choice(frozen_powerset, n_snps)
        output_strains = [[] for i in range(n_strains)]
        for index_list in these_indices:
            for i in range(len(output_strains)):
                if i in index_list:
                    output_strains[i].append(1)
                else:
                    output_strains[i].append(0)

        output_strains = [tuple(strain) for strain in output_strains]
        n_unique_strains = len(set(output_strains))
        # test if it satisfies a few conditions
        # if it does, we return
        # otherwise; loop back and try again
        if n_unique_strains == n_strains:
            break

            #     for freq, strain in zip(strain_freqs, output_strains):
            #         print freq, strain

    sampled_strains = np.random.choice(range(n_strains),
                                       size=coverage * n_strains,
                                       p=strain_freqs)

    output_snp_counts = defaultdict(list)

    for strain_index in sampled_strains:
        snp_index = np.random.choice(range(10))
        snp_value = output_strains[strain_index][snp_index]
        output_snp_counts[snp_index].append(snp_value)

    # print sampled_strains

    output_snp_freqs = {k: float(sum(v)) / len(v) for k, v in
                        output_snp_counts.items()}

    #     for k, v in output_snp_counts.items(): print k, v[:10]
    #     for k, v in output_snp_freqs.items(): print k, v
    return output_snp_counts, output_strains




    ## Solve the problem using given test data

    # generate_test_SNPs()


input_snp_counts, input_strains = generate_test_SNPs(coverage=1000)

for k, v in input_snp_counts.items(): print
k, v[:10]


def get_freqs_from_strains_and_counts(input_snp_counts, input_strains):
    #     print input_snp_count
    snp_freqs = []
    for k in range(10):
        raw_snps = input_snp_counts[k]
        snp_freq = float(sum(raw_snps)) / len(raw_snps)
        snp_freqs.append(snp_freq)

    snp_freqs = np.array(snp_freqs)

    strain_matrix = np.array(input_strains)
    print
    strain_matrix

    strain_freqs = np.linalg.lstsq(strain_matrix.T, snp_freqs)
    for freq, strain in zip(strain_freqs[0], input_strains):
        print
        freq, strain
    return


get_freqs_from_strains_and_counts(input_snp_counts, input_strains)