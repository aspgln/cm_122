from collections import defaultdict, Counter
import cPickle as pickle
from os.path import join, exists, splitext
import time
import os
import sys
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
import numpy as np



def read_reads(reads_fn):
    f = open(reads_fn, 'r')
    first_line = True
    all_reads = []
    for line in f:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        all_reads.append(line.strip())
    return all_reads

def read_strains(strains_fn):
    f = open(strains_fn, 'r')
    first_line = True
    all_strains = []
    for line in f:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        all_strains.append(line.strip())
    return all_strains

def find_snps(strains):
    strain_len = len(strains[0])
    snps = []
    for i in range(strain_len):
        snp = []
        for j in range(4):
            snp.append(strains[j][i])

        # find which strain contains snp
        if snp[0] != snp[1] or snp[0] != snp[2] or snp[0] != snp[3] or snp[1] != snp[2] \
                or snp[1] != snp[3] or snp[2] != snp[3]:

            snp = [i, snp[0], snp[1], snp[2], snp[3]]
            print 'snp :', snp
            snps.append(snp)
    return snps

def find_target_reads(strains, snps):
    target_reads = []
    for snp in snps:
        d = Counter()
        for i in snp[1:5]:
            d[i] += 1
        m = d.most_common(2)
        for j in range (4):
            if snp[j+1] == m[1][0]:
                a = snp[j+1]
                break
        print strains[j-1]
        start = snp[0] - 49
        reads_contains_snp= [snp[0], [snp[a+1], snp[a+2]] ]
        print reads_contains_snp
        for i in range(50):
            for j in range(1,2):
                target_reads[j].append(strains[a+j-1][start+i : start+i+50])
        target_reads.append(reads_contains_snp)
    return target_reads


def calculate_snp_freq(reads, target_reads):
    num_snp = (len(target_reads))
    d = defaultdict(list)

    options = []
    for i in range(num_snp):
        option1 = target_reads[i][1][0]
        option2 = target_reads[i][2][0]
        options.append([option1, option2])
        for j in target_reads[i][1]:
            if j in d.keys():
                d[j].append(i)
                d[j].append(option1)
            else:
                d[j] = [i, option1]
        for k in target_reads[i][2]:
            if k in d.keys():
                d[k].append(i)
                d[k].append(option2)
            else:
                d[k] = [i, option2]
    print options

    num_option1 = [0 for i in range(num_snp)]
    num_option2 = [0 for i in range(num_snp)]
    for read in reads:
        if not (read in d.keys()):
            continue
        cases = d[read]
        n_cases = len(cases) / 2
        if n_cases != 1:
            for i in range(n_cases):
                which_snp = d[read][2*i]
                if cases[2*i+1] == options[which_snp][0]:
                    num_option1[which_snp] += float(1) / n_cases
                else:
                    num_option2[which_snp] += float(1) / n_cases
    print num_option1
    print num_option2

    freq = []
    origin = []
    for i in range(num_snp):
        tmp = float(num_option1[i]) / (num_option1[i] + num_option2[i])
        if tmp <= 0.5:
            freq.append(tmp)
            origin.append(options[i][1])
        else:
            freq.append(1 - tmp)
            origin.append(options[i][0])
    print freq
    print origin

    return freq, origin

def find_strain_matrix(snps, origin):
    matrix = [[], [], [], []]

    for i in range(len(snps)):
        for j in range(4):
            if snps[i][j+1] == origin[i]:   # origin
                matrix[j].append(0)
            else:
                matrix[j].append(1)
    matrix = [tuple(item) for item in matrix]
    print matrix
    return matrix

def solve_equation(strain_matrix, snp_freq):
    strain_matrix = np.array(strain_matrix)
    print strain_matrix
    snp_freq = np.array(snp_freq)
    print snp_freq
    strain_freqs = np.linalg.lstsq(strain_matrix.T, snp_freq)
    strain_freqs = strain_freqs[0]
    print strain_freqs
    return strain_freqs


if __name__ == "__main__":
    input_folder = './hw4_W_1'
    strains_fn = join(input_folder, 'hw4_W_1_strains.txt')
    reads_fn = join(input_folder, 'hw4_W_1_reads.txt')

    strains = read_strains(strains_fn)
    reads = read_reads(reads_fn)

    snps = find_snps(strains)
    target_reads = find_target_reads(strains, snps)
    snp_freq, origin = calculate_snp_freq(reads, target_reads)
    strain_matrix = find_strain_matrix(snps, origin)
    strain_freqs = solve_equation(strain_matrix, snp_freq)

    output_str = ''
    for freq, strain in zip(strain_freqs, strains):
        output_str += str(freq) + ',' + strain + '\n'

    print output_str
    output_fn = join(input_folder, 'ans.txt')
    with(open(output_fn, 'w')) as output_file:
        for line in output_str:
            output_file.write(line + '\n')



