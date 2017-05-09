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
        tmp = []
        for j in range(4):
            tmp.append(strains[j][i])
        if tmp[0] != tmp[1] or tmp[1] != tmp[2] or tmp[2] != tmp[3]:
            tmp = [i, tmp[0], tmp[1], tmp[2], tmp[3]]
            print tmp
            snps.append(tmp)
    return snps

def find_relative_strings(strains, snps):
    strings = []
    for snp in snps:
        a = 0
        for i in range(3):
            if snp[i+1] != snp[i+2]:
                a = i
                break

        strings_for_one_snp = [snp[0], [snp[a+1]], [snp[a+2]]]
        start = snp[0] - 49
        for i in range(50):
            for j in range(2):
                strings_for_one_snp[j+1].append(strains[a+j][start+i : start+i+50])
        print strings_for_one_snp
        strings.append(strings_for_one_snp)
    return strings

def find_snp_freq(reads, strings):
    nsnp = (len(strings))
    hash = defaultdict(list)

    options = []
    for i in range(nsnp):
        option1 = strings[i][1][0]
        option2 = strings[i][2][0]
        options.append([option1, option2])
        for item in strings[i][1]:
            if item in hash.keys():
                hash[item].append(i)
                hash[item].append(option1)
            else:
                hash[item] = [i, option1]
        for item in strings[i][2]:
            if item in hash.keys():
                hash[item].append(i)
                hash[item].append(option2)
            else:
                hash[item] = [i, option2]
    print options

    noption1 = [0 for i in range(nsnp)]
    noption2 = [0 for i in range(nsnp)]
    for read in reads:
        if not (read in hash.keys()):
            continue
        cases = hash[read]
        n_cases = len(cases) / 2
        if n_cases != 1:
            for i in range(n_cases):
                which_snp = hash[read][2*i]
                if cases[2*i+1] == options[which_snp][0]:
                    noption1[which_snp] += float(1) / n_cases
                else:
                    noption2[which_snp] += float(1) / n_cases
    print noption1
    print noption2

    freq = []
    origin = []
    for i in range(nsnp):
        tmp = float(noption1[i]) / (noption1[i] + noption2[i])
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
    relative_strings = find_relative_strings(strains, snps)
    snp_freq, origin = find_snp_freq(reads, relative_strings)
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


            # 0.5 0.3 0.15 0.5



