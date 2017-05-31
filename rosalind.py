'''
    Rosalind problems solving
    http://rosalind.info/problems/
'''

#Problem #1 http://rosalind.info/problems/dna/

# inp = input()
# arr = ["A", "C", "G", "T"]
# ret_str = ""
# for i in range(0, 4):
#     ret_str += str(inp.count(arr[i])) + " "
# print(ret_str)

# works

#problem #2

# inp = input()
# ret = inp.replace("T", "U")
# print(ret)

#good

#problem #3

# inp = input()
# rev_inp = inp[::-1]
# def rep(in_str, l1, l2):
#     proc = in_str.replace(l1, "K")
#     proc = proc.replace(l2, "R")
#     proc = proc.replace("K", l2)
#     proc = proc.replace("R", l1)
#     return proc
# ret_str = rep(rev_inp, "A", "T")
# ret_str = rep(ret_str, "C", "G")
# print(ret_str)

#solved

# # Problem #4 Mendel's First Law
#
# # implementation using combinations = simulating actual mating
#
# import itertools as it
#
# k, m, n = [int(x) for x in input().split()]
#
# population = []
# population.extend(k * ['YY'])
# population.extend(m * ['Yy'])
# population.extend(n * ['yy'])
# kids = []
#
# for x in it.combinations(population, 2):
#     kids.append("".join(x[:]))
#
# print(population)
# print(kids)
# rec = kids.count("yyyy") + kids.count("Yyyy") / 2 + kids.count('YyYy') / 4 + kids.count("yyYY") / 2
# print(rec)
# print(len(kids))
# print(1 - rec/len(kids))
# t = n + m + k
# # some stolen code that works
# print(((k*k - k) + 2*(k*m) + 2*(k*n) + (.75*(m*m - m)) + 2*(.5*m*n))/((k + m + n)*(k + m + n -1)))
# # my impelementation via probability that works
# print(1 - (((n/t)*((n-1)/(t-1))) + ((m/t)*((m-1)/(t-1))) *0.25 + ((n/t)*(m/(t-1)))*0.5 + ((m/t)*(n/(t-1)*0.5))))
#
# # done

# Counting Point Mutations
# s = input()
# t = input()
#
# counter = 0
#
# for i in range(len(s)):
#     if s[i] != t[i]:
#         counter += 1
# print(counter)
#
# # done

# Enumerating Gene Orders

# import itertools as it
#
# n = int(input())
# def factorial(n):
#     return 1 if n < 1 else n * factorial(n-1)
#
# print(factorial(n))
#
# input_array = "".join([format(x, '01d') for x in range(1, n+1) ])
#
# arr = list(map(" ".join, it.permutations(input_array)))
# for i in range(len(arr)):
#     print(arr[i])

# done

# Rabbits and Recurrence Relations

# def rabbify(n, k):
#     return 1 if n < 3 else rabbify(n-1, k) + k * rabbify(n-2, k)
#
# print(rabbify(30,3))

# done

# Finding a Motif in DNA

# s = input()
# t = input()
# beg = 0
# return_string = ''
#
# while(s.find(t, beg) != -1):
#     cur = s.find(t, beg)
#     return_string += str(cur+1) + " "
#     beg = cur + 1
#
# print(return_string)
#
# # done

# Computing GC Content

# from Bio import SeqIO
#
# input_file = "rosalind_gc.txt"
#
# def count_CG(data):
#     c_occ = data.count('C')
#     g_occ = data.count('G')
#     percentage = float((g_occ + c_occ) / len(data))
#     return(100 * percentage)
#
# fasta_sequences = list(SeqIO.parse(open(input_file),'fasta'))
#
# max_item = fasta_sequences[1]
#
# for fasta in fasta_sequences:
#     if count_CG(str(fasta.seq)) > count_CG(str(max_item.seq)):
#         max_item = fasta
#
# print(max_item.id, count_CG(max_item.seq), sep="\n")

# done

# Mortal Fibonacci Rabbits

# n, m = input().split()
# n, m = int(n), int(m)
# generations = [1, 1]
#
# def fib(i, j):
#     count = 2
#     while (count < i):
#         if (count < j):
#             generations.append(generations[-2] + generations[-1]) #recurrence relation before rabbits start dying (simply fib seq Fn = Fn-2 + Fn-1)
#         elif (count == j or count == j+1):
#             print ("in base cases for newborns (1st+2nd gen. deaths)") #Base cases for subtracting rabbit deaths (1 death in first 2 death gens)
#             generations.append((generations[-2] + generations[-1]) - 1)#Fn = Fn-2 + Fn-1 - 1
#         else:
#             generations.append((generations[-2] + generations[-1]) - (generations[-(j+1)])) #Our recurrence relation here is Fn-2 + Fn-1 - Fn-(j+1)
#         count += 1
#     return (generations[-1])
#
# print (fib(n, m))
# print ("Here's how the total population looks by generation: \n" + str(generations))

# stolen code

# Translating RNA into Protein

# RNA = input()
#
# def transcript(RNA, step=3):
#     codon_table = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
#                    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
#                    "UAU": "Y", "UAC": "Y", "UAA": "", "UAG": "",
#                    "UGU": "C", "UGC": "C", "UGA": "", "UGG": "W",
#                    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
#                    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
#                    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
#                    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
#                    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
#                    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
#                    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
#                    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
#                    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
#                    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
#                    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
#                    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G", }
#     cur = 0
#     protein = ""
#     while cur < len(RNA):
#         codon = RNA[cur:cur+3]
#         protein += codon_table[codon]
#         cur += step
#     return protein
#
# print(transcript(RNA))

# done

# Calculating Protein Mass

# protein = input()
#
# mass_table= {
#     'A': 71.03711,
#     'C': 103.00919,
#     'D': 115.02694,
#     'E': 129.04259,
#     'F': 147.06841,
#     'G': 57.02146,
#     'H': 137.05891,
#     'I': 113.08406,
#     'K': 128.09496,
#     'L': 113.08406,
#     'M': 131.04049,
#     'N': 114.04293,
#     'P': 97.05276,
#     'Q': 128.05858,
#     'R': 156.10111,
#     'S': 87.03203,
#     'T': 101.04768,
#     'V': 99.06841,
#     'W': 186.07931,
#     'Y': 163.06333,
# }
#
# mass = 0.0
#
# for acid in protein:
#     mass += mass_table[acid]
#
# print(mass)

# done

#
# from Bio import SeqIO
#
# input_file = "rosalind_lcsm.txt"
#
# strings_list = []
#
# fasta_sequences = list(SeqIO.parse(open(input_file),'fasta'))
#
# for item in fasta_sequences:
#     strings_list.append(item.seq)
#
# def long_substr(data):
#     substr = ''
#     if len(data) > 1 and len(data[0]) > 0:
#         for i in range(len(data[0])):
#             for j in range(len(data[0])-i+1):
#                 if j > len(substr) and is_substr(data[0][i:i+j], data):
#                     substr = data[0][i:i+j]
#     return substr
#
# def is_substr(find, data):
#     if len(data) < 1 and len(find) < 1:
#         return False
#     for i in range(len(data)):
#         if find not in data[i]:
#             return False
#     return True
#
# print (long_substr(strings_list))

# stolen code
# https://en.wikipedia.org/wiki/Longest_common_substring_problem

'''Needs Solving'''
# Overlap Graphs

# input_file = "rosalind_grph.txt"
#
# k = 3
#
# def parse_list_from_fasta(file):
#     from Bio import SeqIO
#     strings_list = []
#     fasta_sequences = list(SeqIO.parse(open(input_file), 'fasta'))
#     return fasta_sequences
#
# def check_overlapping(fasta_list, suffix_length):
#     overlapping_graphs = []
#     for outer in fasta_list:
#         for inner in fasta_list:
#             if outer.id != inner.id and str(inner.seq).endswith(str(outer.seq)[:suffix_length]):
#                overlapping_graphs.append((outer.id, inner.id))
#     return(overlapping_graphs)
#
# strings_list = parse_list_from_fasta(input_file)
# for pair in check_overlapping(strings_list, 3):
#     print(pair[1], pair[0])

# # Calculating Expected Offspring
#
# k, l, m, n, o, p = [float(x) for x in input().split()]
#
# prob = 2 * (k + l + m + 0.75*n + 0.5*o + 0*p)
# print(prob)
#
# # done

# Enumerating k-mers Lexicographically

# import itertools as it
#
# f = open('rosalind_lexf.txt', 'r')
#
# array = "".join(list(f.readline().strip().split(" ")))
# k = int(f.readline())
#
# sorted_perms = []
#
# for x in it.product(array, repeat=k):
#     line = ""
#     for i in range(k):
#         line += x[i]
#     sorted_perms.append(line)
#
#
# sorted_perms = sorted(sorted_perms)
#
# for x in sorted_perms:
#     print(x)

# done

# RNA Splicing

# from Bio import SeqIO
#
# input_file = "rosalind_splc.txt"
#
# fasta_sequences = list(SeqIO.parse(open(input_file), 'fasta'))
#
# pre_mrna = fasta_sequences[0]
# rna = str(pre_mrna.seq)
#
# for item in fasta_sequences:
#     start = 1
#     if item.id != pre_mrna.id:
#         # rna = "".join(rna.split(str(item.seq)))
#         rna = rna.replace(str(item.seq), "")
# print(rna)
#
# rna = rna.replace("T", "U")
#
# def transcript(RNA, step=3):
#     codon_table = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
#                    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
#                    "UAU": "Y", "UAC": "Y", "UAA": "", "UAG": "",
#                    "UGU": "C", "UGC": "C", "UGA": "", "UGG": "W",
#                    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
#                    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
#                    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
#                    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
#                    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
#                    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
#                    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
#                    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
#                    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
#                    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
#                    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
#                    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G", }
#     cur = 0
#     protein = ""
#     while cur < len(RNA):
#         codon = RNA[cur:cur+3]
#         protein += codon_table[codon]
#         cur += step
#     return protein
#
# print(transcript(rna))

# done

# Independent Alleles

# def factorial(n):
#     return 1 if n < 1 else n * factorial(n-1)
#
# def count_c(n, k):
#     return(factorial(n) / (factorial(k) * factorial(n-k)))
#
# k, n = [int(x) for x in input().split()]
#
# N = pow(2, k)
#
# acc = 0
#
# for i in range(n, N + 1):
#     combinations = count_c(N, i)
#     success = pow(0.25, i)
#     failure = pow(0.75, (N - i))
#     acc += (combinations * success * failure)
#
# print(acc)
# done

# Transitions and Transversions TRAN

# input_file = "rosalind_tran.txt"
#
# def parse_list_from_fasta(file):
#     from Bio import SeqIO
#     strings_list = []
#     fasta_sequences = list(SeqIO.parse(open(input_file), 'fasta'))
#     return fasta_sequences
#
# fasta_sequences = parse_list_from_fasta(input_file)
#
# transitions = {
#     'A': ['G'],
#     'G': ['A'],
#     'C': ['T'],
#     'T': ['C']
# }
#
# transversions = {
#     'A': ['T', 'C'],
#     'G': ['T', 'C'],
#     'T': ['A', 'G'],
#     'C': ['A', 'G']
# }
#
# def count_tt_ratio(original, mutated):
#     tv_count = 0
#     ts_count = 0
#     for i in range(min(len(original), len(mutated))):
#         if (mutated[i] in transitions.get(original[i])):
#             ts_count += 1
#         elif (mutated[i] in transversions.get(original[i])):
#             tv_count += 1
#     if(tv_count):
#         return(ts_count/tv_count)
#     else:
#         return 0
#
# print(count_tt_ratio(str(fasta_sequences[0].seq), str(fasta_sequences[1].seq)))

# done

# Ordering Strings of Varying Length Lexicographically

# f = open('rosalind_lexv.txt', 'r')
#
# array = f.readline().strip().split(" ")
# k = int(f.readline())
#
# alphabet = "".join(array)
#
# def sorted_with_alphabet(iterable_, alphabet):
#
#     return sorted(iterable_, key=lambda word: [alphabet.index(c) for c in word])
#
# def combine(word, length, base='', res=[]):
#     if length>0:
#         for char in word:
#             res.append(base+char)
#             combine(word, length-1, base+char, res)
#     return res
#
# all_strings = combine(alphabet, length=k)
#
# for i in sorted_with_alphabet(all_strings, alphabet):
#     print(i)

#solved

# Finding a Spliced Motif SSEQ
#
# input_file = "rosalind_sseq.txt"
#
# def parse_list_from_fasta(file):
#     from Bio import SeqIO
#     fasta_sequences = list(SeqIO.parse(open(input_file), 'fasta'))
#     return fasta_sequences
#
# def find_indicies(master, sub):
#     result = []
#     i = 0
#     j = 0
#     while i < len(master) and j < len(sub):
#         if master[i] == sub[j]:
#             result.append(i+1)
#             j += 1
#         i += 1
#
#     return(result)
#
# fasta_seqs= parse_list_from_fasta(input_file)
# indicies = find_indicies(fasta_seqs[0].seq, fasta_seqs[1].seq)
#
#
# print(' '.join([str(x) for x in indicies]))

# done

#Perfect Matchings and RNA Secondary Structures pmch
#
# from math import factorial
#
# def parse_list_from_fasta(input_file):
#     from Bio import SeqIO
#     fasta_sequences = list(SeqIO.parse(open(input_file), 'fasta'))
#     return fasta_sequences
#
# input_file = 'rosalind_pmch.txt'
#
# fasta_sequences = parse_list_from_fasta(input_file)
#
# rna = fasta_sequences[0].seq
#
# perfect_matches = factorial(rna.count('A')) * factorial(rna.count('C'))
#
# print(perfect_matches)

# works

'''needs solving'''
# Locating Restriction Sites
## revp

# def parse_list_from_fasta(input_file):
#     from Bio import SeqIO
#     fasta_sequences = list(SeqIO.parse(open(input_file), 'fasta'))
#     return fasta_sequences
#
# input_file = 'rosalind_revp.txt'
#
# transitions = {
#     'A': 'G',
#     'G': 'A',
#     'C': 'T',
#     'T': 'C'
# }
#

# Pitfalls of Reversing Translation
## mrna

# prot_trans = {}
#
# with open('codons.txt', 'r') as f:
#     for line in f:
#         rna, protein = line.split()
#         prot_trans.setdefault(str(protein), []).append(str(rna))
#
# def count_dict_list(dict_):
#     ret = {}
#     for key in dict_.keys():
#         ret[str(key)] = len(dict_[str(key)])
#     return(ret)
#
# prot_count = count_dict_list(prot_trans)
#
# def count_problem(str, dict_):
#     mult = 1
#     for char in str:
#         mult *= int(dict_[char])
#     mult *= dict_["Stop"]
#     return(mult)
#
# prot_string = str(input())
#
# result = count_problem(prot_string, prot_count) % 1000000
#
# print(result)

# done

# Partial Permutations
## pper
# from math import factorial
# def part_perm(n, k):
#     import math
#     if n == k:
#         return(1)
#     elif k == 1:         # see georg's comment
#         return(n)
#     elif k > n:          # will be executed only if y != 1 and y != x
#         return(0)
#     else:                # will be executed only if y != 1 and y != x and x <= y
#         a = math.factorial(n)
#         c = math.factorial(n-k)  # that appears to be useful to get the correct result
#         res = a / c
#         return(int(res))
#
# n, k = [int(x) for x in input().split()]
#
# print(part_perm(n, k) % 1000000)

# done

# Consensus and Profile
## cons

# def parse_list_from_fasta(input_file):
#     from Bio import SeqIO
#     fasta_sequences = list(SeqIO.parse(open(input_file), 'fasta'))
#     return fasta_sequences
#
# input_file = "rosalind_cons.txt"
#
# fasta_seqs = parse_list_from_fasta(input_file)
#
# str_len = len(fasta_seqs[0].seq)
#
# profile = {
#     'A': [0]*str_len,
#     'C': [0]*str_len,
#     'G': [0]*str_len,
#     'T': [0]*str_len
# }
#
# for item in fasta_seqs:
#     for i in range(str_len):
#         profile[
#                 item[i]
#         ][i] += 1
#
# consensus = ''
#
# for i in range(str_len):
#     max = ('A', profile['A'][i])
#     for key in profile.keys():
#         if profile[key][i] > max[1]:
#             max = (key, profile[key][i])
#     consensus += str(max[0])
#
# print(consensus)
#
# for key in profile.keys():
#     print(key, ': ', sep='', end='')
#     print(*profile[key], sep=' ')

# prob

# import math
#
# fid = open('rosalind_prob.txt', 'r')
# s = fid.readline().strip()
# A = [float(x) for x in fid.readline().split()]
#
# cg = s.count('C') + s.count('G')
# at = len(s) - cg
#
# B = []
# for x in range(len(A)):
#     B.append(math.log10(
#         (A[x] / 2) ** cg *
#                             (.5 - A[x] / 2) ** at)
#     )
#
# for item in B:
#     print(item, end=' ')
'''What?'''

# completing a tree
## tree

# def problem(n, edges):
#     return n - len(edges) - 1
#
#
# dataset = open("rosalind_tree.txt").readlines()
#
# n = int(dataset[0])
#
#
# edges = []
#
# for i in range(1, len(dataset)):
#     edges.append(map(int, dataset[i].split()))
#
# print(edges)
# print(problem(n, edges))
'''What??'''
##inod
# n-2
##sset
# print(2 ** 818 % 1000000)
# read
import re
n, A, B = [re.sub(r"{|}|,", "", x) for x in open('rosalind_seto.txt').readlines()]

# preprocess
n = int(n)
A = set([int(x) for x in str.split(A)])
B = set([int(x) for x in str.split(B)])

universal_set = set(range(1, n+1))

union = set()
intersection = set()
a_sub_b = set.copy(A)
b_sub_a = set.copy(B)
a_comp = set.copy(universal_set)
b_comp = set.copy(universal_set)

for item in B:
    union.add(item)
    b_comp.remove(item)
    if item in A:
        b_sub_a.remove(item)

for item in A:
    union.add(item)
    a_comp.remove(item)
    if item in B:
        intersection.add(item)
        a_sub_b.remove(item)

list_ = [union, intersection, a_sub_b, b_sub_a, a_comp, b_comp]

[print(x) for x in list_]