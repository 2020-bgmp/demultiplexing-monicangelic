#!/usr/bin/env python

import argparse
import gzip
import itertools
import re

def get_args():
    '''Get names of input files, qscores, and lengths'''
    parser = argparse.ArgumentParser(description = "analyzes and returns the mean quality scores of all four files and outputs results to graph")
    parser.add_argument("-i", "--indexfile", type=str, help="File of 24 indexes", required = True)
    #parser.add_argument("-q", "--qscore", type=int, help="Qscore cutoff", required = True)
    parser.add_argument("-r1", "--read1", type=str, help="Read 1 fastq file", required = True)
    parser.add_argument("-r2", "--read2", type=str, help="Index 1 fastq file", required = True)
    parser.add_argument("-r3", "--read3", type=str, help="Index 2 fastq file", required = True)
    parser.add_argument("-r4", "--read4", type=str, help="Read 2 fastq file", required = True)

    return parser.parse_args()
args = get_args()

#Setting global variables
index_file = args.indexfile
#q_cutoff = args.qscore
read1fq = args.read1
index1fq = args.read2
index2fq = args.read3
read2fq = args.read4


#Dictionary set up for all index information. Sequences = keys, values = remaining information

with open(index_file, 'r') as fh:
    index_dict = {}
    LN = 0
    for line in fh:
        line = line.strip("\n")
        line = line.split("\t")
        LN+=1
        if LN > 1:
            pos, index = line[3], line[4]
            index_dict[pos] = index

#Function1
def rev_comp(sequence):
    '''To get the reverse compliment of a sequence'''
    revseq = ""
    for i in sequence:
        if i == 'A':
            revseq += 'T'
        if i == 'T':
            revseq += 'A'
        if i == 'C':
            revseq += 'G'
        if i == 'G':
            revseq += 'C'
        if i == 'N':
            revseq += 'N'
    revcomp = revseq[::-1]

    return revcomp

#Function2
#def index_file_dict(index_dict):
 #   in_file_dict = {}
  #  for name in index_dict:
   #     seq = index_dict[name]
    #    index_file_dict[seq] = (str(name)+"_fwd.fastq",str(name)+"_rev.fastq")
    #return index_file_dict

#in_file_dict = index_file_dict(index_dict)

#Function3
def convert_phred(letter):
    """Converts a single character into a phred score"""
    score =ord (letter)-33
    return score

#Function3
def qmean_score(line):
    sum = 0
    for pos in line:
        score = convert_phred(pos)
        sum += score
    return(sum/len(line))


#Create another dictionary of all permutations of index combos
perm_dict = {}

for i in itertools.product(index_dict.keys(), repeat = 2):
    perm_dict[i] = 0
#print(perm_dict)

#Create dictionary that has matched counts
match_counts = {}

for sequence in index_dict.values():
    match_counts[sequence] = 0


#Create and open output files

#Make dictionary for R1 and R2 files, but not for the other files
#Function3
def created_matched_files(index_dict):
    fdict_r1 = {}
    fdict_r4 = {}
    for index in index_dict:
        seq = index_dict[index]
        fdict_r1[seq] =open(str(index)+"_matched_r1.fastq", 'w')
        fdict_r4[seq] =open(str(index)+"_matched_r4.fastq", 'w')
        #'num_%s_matched_r4.fastq.gz' % num, 'w')
    return fdict_r1, fdict_r4

fdict_r1, fdict_r4 = created_matched_files(index_dict)


hop_R1 =open('indexhopped_r1.fastq', 'w')
hop_R4 =open('indexhopped_r4.fastq', 'w')
unk_R1 =open('unknown_r1.fastq', 'w')
unk_R4 =open('unknown_r4.fastq', 'w')

#Open all 4 fastq input files 
R1_File = gzip.open(read1fq, 'rt')
R2_File = gzip.open(index1fq, 'rt')
R3_File = gzip.open(index2fq, 'rt')
R4_File = gzip.open(read2fq, 'rt')

#Create lists to hold a record from each input file
r1_list = []
r2_list = []
r3_list = []
r4_list = []

#Create counters for output files
LN = 0
unk_count = 0
hop_count = 0
total_mcount = 0
total_r = 0

#Start the main algorithm

for liner1, liner2, liner3, liner4 in zip(R1_File,R2_File,R3_File,R4_File):
    LN+=1
    r1_list.append(liner1.rstrip())
    r2_list.append(liner2.rstrip())
    r3_list.append(liner3.rstrip())
    r4_list.append(liner4.rstrip())
    #Get reverse compliment of r2
    #in1_seq = r2_list[1]
    #revcom_in2_seq = rev_comp(r3_list[1])
    if len(r1_list)==4:
        r1_list[0] = r1_list[0] + ' ' + r2_list[1] + '-' + r3_list[1]
        r4_list[0] = r4_list[0] + ' ' + r2_list[1] + '-' + r3_list[1]
        total_r+=1
        if 'N' in r2_list[1] or 'N' in r3_list[1] or qmean_score(r2_list[3]) < 30 or qmean_score(r3_list[3]) < 30 or rev_comp(r3_list[1]) not in index_dict or r2_list[1] not in index_dict:
            unk_R1.write('\n'.join(r1_list)+'\n')
            unk_R4.write('\n'.join(r4_list)+'\n')
            unk_count+=1
        elif r2_list[1] == rev_comp(r3_list):
            variable = index_dict[r2_list[1]]
            fdict_r1[variable].write('\n'.join(r1_list)+'\n')
            fdict_r4[variable].write('\n'.join(r4_list)+'\n')
            total_mcount += 1
            if variable in match_counts:
                match_counts[variable] += 1
        else:
            hop_R1.write('\n'.join(r1_list)+'\n')
            hop_R4.write('\n'.join(r4_list)+'\n')
            hop_count +=1
            if (r2_list[1], rev_comp(r3_list[1])) in perm_dict:
                perm_dict[r2_list[1], rev_comp(r3_list[1])]+=1
        
        #Clear the records
        r1_list = []
        r2_list = []
        r3_list = []
        r4_list = []

#Close input and output files
for pos in index_dict.values():
    fdict_r1[pos].close()
    fdict_r4[pos].close()

unk_R1.close()
unk_R4.close()
hop_R1.close()
hop_R4.close()
R1_File.close()
R2_File.close()
R3_File.close()
R4_File.close()

print('All files closed')

#Writing out stat summary to a text file
with open('demultiplex_summary.txt','w') as fh:
    fh.write("Total reads: " + str(total_r) + '\n')
    fh.write("Number of matched reads: " + str(total_mcount) + '\n')
    fh.write("Number of index hopped reads: " + str(hop_count) + '\n')
    fh.write("Number of unknown reads: " + str(unk_count) + '\n')

#Writing out table of index permutations
with open('indexperm_counts.tsv', 'w') as fh1:
    fh1.write('index pair combination' + '\t' + 'number of indexes in category' + '\n')
    for pair, count in perm_dict.items():
        fh1.write(str(pair) + '\t' + str(count) + '\t' + str(round(count/hop_count*100, 2)) + '%''\n')