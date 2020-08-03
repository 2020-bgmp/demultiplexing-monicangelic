#!/usr/bin/env python

import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description = "analyzes and returns the mean quality scores of all four files and outputs results to graph")
    parser.add_argument("-f", "--file", help="The fastq file being used", required = True)
    parser.add_argument("-l", "--length", type=int, help="The length of sequence", required = True)
    parser.add_argument("-o", "--output", help="The name of the output graph", required = True)

    return parser.parse_args()
args = get_args()

inputfile = args.file
sequence_length = args.length
outputfile = args.output

def init_list(array, value=0.0):
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    
    for i in range(sequence_length):
        array.append(value)
    return array


def convert_phred(letter):
    """Converts a single character into a phred score"""
    
    score =ord(letter)-33
    return score

qscore_sum = []
qscore_sum = init_list(qscore_sum)
print(qscore_sum)

def pop_sumlist(inputfile):
    '''Opens fastq file and takes qscore line which calculates the sum 
        of the quality scores at each position in the seq'''

    line_c = 0

    with open(inputfile, "r") as fh:
        for line in fh:
            line_c +=1
            if line_c%100000000 == 0:
                print("completed: ", line_c)
            if line_c%4 == 0:
                letter_c=0
                for letter in line.rstrip("\n"):
                    score = convert_phred(letter)
                    qscore_sum[letter_c]+=score
                    letter_c += 1
    return qscore_sum, letter_c, line_c
qscore_sum, letter_c, line_c = pop_sumlist(inputfile)
print(letter_c)


def calculate_avg(inputfile):
    '''Takes qscores and gives the average for each position'''
        
    average_scores = []
    for running_sum in qscore_sum:
            record = line_c/4
            average = running_sum/record
            average_scores.append(average)
    
    return average_scores
average_scores = calculate_avg(inputfile)

import matplotlib.pyplot as plt

plt.bar(range(0,101),average_scores)
plt.xlabel("Base Pair")
plt.ylabel("Mean Quality Score")
plt.title("Distribution of Mean Quality Score")
plt.savefig(outputfile)