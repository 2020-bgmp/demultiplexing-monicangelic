# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 |
| 1294_S1_L008_R2_001.fastq.gz | index1 |
| 1294_S1_L008_R3_001.fastq.gz | index2 |
| 1294_S1_L008_R4_001.fastq.gz | read2 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. ```Your answer here```
    3. ```Your answer here```
3. N base calls
R2  3328051
R3  3328051
code: zcat 1294_S1_L008_R2_001.fastq.gz 1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep "N" | wc -l

## Part 2
1. Define the problem
Overall, we want to identify which reads had indexes that matched, which had index-hopping occur, and which have undetermined base calls. We also need to sort all of the reads into their own files, making sure all the information is contained within them. And we want to keep a count of how many indices match and which don't; before and after quality filtering of index reads. We want to make sure there are an additional four files: two FASTQ files for non-matching index-pairs and two FASTQ files when one or both index reads are unknown or of low quality 
2. Describe output
Properly label and match indices, and tag which reads are reverse compliments
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [4 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
