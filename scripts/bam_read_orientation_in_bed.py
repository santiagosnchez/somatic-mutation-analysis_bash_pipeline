#!/bin/env python3

import pysam
import sys
import math

def read_bed(file):
    bed = []
    with open(file) as f:
        for line in f:
            bed.append(line.rstrip().split("\t"))
    return bed

def make_counter(keys):
    counter = {keys[0]:0, keys[1]:0}
    return counter

def count_reads(read, counter):
    if read.is_proper_pair:
        if not read.is_reverse:
            if read.is_read1 and read.pos < read.mpos:
                counter['f1r2'] += 1
            else:
                counter['f2r1'] += 1
        else:
            if read.is_read2 and read.pos > read.mpos:
                counter['f1r2'] += 1
            else:
                counter['f2r1'] += 1
    return counter

message="""

Counts read in BAM that are F1R2 or F2R1 within intervals (BED)

read_orientation.in_bed.py <bam_file> <bed_file> [ <sample_name> ]

examples:
read_orientation.in_bed.py alignment.bam intervals.bed > read_orientation.tsv
read_orientation.in_bed.py alignment.bam intervals.bed sample1

Columns are:
chrom   start    end    gene     F1R2     F2R1    log2(F1R2+1/F2R1+1)

"""

# progess args
if not (len(sys.argv) == 3 or len(sys.argv) == 4):
    print(message)
    sys.exit(1)
elif any([ a == '-h' or a == '--help' for a in sys.argv ]):
    print(message)
    sys.exit(1)

# get file names
bamfile = sys.argv[1]
bedfile = sys.argv[2]
sample = None
if len(sys.argv) == 4:
    sample  = sys.argv[3]
    out_all = open(sample + ".read_orientation.target_count.tsv", "w")
    out_sum = open(sample + ".read_orientation.summary.tsv", "w")
    # for summary stats
    log2ratio_all = []
    sum_counts = (0,0)

# read data objects
bam = pysam.AlignmentFile(bamfile)
bed = read_bed(bedfile)

# determine col 4
col1 = 'range'
if len(bed[0]) == 4:
    col1 = 'gene'

# print header
if sample:
    print('chrom','start','end','target','F1R2','F2R1','log2_F1R2_over_F2R1',sep="\t", file=out_all)
else:
    print('chrom','start','end','target','F1R2','F2R1','log2_F1R2_over_F2R1',sep="\t")

# parse intervals
for interval in bed:
    ch = interval[0]
    st = interval[1]
    en = interval[2]
    if col1 == 'gene':
        ge = interval[3]
    else:
        # make an interval name
        ge = ch + ":" + st + "-" + en
    # make dict counter
    counter = make_counter(['f1r2','f2r1'])
    # count
    for read in bam.fetch(ch, int(st), int(en)):
        counter = count_reads(read, counter)
    log2ratio = math.log2((counter['f1r2']+1)/(counter['f2r1']+1))
    if sample:
        # only targets with at least 50 reads
        if sum([counter['f1r2'], counter['f2r1']]) >= 50:
            log2ratio_all.append(log2ratio)
        sum_counts = (counter['f1r2']+sum_counts[0], counter['f2r1']+sum_counts[1])
        print(ch, st, en, ge, counter['f1r2'], counter['f2r1'], log2ratio, sep="\t", file=out_all)
    else:
        print(ch, st, en, ge, counter['f1r2'], counter['f2r1'], log2ratio, sep="\t")

if sample:
    # get stats
    mean_log2ratio_all = sum(log_dif)/len(log_dif)
    var_log2ratio_all = sum([ (mean_log2ratio_all - x)**2 for x in log_dif ])/len(log_dif)
    std_log2ratio_all = math.sqrt(var_log2ratio_all)
    total_log2_ratio = math.log2(sum_counts[0]/sum_counts[1])
    # header
    print("total_F1R1","total_F2R1","total_log2_ratio","mean_log2_ratio_F1R2_over_F2R1_50plus","std_log2_ratio_F1R2_over_F2R1_50plus", sep="\t", file=out_sum)
    # print data
    print(sum_counts[0], sum_counts[1], total_log2_ratio, mean_log2ratio_all, std_log2ratio_all, sep="\t", file=out_sum)
