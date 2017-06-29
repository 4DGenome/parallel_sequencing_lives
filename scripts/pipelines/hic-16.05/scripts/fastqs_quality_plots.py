#!/usr/bin/python


# Created on: Aug 6th, 2015
# Usage: python fastqs_quality_plots.py
# Goal: generate per-FASTQ plots showing quality metrics for the reads


# Import python modules/functions
import sys
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from matplotlib import pyplot as plt 
from pytadbit.utils.fastq_utils import quality_plot

QUALITY_PLOTS = sys.argv[1]
paired1 = sys.argv[2]
paired2 = sys.argv[3]
reads_nummber_qc = sys.argv[4]
restriction_enzyme = sys.argv[5]

# Generate quality plots for each processed FASTQ
# and print out the percentage of dangling-ends and ligated sites for each read
values = [] 
for infile in [paired1, paired2]:
    bname = infile.split("/")[-1].replace(".fastq.gz", "")
    outfile = '%s/%s_processed_reads_quality.png' % (QUALITY_PLOTS, bname)
    a, b = quality_plot(infile, nreads=int(reads_nummber_qc), r_enz=restriction_enzyme, savefig=outfile)    
    values.append(a)
    values.append(b)
print ';'.join([str(i) for i in values])