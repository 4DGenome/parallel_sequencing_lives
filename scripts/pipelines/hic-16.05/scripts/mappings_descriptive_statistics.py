#!/usr/bin/python


# Created on: Aug 6th, 2015
# Usage: python mappings_descriptive_statistics.py
# Goal: generate plots summarizing quality metrics of the mapped reads



# Import python modules/functions
import sys
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from matplotlib import pyplot as plt 
from matplotlib import rcParams
import pandas as pd

import numpy as np
from pytadbit.mapping.analyze import plot_iterative_mapping
from pytadbit.mapping.analyze import insert_sizes
from pytadbit.mapping.analyze import plot_distance_vs_interactions
from pytadbit.mapping.analyze import hic_map
from collections import OrderedDict
from module_plotting_functions import plot_genomic_distribution
#from module_plotting_functions import plot_proportion_mapped_reads

# Get variables from script parameters
PROCESSED = sys.argv[1]
POSTMAPPING_PLOTS = sys.argv[2]
COVERAGES = sys.argv[3]
maps1 = sys.argv[4]
maps2 = sys.argv[5]
n_reads_trimmed = int(sys.argv[6])
genomic_coverage_resolution = sys.argv[7]

# Plotting parameters
plt.rcParams['font.size'] = 20 
plt.rcParams['font.weight'] = 'medium'
#plt.rcParams['font.family'] = 'sans-serif' 
#plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['lines.linewidth'] = 2.0
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['legend.frameon'] = False
plt.rcParams['savefig.bbox'] = 'tight'

# Plot: fraction of mapped reads
infiles = [maps1, maps2]
pair_id = infiles[0].split("/")[-1].replace("_read1_map.tsv", "")
outfile = '%s/%s_plot_proportion_mapped_reads.png' % (POSTMAPPING_PLOTS, pair_id)
reads_mapped_per_iteration = plot_iterative_mapping(fnam1=infiles[0], fnam2=infiles[1], total_reads=n_reads_trimmed, savefig=outfile)
reads_mapped_per_iteration = pd.DataFrame.from_dict(reads_mapped_per_iteration)
reads_mapped_per_iteration.columns = ['read1', 'read2']
fraction_mapped_read1 = list(reads_mapped_per_iteration['read1'])[-1] / float(n_reads_trimmed)
fraction_mapped_read2 = list(reads_mapped_per_iteration['read2'])[-1] / float(n_reads_trimmed)
fraction_mapped_str = ",".join([str(i) for i in [fraction_mapped_read1, fraction_mapped_read2]])

# Plot: distribution of dangling-end lengths
plt.rcParams['font.size'] = 12
infile = '%s/%s_both_map.tsv' % (PROCESSED, pair_id)
outfile = '%s/%s_plot_distribution_dangling_ends_lengths.png' % (POSTMAPPING_PLOTS, pair_id)
insert_sizes(infile, xlog=False, max_size=99.9, savefig=outfile)

# Plot: Decay of interaction counts with genomic distamce
plt.rcParams['font.size'] = 12 
outfile = '%s/%s_plot_decay_interaction_counts_genomic_distance.png' % (POSTMAPPING_PLOTS, pair_id)   
myvalues = plot_distance_vs_interactions(infile, max_diff=50000000, resolution=10000, savefig=outfile)
slope = str(myvalues[1][0])

# Plot: sequencing coverage along chromosomes
outfile = '%s/%s_plot_genomic_coverage_mapped_%s.png' % (POSTMAPPING_PLOTS, pair_id, genomic_coverage_resolution)
plt.rcParams['font.size'] = 20
coverages = plot_genomic_distribution(infile, name='mapped', savefig=outfile, resolution=genomic_coverage_resolution, pair_id=pair_id)
outfile = '%s/%s_plot_genomic_coverage_mapped_%s.bed' % (COVERAGES, pair_id, genomic_coverage_resolution)
coverages.to_csv(outfile, sep='\t', index=False)

# Plot: interaction matrix
bin_size = {}
bin_size['Mb'] = 1e6
outfile = '%s/%s_plot_interaction_matrix_mapped_%s.png' % (POSTMAPPING_PLOTS, pair_id, genomic_coverage_resolution)
plt.rcParams['font.size'] = 12
hic_map(infile, resolution=int(bin_size[genomic_coverage_resolution]), savefig=outfile, decay=False, cmap='jet')

print ";".join([fraction_mapped_str, slope])