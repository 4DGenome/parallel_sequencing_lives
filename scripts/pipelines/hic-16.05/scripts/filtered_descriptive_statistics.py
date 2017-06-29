#!/usr/bin/python


# Created on: Aug 6th, 2015
# Usage: python filtered_descriptive_statistics.py
# Goal: generate plots summarizing quality metrics of the filtered/excluded reads



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
from pytadbit import load_hic_data_from_reads
from collections import OrderedDict
from module_plotting_functions import plot_genomic_distribution
from module_plotting_functions import plot_proportion_mapped_reads

# Get variables from the script parameters
filtered_reads = sys.argv[1]
dangling_ends = sys.argv[2]
self_circle = sys.argv[3]
POSTFILTERING_PLOTS = sys.argv[4]
COVERAGES = sys.argv[5]
genomic_coverage_resolution = sys.argv[6]

# Plotting parameters
plt.rcParams['font.size'] = 20 
plt.rcParams['font.weight'] = 'medium' 
plt.rcParams['lines.linewidth'] = 2.0
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['legend.frameon'] = False
plt.rcParams['savefig.bbox'] = 'tight'

bin_size = {}
bin_size['Mb'] = 1e6

pair_id = filtered_reads.split("/")[-1].replace("_filtered_map.tsv", "")

outfile = '%s/%s_plot_genomic_coverage_filtered_%s.png' % (POSTFILTERING_PLOTS, pair_id, genomic_coverage_resolution)
plt.rcParams['font.size'] = 20
coverages = plot_genomic_distribution(filtered_reads, name='filtered', savefig=outfile, resolution=genomic_coverage_resolution, pair_id=pair_id)
outfile = '%s/%s_genomic_coverage_filtered_%s.bed' % (COVERAGES, pair_id, genomic_coverage_resolution)
coverages.to_csv(outfile, sep='\t', index=False)
# Filtered reads: interaction matrix
outfile = '%s/%s_plot_interaction_matrix_filtered_%s.png' % (POSTFILTERING_PLOTS, pair_id, genomic_coverage_resolution)
plt.rcParams['font.size'] = 12
hic_map(filtered_reads, resolution=int(bin_size[genomic_coverage_resolution]), savefig=outfile, decay=False, cmap='jet')  

# Dangling ends: sequencing coverage along chromosomes
outfile = '%s/%s_plot_genomic_coverage_dangling_ends_%s.png' % (POSTFILTERING_PLOTS, pair_id, genomic_coverage_resolution)
plt.rcParams['font.size'] = 20
coverages = plot_genomic_distribution(dangling_ends, name='dangling_ends', savefig=outfile, resolution=genomic_coverage_resolution, pair_id=pair_id)
outfile = '%s/%s_genomic_coverage_dangling_ends_%s.bed' % (COVERAGES, pair_id, genomic_coverage_resolution)
coverages.to_csv(outfile, sep='\t', index=False)

# self-circle ends: sequencing coverage along chromosomes
outfile = '%s/%s_plot_genomic_coverage_self_circle_%s.png' % (POSTFILTERING_PLOTS, pair_id, genomic_coverage_resolution)
plt.rcParams['font.size'] = 20
coverages = plot_genomic_distribution(self_circle, name='self_circle', savefig=outfile, resolution=genomic_coverage_resolution, pair_id=pair_id)
outfile = '%s/%s_genomic_coverage_self_circle_%s.bed' % (COVERAGES, pair_id, genomic_coverage_resolution)
coverages.to_csv(outfile, sep='\t', index=False)
