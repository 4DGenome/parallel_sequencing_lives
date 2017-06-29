'''
Python modules for the analysis of Hi-C data (created on: July 7th, 2015)

'''


import sys
import pandas as pd
import numpy as np
import multiprocessing as mu
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from matplotlib import pyplot as plt 
from collections import OrderedDict



def plot_proportion_mapped_reads(reads1=None, reads2=None, mappable_reads=None, pair_id=None, savefig=None):

    # Plot 'scaffold'
    plt.close('all')
    fig, ax1 = plt.subplots(figsize=(10, 10))
    colors = ['olive', 'darkcyan']
    labels = ['Read1', 'Read2']
    # Extract the number of mapped reads in each iteration
    fraction_mapped = []
    for i, f in enumerate([reads1, reads2]):
        iterations = []
        mapped = []
        with open(f) as f:
            for line in f:
                if line.startswith('#'):
                    if line.startswith('# MAPPED'):
                        it, m = line.split()[2:]
                        iterations.append(int(it))
                        mapped.append(int(m))
                else:
                    break
        # Plot data points
        ax1.plot(iterations, np.cumsum(mapped) / (1. * mappable_reads),
            color=colors[int(i)-1],
            label=labels[int(i)-1])
        fraction_mapped.append(round(np.cumsum(mapped)[-1] / (1. * mappable_reads), ndigits=2))
    ax1.set_xlabel('Iteration number', labelpad = 25)
    ax1.set_ylabel('Proportion of mapped reads', labelpad = 25)
    ax1.set_ylim(0, 1)
    ax1.legend(loc='lower right')
    ax1.text(0.5, 1.10, 'FASTQ ID = %s' % pair_id, transform=ax1.transAxes, horizontalalignment='center')
    ax1.text(0.5, 1.05, 'No. mappable reads per read = %i' % mappable_reads, transform=ax1.transAxes, horizontalalignment='center')
    plt.savefig(savefig, bbox_inches = 'tight', dpi = 300)
    return fraction_mapped



def plot_genomic_distribution(fnam, name=None, first_read=True, resolution=None,
                              axe=None, ylim=None, savefig=None,
                              chr_names=None, nreads=None, pair_id=None):
    """
    :param fnam: input file name
    :param True first_read: uses first read.
    :param 100 resolution: group reads that are closer than this resolution
       parameter
    :param None axe: a matplotlib.axes.Axes object to define the plot
       appearance
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    :param None chr_names: can pass a list of chromosome names in case only some
       them the need to be plotted (this option may last even more than default)
    
    """

    bin_size = {}
    bin_size['Mb'] = 1e6

    name_to_color = {}
    name_to_color['mapped'] = 'gray'
    name_to_color['filtered'] = 'blue'
    name_to_color['dangling_ends'] = 'red'
    name_to_color['self_circle'] = 'red'

    # Parse input file
    distr = {}
    idx1, idx2 = (1, 3) if first_read else (7, 9)
    genome_seq = OrderedDict()
    fhandler = open(fnam)
    line = fhandler.next()
    if chr_names:
        chr_names = set(chr_names)
        cond1 = lambda x: x not in chr_names
    else:
        cond1 = lambda x: False
    if nreads:
        cond2 = lambda x: x >= nreads
    else:
        cond2 = lambda x: False
    cond = lambda x, y: cond1(x) and cond2(y)
    count = 0
    while line.startswith('#'):
        if line.startswith('# CRM '):
            crm, clen = line[6:].split('\t')
            genome_seq[crm] = int(clen)
        line = fhandler.next()
    try:
        while True:
            crm, pos = line.strip().split('\t')[idx1:idx2]
            count += 1
            if cond(crm, count):
                line = fhandler.next()
                if cond2(count):
                    break
                continue
            pos = int(pos) / int(bin_size[resolution])
            try:
                distr[crm][pos] += 1
            except KeyError:
                try:
                    distr[crm][pos] = 1
                except KeyError:
                    distr[crm] = {pos: 1}
            line = fhandler.next()
    except StopIteration:
        pass
    fhandler.close()

    # Data frame that will store read count values
    my_columns = ['chrom', 'start', 'end', 'read_counts']
    df = pd.DataFrame(columns=my_columns)

    # Plot scaffold
    plt.close('all')
    fig, ax = plt.subplots(nrows=len(genome_seq), figsize=(15,len(genome_seq)), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.1)  
    max_ylim = max([np.percentile(distr[crm].values(), 90) for crm in distr.keys()])
    max_xlim = max([max(genome_seq.values()) for crm in distr.keys()]) / int(bin_size[resolution]) * 1.1
    fig.text(0.05, 0.5, 'Number of reads\n(%s-bins)' % resolution, rotation=90, horizontalalignment='center')
    fig.text(0.5, 0.9125, 'FASTQ ID = %s (%s)' % (pair_id, name), horizontalalignment='center')
    for i, crm in enumerate(chr_names if chr_names else genome_seq
                            if genome_seq else distr):

        #if crm == "chrM":
        #    continue

        # Add genomic values
        if crm in distr:
            ax[i].plot(range(max(distr[crm])), [distr[crm].get(j, 0) for j in xrange(max(distr[crm]))],
                        color=name_to_color[name], lw=2, alpha=0.70)
        # Axes and labels        
        ax[i].set_xlim(0, max_xlim)
        ax[i].set_ylim(0, max_ylim)
        try:
            tmp = max(range(max(distr[crm])))
        except ValueError:
            tmp = 1 # esto engloba los casos de coordenada maxima = 0
        except KeyError:
            tmp = 0 # esto en caso de que no este el crm dentro de distr
        #if crm in distr:
        #    tmp = max(range(max(distr[crm]))) if len(distr[crm]) > 1 else 1
        #else:
        #    tmp = 0
        ax[i].text(tmp+15, max_ylim/2, crm, verticalalignment='center',
                horizontalalignment='center', fontsize = 20, weight='bold')
        if i == (len(genome_seq)-1):
            ax[i].set_xlabel('Position on chromosome (%s)' % resolution, labelpad=25)
        for tick in ax[i].yaxis.get_major_ticks():
            tick.label.set_fontsize(10)
        # Add read counts to data frame
        if crm in distr:
            df_crm = pd.DataFrame(distr[crm].items(), columns=['size_bin', 'read_counts'])
            df_crm['chrom'] = crm
            df_crm['start'] = df_crm['size_bin'] * bin_size[resolution]
            df_crm['end'] = df_crm['start'] + bin_size[resolution] - 1
            df = pd.concat([df, df_crm[my_columns]])
    plt.locator_params(nbins=4)
    if savefig:
       plt.savefig(savefig, bbox_inches = 'tight', dpi = 300)
    
    df[['start', 'end', 'read_counts']] = df[['start', 'end', 'read_counts']].astype(int)    
    return df