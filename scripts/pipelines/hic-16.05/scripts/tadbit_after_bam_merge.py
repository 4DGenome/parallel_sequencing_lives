#!/usr/bin/env python

#
# Gets the contacts in BAM format, the filtering flags adn the desired resolutions (A/B and TADs)
#  and produce normalized and correlation matrices plus A/B compartments and TADs
#

##  Dependencies

from pytadbit                     import load_hic_data_from_reads
from pytadbit.utils.file_handling import mkdir
from pytadbit                     import tadbit
from pytadbit.parsers.tad_parser  import parse_tads
from pytadbit.utils.normalize_hic import iterative, expected
from pytadbit                     import HiC_data
from os                           import path, system
import subprocess, sys, os, re, numpy, math, pysam, datetime, time
from collections import OrderedDict
from itertools import izip

# reopen stdout file descriptor with write mode 
# and 0 as the buffer size (unbuffered) 

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) 

## Main

def main():

    # get arguments

    print printime() + "Main settings ..."

    filter_exclude = int(sys.argv[1])
    filter_include = int(sys.argv[2])
    perc_zero = int(sys.argv[3])
    outdir = sys.argv[4]
    ncpus = int(sys.argv[5])
    resolution_list = [int(sys.argv[6]), int(sys.argv[7])]
    bam_list = [i for i in sys.argv[8:]]


    # create dir

    # mkdir(outdir)

    # info message

    print "Take reads from"
    for i in bam_list:
        print(i)
    print "include contacts with flag = " + str(filter_include)
    print "exclude contacts with flag = " + str(filter_exclude)
    print "save results at " + str(outdir)
    print "(using " + str(ncpus) + " cpus and perc_zero = " + str(perc_zero) + ")"
    print "process A/B compartments at resolution " + str(resolution_list[0])
    print "process TADs at resolution " + str(resolution_list[1])

    # load data

    print printime() + "Loading data ..."

    dat_list, bin_list = bam_to_hic_data(bam_list, resolution_list, filter_exclude, filter_include)

    print printime() + "... done!"

    # process A/B compartments

    print printime() + "Processing A/B compartments ..."

    process_AB(dat_list[0], perc_zero, resolution_list[0], outdir, bin_list[0])

    print printime() + "... done!"

    # process TADs

    print printime() + "Processing TADs compartments ..."

    process_TAD(dat_list[1], perc_zero, resolution_list[1], ncpus, outdir, bin_list[1])

    print printime() + "... done!"

    print printime() + "FINISHED!"

# Define functions

def compress(a, outfile):
    of = open(outfile + ".tmp", "w")
    of.write('\n'.join(a))
    of.close()
    subprocess.check_call(["sort", "-k", "1,1", "-k", "2,2n", "-o", outfile, outfile + ".tmp"])
    subprocess.check_call(["rm", outfile + ".tmp"])
    subprocess.check_call(["/software/mb/bin/bgzip", "-f", outfile])
    subprocess.check_call(["/software/mb/bin/tabix", "-f", "-s", "1", "-b", "2", "-e", "2",
                           outfile + ".gz"])

def printime():
    return str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S\t'))

def nice(reso):
    if reso >= 1000000:
        return '%dMb' % (reso / 1000000)
    return '%dkb' % (reso / 1000)

def load_tad_height(tad_def, size, beg, end, hic_data):
    """
    to calculate tad densities
    """
    bias, zeros = hic_data.bias, hic_data.bads
    tads, _ = parse_tads(tad_def)
    diags = []
    for k in xrange(1, size):
        try:
            diags.append(sum(
                hic_data[i, i + k] / bias[i + k] / bias[i]
                for i in xrange(beg, end - k) if not i in zeros and not i + k in zeros
                ) / float(sum(1 for i in range(beg, end - k)
                              if not i in zeros and not i + k in zeros)))
        except ZeroDivisionError:
            diags.append(0.)
    for tad in tads:
        start, final = (int(tads[tad]['start']) + 1,
                        int(tads[tad]['end']) + 1)
        matrix = sum(
            hic_data[i, j] / bias[j] / bias[i]
            for i in xrange(beg + start - 1, beg + final - 1) if not i in zeros
            for j in xrange(i + 1          , beg + final - 1) if not j in zeros)
        try:
            height = float(matrix) / sum(
                [diags[i - 1] * (final - start - i)
                 for i in xrange(1, final - start)])
        except ZeroDivisionError:
            height = 0.
        tads[tad]['height'] = height
    return tads

def init_hic_data(dat, resolution):
    n = sum(dat.lengths) / resolution
    sections = OrderedDict(zip(dat.references,
                               [x / resolution + 1 for x in dat.lengths]))
    bins = []
    for crm in sections:
        len_crm = sections[crm]
        bins.extend([(crm, i) for i in xrange(len_crm)])
    bins_dict = dict([(j, i) for i, j in enumerate(bins)])
    dat = HiC_data((),
                   sum(sections.values()),
                   sections,
                   bins_dict,
                   resolution = resolution)
    return dat

def get_bins(dat, resolution):
    n = sum(dat.lengths) / resolution
    sections = OrderedDict(zip(dat.references,
                               [x / resolution + 1 for x in dat.lengths]))
    bins = []
    for crm in sections:
        len_crm = sections[crm]
        bins.extend([(crm, i) for i in xrange(len_crm)])
    bins_dict = dict([(j, i) for i, j in enumerate(bins)])
    return bins_dict

def bam_to_hic_data(bam_list, resolution_list, filter_exclude, filter_include):
    # get first bam
    inbam = bam_list[0]
    # open bam file
    bamfile = pysam.AlignmentFile(inbam, "rb")
    # get sections
    sections = OrderedDict(zip(bamfile.references,
                               bamfile.lengths))
    # init HiC_data objects
    dat_list = [init_hic_data(bamfile, resolution) for resolution in resolution_list]
    bin_list = [get_bins(bamfile, resolution) for resolution in resolution_list]
    # close bam file    
    bamfile.close()
    # access bam file per chromosome
    for chrom in sections.keys():
        print "Processing chromosome " + str(chrom)
        # per file
        for inbam in bam_list:
            j = 0
            print "Processing file " + str(inbam)
            # read line by line
            for line in pysam.view("-F", str(filter_exclude),
                                   "-f", str(filter_include),
                                   inbam,
                                   chrom):
                # get info
                info = line.strip().split("\t")
                chrom1 = info[2]
                pos1 = int(info[3])
                chrom2 = info[6]
                pos2 = int(info[7])
                if chrom2 == "=":
                    chrom2 = chrom1
                # get bins and add counts
                try:
                    for i in xrange(len(resolution_list)):
                        b1 = bin_list[i][(chrom1, int(pos1) / resolution_list[i])]
                        b2 = bin_list[i][(chrom2, int(pos2) / resolution_list[i])]
                        dat_list[i][b1, b2] += 1
                except KeyError:
                    pass
                j += 1
            print str(j) + " lines processed"
    return dat_list, bin_list

# functions to write raw/normalized/correlation matrices as compressed and indexed bed file

def write_matrix_tabix(hic_object, norm, outfile, reso):
    out = open(outfile, 'w')
    for sec in hic_object.chromosomes:
        lline = (hic_object.section_pos[sec][1] - hic_object.section_pos[sec][0]) * reso
        out.write(''.join([''.join(
            ['{}\t{}\t{}\t{}\n'.format(sec, i, j, val)
             for j, val in izip(xrange(i, lline, reso), line[i/reso:]) if val])
                           for i, line in izip(xrange(0, lline, reso),
                                               hic_object.yield_matrix(
                                                   focus=sec, normalized=norm))]))
    out.close()
    subprocess.check_call(["/software/mb/bin/bgzip", "-f", outfile])
    subprocess.check_call(["/software/mb/bin/tabix", "-f", "-p", "bed", outfile + ".gz"])


def write_correlation_matrix_tabix(hic_object, outfile, reso):

    out = open(outfile, 'w')
    for chromosome in hic_object.chromosomes.keys():
        if "M" in chromosome:
            continue
        mat = [[v / hic_object.expected[abs(i - j)] for j, v in enumerate(line)]
               for i, line in enumerate(hic_object.yield_matrix(normalized=True,
                                                                focus=chromosome))]
        if len(mat) == 1:
            continue
        mat = numpy.corrcoef(mat)

        lline = (hic_object.section_pos[chromosome][1] -
                 hic_object.section_pos[chromosome][0]) * reso
        out.write(''.join([''.join(
            ['{}\t{}\t{}\t{}\n'.format(chromosome, i, j, v)
             for j, v in izip(xrange(i, lline, reso), line[i/reso:])])
                           for i, line in izip(xrange(0, lline, reso), mat)]))
    out.close()

    subprocess.check_call(["/software/mb/bin/bgzip", "-f", outfile])
    subprocess.check_call(["/software/mb/bin/tabix", "-f", "-p", "bed", outfile + ".gz"])


def write_matrices(hic_data, outdir, reso):

    # Store matrices
    start = 0
    end = len(hic_data)
    print "Store raw data\n"
    write_matrix_tabix(hic_data, False, outdir + 'raw_%s.tsv' % nice(reso), reso)
    print "Store normalized data\n"
    write_matrix_tabix(hic_data, True,
                       outdir + 'normalized_%s.tsv' % nice(reso), reso)
    print "Compute correlation matrix and store it\n"
    write_correlation_matrix_tabix(hic_data, 
                                   outdir + 'correlation_%s.tsv' % nice(reso), reso)

def process_AB(hic_data, perc_zero, reso, outdir, bins):

    # Get poor bins

    print 'Get poor bins...'
    try:
        hic_data.filter_columns(perc_zero=perc_zero, by_mean=True)
    except ValueError:
        perc_zero = 100
        hic_data.filter_columns(perc_zero=perc_zero, by_mean=True)

    binsrev = {y:x for x,y in bins.iteritems()}

    bad_file = outdir + 'bad_rows_%s_%d.tsv' % (nice(reso), perc_zero)
    bads = [binsrev[i][0] + "\t" + str(binsrev[i][1] * reso) + "\t" + str(i)
            for i in hic_data.bads.keys()]

    compress(bads, bad_file)

    # Identify biases

    print 'Get biases using ICE...'

    hic_data.normalize_hic(silent=False, max_dev=0.1, iterations=0,
                           factor=1) # cells of the matrix have a mean of 1

    bias_file = outdir + 'bias_%s.tsv' % nice(reso)
    bias = [binsrev[i][0] + "\t" + str(binsrev[i][1] * reso) + "\t" + '%d\t%f' % (i, hic_data.bias[i])
            for i in hic_data.bias]

    compress(bias, bias_file)

    # percentage of cis interactions

    print 'Getting percentage of cis interactions...'

    cis_trans_N_D = hic_data.cis_trans_ratio(normalized=True , diagonal=True )
    cis_trans_n_D = hic_data.cis_trans_ratio(normalized=False, diagonal=True )
    cis_trans_N_d = hic_data.cis_trans_ratio(normalized=True , diagonal=False)
    cis_trans_n_d = hic_data.cis_trans_ratio(normalized=False, diagonal=False)

    cistrans_file = outdir + 'cis_trans_ratio_%s.tsv' % nice(reso)
    out_cistrans = open(cistrans_file, "w")

    out_cistrans.write("Cis/trans_ratio\tnormalized\twith_diagonal\t" + str(cis_trans_N_D) + "\n")
    out_cistrans.write("Cis/trans_ratio\tnormalized\twithout_diagonal\t" + str(cis_trans_N_d) + "\n")
    out_cistrans.write("Cis/trans_ratio\traw\twith_diagonal\t" + str(cis_trans_n_D) + "\n")
    out_cistrans.write("Cis/trans_ratio\traw\twithout_diagonal\t" + str(cis_trans_n_d) + "\n")
    out_cistrans.close()

    # Compute expected

    print 'Get expected counts ...'

    hic_data.expected = expected(hic_data, bads = hic_data.bads)

    # store matrices

    print 'Store matrices'
    write_matrices(hic_data, outdir, reso)

    # getting compartments

    print 'Searching compartments'

    ev = hic_data.find_compartments()
    ev_file = outdir + 'ev_%s.tsv' % nice(reso)

    out = []
    chroms = ev.keys()
    chroms.sort()
    for ch in chroms:
        for i in xrange(len(ev[ch][0]) - 1):
            out.append("\t".join((ch, str(i * reso), str(ev[ch][0][i]), str(ev[ch][1][i]))))

    compress(out, ev_file)

    cmprt_file = outdir + 'compartments_%s.tsv' % nice(reso)
    hic_data.write_compartments(cmprt_file)


def process_TAD(hic_data, perc_zero, reso, cpus, outdir, bins):

    # Get poor bins

    print 'Get poor bins...'

    try:

        hic_data.filter_columns(perc_zero=perc_zero, by_mean=True)

    except ValueError:

        perc_zero = 100
        hic_data.filter_columns(perc_zero=perc_zero, by_mean=True)

    binsrev = {y:x for x,y in bins.iteritems()}

    bad_file = outdir + 'bad_rows_%s_%d.tsv' % (nice(reso), perc_zero)
    bads = [binsrev[i][0] + "\t" + str(binsrev[i][1] * reso) + "\t" + str(i) for i in hic_data.bads.keys()]

    compress(bads, bad_file)

    # Identify biases

    print 'Get biases using ICE...'

    hic_data.normalize_hic(silent=False, max_dev=0.1, iterations=0,
                           factor=1) # cells of the matrix have a mean of 1

    bias_file = outdir + 'bias_%s.tsv' % nice(reso)
    bias = [binsrev[i][0] + "\t" + str(binsrev[i][1] * reso) + "\t" + '%d\t%f' % (i, hic_data.bias[i]) for i in hic_data.bias]

    compress(bias, bias_file)

    # percentage of cis interactions

    print 'Getting percentage of cis interactions...'

    cis_trans_N_D = hic_data.cis_trans_ratio(normalized=True , diagonal=True )
    cis_trans_n_D = hic_data.cis_trans_ratio(normalized=False, diagonal=True )
    cis_trans_N_d = hic_data.cis_trans_ratio(normalized=True , diagonal=False)
    cis_trans_n_d = hic_data.cis_trans_ratio(normalized=False, diagonal=False)

    cistrans_file = outdir + 'cis_trans_ratio_%s.tsv' % nice(reso)

    out_cistrans = open(cistrans_file, "w")
    out_cistrans.write("Cis/trans_ratio\tnormalized\twith_diagonal\t" + str(cis_trans_N_D) + "\n")
    out_cistrans.write("Cis/trans_ratio\tnormalized\twithout_diagonal\t" + str(cis_trans_N_d) + "\n")
    out_cistrans.write("Cis/trans_ratio\traw\twith_diagonal\t" + str(cis_trans_n_D) + "\n")
    out_cistrans.write("Cis/trans_ratio\traw\twithout_diagonal\t" + str(cis_trans_n_d) + "\n")
    out_cistrans.close()

    # Compute expected

    print 'Get expected counts ...'

    hic_data.expected = expected(hic_data, bads = hic_data.bads)

    # store matrices

    print 'Store matrices'

    write_matrices(hic_data, outdir, reso)

    # getting TAD borders

    print 'Searching TADs'

    for crm in hic_data.chromosomes:

        print '  - %s' % crm

        matrix = hic_data.get_matrix(focus=crm)
        beg, end = hic_data.section_pos[crm]
        size = len(matrix)

        if size < 10:
            print "     Chromosome too short (%d bins), skipping..." % size
            continue

        # transform bad column in chromosome referential

        remove = tuple([1 if i in hic_data.bads else 0 for i in xrange(beg, end)])

        # maximum size of a TAD

        max_tad_size = size

        result = tadbit([matrix], remove=remove,
                        n_cpus=cpus, verbose=False,
                        max_tad_size=max_tad_size,
                        no_heuristic=0)
        
        tads = load_tad_height(result, size, beg, end, hic_data)

        table = ''
        table += '%s\t%s\t%s\t%s%s\n' % ('#', 'start', 'end', 'score', 'density')

        for tad in tads:

            table += '%s\t%s\t%s\t%s%s\n' % (
                tad, int(tads[tad]['start'] + 1), int(tads[tad]['end'] + 1),
                abs(tads[tad]['score']), '\t%s' % (round(
                    float(tads[tad]['height']), 3)))

        out_tad = outdir + 'tads_%s_%s.tsv' % (
            crm, nice(reso))

        out = open(out_tad, 'w')
        out.write(table)
        out.close()

##  Run!

if __name__ == "__main__":
    exit(main())

