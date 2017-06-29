#!/usr/bin/python


# Created on: Aug 5th, 2015
# Usage: python map_process_merge.py
# Goal: Map, Process mapped reads according to restriction enzyme fragments and Merge mapped "read1" and "read2"



# Import python modules/functions
import subprocess, sys, os
from pytadbit.mapping.full_mapper import full_mapping


gem_index = sys.argv[1]
SAMPLE = sys.argv[2]
species = sys.argv[3]
read_length = sys.argv[4]
paired1 = sys.argv[5]
paired2 = sys.argv[6]
restriction_enzyme = sys.argv[7]
fasta = sys.argv[8]
slots = sys.argv[9]
frag_map = sys.argv[10]
version = sys.argv[11]

print gem_index
print SAMPLE
print species
print read_length
print paired1
print paired2
print restriction_enzyme
print fasta
print slots
print frag_map
print version


# ========================================================================================
# Mapping
# ========================================================================================

# Generate GEM mapper index if it does not exist
if not os.path.isfile(gem_index):
#    fasta = '%s/db/reference_genome/%s/%s/%s.fa' % (PROJECT, species, assembly_version, assembly_version)
    subprocess.call(["gemtools", "index", "-i", fasta, "-t", "8"])

# Output directories
MAP_DIR = '%s/mapped_reads/%s' % (SAMPLE, version)
if not os.path.exists(MAP_DIR):
    os.makedirs(MAP_DIR)

# Make windows according to the mapping mode
if frag_map == 'True':
    frag_map = True
    windows = None
elif frag_map == 'False':
    frag_map = False
    range_stop = range(20, int(read_length)+1, 5)
    range_start = [1] * len(range_stop)
    windows = (zip(*(range_start, range_stop)))

# call mapping function for read1 and read2
for infile in [paired1, paired2]:
    bname = infile.split("/")[-1].replace(".fastq.gz", "")
    maps = full_mapping(
             gem_index_path      = gem_index,
             fastq_path          = infile,
             out_map_dir         = '%s/%s/' % (MAP_DIR, bname),
             r_enz               = restriction_enzyme,
             windows             = windows,
             temp_dir            = '%s/tmp_dir_%s/' % (MAP_DIR, bname),
             frag_map            = frag_map,
             nthreads            = slots)



# ========================================================================================
# Process mapped reads according to restriction enzyme fragments, Merging mapped "read1" and "read2"
# ========================================================================================

# Import python modules/functions
import glob
from pytadbit.parsers.map_parser    import parse_map
from pytadbit.parsers.genome_parser import parse_fasta
from pytadbit.mapping        import get_intersection

# Load the genome
genome_seq = parse_fasta(fasta)

# Output directory
RESULTS = '%s/results/%s/processed_reads' % (SAMPLE, version)
if not os.path.exists(RESULTS):
    os.makedirs(RESULTS)

infiles = []
outfiles = []
for infile in [paired1, paired2]:
    bname = infile.split("/")[-1].replace(".fastq.gz", "")
    maps = glob.glob('%s/%s/*' % (MAP_DIR, bname))
    infiles.append(maps)
    outfiles.append('%s/%s_map.tsv' % (RESULTS, bname))

parse_map(infiles[0], infiles[1], outfiles[0], outfiles[1], genome_seq, restriction_enzyme, verbose = True, ncpus=slots)
final_output = outfiles[0].replace('read1', 'both')
get_intersection(outfiles[0], outfiles[1], final_output, verbose = True)
