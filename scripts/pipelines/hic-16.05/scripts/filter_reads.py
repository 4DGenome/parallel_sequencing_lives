#!/usr/bin/python


# Created on: Aug 6th, 2015
# Usage: python filter_reads.py
# Goal: filter aligned reads



# Import python modules/functions
import glob, sys
import pandas as pd
from pytadbit.mapping.filter import filter_reads
from pytadbit.mapping.filter import apply_filter

# Get variables from script parameters
PROCESSED = sys.argv[1]
FILTERED = sys.argv[2]
DANGLING = sys.argv[3]
SELF_CIRCLE = sys.argv[4]
SUMMARY_EXCLUDED = sys.argv[5]
max_molecule_length = int(sys.argv[6])
over_represented = float(sys.argv[7])
min_frag_size = int(sys.argv[8])
max_frag_size = int(sys.argv[9])
re_proximity = int(sys.argv[10])
both_reads_mapped = int(sys.argv[11])

# Count number of reads excluded by each filter (note that a read can be included in more than one filter!)
infile = glob.glob('%s/*_both_map.tsv' % PROCESSED)[0]
pair_id = infile.split("/")[-1].replace("_both_map.tsv", "")
masked = filter_reads(infile, max_molecule_length=max_molecule_length,
                                    over_represented=over_represented,
                                    min_frag_size=min_frag_size,
                                    max_frag_size=max_frag_size,
                                    re_proximity=re_proximity,
                                    verbose=False)
filters_applied_numeric = [1,2,3,4,9,10]
is_applied = []
my_columns = ['filter_index', 'exclusion', 'reads_number', 'reads_fraction']
excluded = pd.DataFrame(columns=my_columns)
for k in xrange(1, len(masked) + 1):
    df = pd.DataFrame([k, masked[k]['name'], masked[k]['reads']]).transpose()
    df.columns = my_columns[:-1]
    df['reads_fraction'] = df['reads_number'] / both_reads_mapped
    excluded = pd.concat([excluded, df])
    if k in filters_applied_numeric:
    	is_applied.append(1)
    else:
    	is_applied.append(0)
    excluded['is_applied'] = is_applied
excluded['exclusion'] = [e.replace(' ', '_') for e in excluded['exclusion']]
outfile = '%s/%s_summary_excluded_per_filter.txt' % (SUMMARY_EXCLUDED, pair_id)
excluded.to_csv(outfile, sep='\t', index=False)

# Output reads
filtered_reads = '%s/%s_filtered_map.tsv' % (FILTERED, pair_id)
n_filtered_reads = apply_filter(infile, filtered_reads, masked, filters=filters_applied_numeric, verbose=False)
dangling_ends = '%s/%s_dangling_ends_map.tsv' % (DANGLING, pair_id)
apply_filter(infile, dangling_ends, masked, filters=[2], reverse=True, verbose=False)
self_circle = '%s/%s_self_circle_map.tsv' % (SELF_CIRCLE, pair_id)
apply_filter(infile, self_circle, masked, filters=[1], reverse=True, verbose=False)                                    

# Print the number of reads that are kept after fileting
print n_filtered_reads