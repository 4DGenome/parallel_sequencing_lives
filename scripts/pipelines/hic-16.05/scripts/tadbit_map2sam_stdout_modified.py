#!/usr/bin/env python

#
# Gets the *_both_filled_map.tsv contacts from TADBIT (and the corresponding filter files)
#  and outputs a modified SAM with the following fields:
#
#    - read ID
#    - filtering flag (see codes in header)
#    - chromosome ID of the first pair of the contact
#    - genomic position of the first pair of the contact
#    - mapped length of the first pair of the contact
#    - pseudo CIGAR with strand orientation of the first pair (0M = -, 1M = +)
#    - chromosome ID of the second pair of the contact
#    - genomic position of the second pair of the contact
#    - SIGNED mapped length of the first pair of the contact (sign ~ orientation)
#    - sequence is missing (*)
#    - quality is missing (*)
#    - TC tag indicating single (1) or multi (2) contact
#
# Each pair of contacts porduces two lines in the output SAM
#


# dependencies

import sys
import os
import collections

# get arguments

input = sys.argv
infile = os.path.realpath(input[1])

# infile = "/users/project/4DGenome/data/MB_1_1/results/processed_reads/" + \
#          "MB_1_1_10082_CGATGT_both_filled_map.tsv"

# outfile = "/mnt/cluster_4DGenome/data/MB_1_1/results/processed_reads/" + \
#          "MB_1_1_10082_CGATGT_both_filled_map.sam"

# get full paths

basename = os.path.basename(infile)
dirname = os.path.dirname(infile)

# define filter codes

filter_keys = collections.OrderedDict()
filter_keys['self_circle'] = 2 ** 0
filter_keys['dangling_end'] = 2 ** 1
filter_keys['error'] = 2 ** 2
filter_keys['extra_dangling_end'] = 2 ** 3
filter_keys['too_close_from_RES'] = 2 ** 4
filter_keys['too_short'] = 2 ** 5
filter_keys['too_large'] = 2 ** 6
filter_keys['over_represented'] = 2 ** 7
filter_keys['duplicated'] = 2 ** 8
filter_keys['random_breaks'] = 2 ** 9
filter_keys['trans'] = 2 ** 10

# translate map + flag into hic-sam (two lines per contact)

def map2sam (line, flag):

    (qname,
     rname, pos, s1, l1, e1, e2,
     rnext, pnext, s2, l2, e3, e4) = line.strip().split('\t')

    # store mapped length and strand in a single value

    mapq = str(int(l1) * ((-1) ** (1 + int(s1))))
    tlen = str(int(l2) * ((-1) ** (1 + int(s2))))

    # multicontact?

    tc = str(int("~" in qname) + 1)
    
    # trans contact?

    if(rname != rnext):

        flag += filter_keys['trans'] 

    return ((qname, str(flag), rname, pos, str(abs(int(mapq))), s1 + "M",
             rnext, pnext, tlen, "*", "*",
             "TC:i:" + tc, "E1:i:" + e1, "E2:i:" + e2, "E3:i:" + e3, "E4:i:" + e4),
            (qname, str(flag), rnext, pnext, str(abs(int(tlen))), s2 + "M",
             rname, pos, mapq, "*", "*",
             "TC:i:" + tc, "E3:i:" + e3, "E4:i:" + e4, "E1:i:" + e1, "E2:i:" + e2))


# get all filters

filter_files = {}

for file in os.listdir(dirname):

    if file.startswith(basename + "_"):

        key = file.replace(basename + "_", "").replace(".tsv", "")
        filter_files[key] = dirname + "/" + file

# write header

print "\t".join(("@HD" ,"VN:1.5", "SO:queryname"))

fhandler = open(infile)
line = fhandler.next()

# chromosome lengths

while line.startswith('#'):

    (_, _, cr, ln) = line.replace("\t", " ").strip().split(" ")

    print "\t".join(("@SQ", "SN:" + cr, "LN:" + ln))

    line = fhandler.next()

# filter codes

for i in filter_keys:

    print "\t".join(("@CO", "filter:" + i, "flag:" + str(filter_keys[i])))

# tags

print "\t".join(("@CO" ,"TC:i", "Multicontact? 0 = no 1 = yes"))
print "\t".join(("@CO" ,"E1:i", "Position of the left RE site of first read"))
print "\t".join(("@CO" ,"E2:i", "Position of the right RE site of first read"))
print "\t".join(("@CO" ,"E3:i", "Position of the left RE site of second read"))
print "\t".join(("@CO" ,"E4:i", "Position of the right RE site of second read"))


# open and init filter files

filter_handler = {}
filter_line = {}

for i in filter_files:
    filter_handler[i] = open(filter_files[i])
    try:
        filter_line[i] = filter_handler[i].next().strip()
    except StopIteration:
        filter_line[i] = ""

n = 0

try:

    while True:

        flag = 0

        # check if read matches any filter

        for i in filter_line:

            if filter_line[i] == line.split("\t")[0]:

                flag += filter_keys[i.replace("-", "_")]

                try:

                    filter_line[i] = filter_handler[i].next().strip()

                except StopIteration:

                    pass

        # get output in sam format

        (out1, out2) = map2sam(line, flag)

        print "\t".join(out1)
        print "\t".join(out2)

        line = fhandler.next()
        n += 1

except StopIteration:

    pass

# close file handlers

fhandler.close()
for i in filter_files:
    filter_handler[i].close()
