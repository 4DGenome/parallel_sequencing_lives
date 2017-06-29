#!/usr/bin/python


#==================================================================================================
# Created on: 2017-06-28
# Usage: ./sample_id_generator.py
# Author: Javier Quilez (GitHub: jaquol)
# Goal: generates unique sample identifier (SAMPLE_ID) based on the samples metadata
#==================================================================================================



#==================================================================================================
# CONFIGURATION VARIABLES AND IMPORT DATA
#==================================================================================================

# Import python modules and functions
import pandas as pd
import sys
import hashlib

# Get variables from script parameters
spreadsheet = sys.argv[1]
mode = sys.argv[2]

# read metadata
df = pd.read_table(spreadsheet)

# required fields
biological_fields = ['CELL_TYPE', 'PRE_TREATMENT', 'PRE_TREATMENT_TIME', 'TREATMENT', 'TREATMENT_TIME', 'CONTROL']
experiment_fields = ['EXPERIMENT_ID']
technical_fields = ['HIC', 'SEQUENCING_PRE_LIBRARY', 'SEQUENCING_LIBRARY', 'SEQUENCING_CORE', 'SEND_FOR_SEQUENCING_ON']



#==================================================================================================
# FUNCTIONS
#==================================================================================================

# checks if there are samples with missing values in the fields used to generate unique SAMPLE_ID values
def warn_missing_values(df):
	df_nonmissing = df[biological_fields + experiment_fields + technical_fields].dropna()
	n1 = len(df)
 	n2 = len(df_nonmissing)
 	if n1 != n2:
 		print '[warn]\twarning: %i samples with missing values in required fields used to generate SAMPLE_ID' % (n1 - n2)
 		print '[warn]\tunique SAMPLE_ID values may be inaccurate!'


# concatenates the values (as strings and semicolon-separated) from the selected list of fields
def make_string_from_fields(row, list_fields):
	list_fields_converted = []
	for l in list_fields:
		value = row[l]
		if type(value) is str:
			list_fields_converted.append(value)
		elif type(value) == int:
			list_fields_converted.append(str(value))
		elif type(value) == float:
			list_fields_converted.append(str(int(value)))
	return ";".join(list_fields_converted)


# converts any-size string to unique string of length 11 characters
# function tested/borrowed from:
# http://www.peterbe.com/plog/best-hashing-function-in-python
# I use h11 instead of h6 (which produces shorter strings) as the second may produce non-ASCII characters e.g. '=', which we want to avoid
#def h6(w):
#    h = hashlib.md5(w)
#    return h.digest().encode('base64')[:6]
def h11(w):
    return hashlib.md5(w).hexdigest()[:9]


# apply to samples which do not have yet a SAMPLE_ID
def assign_sample_id(df):

	samples = df.copy()

	# check if there are missing values
	warn_missing_values(samples)

    # generate unique SAMPLE_ID, made of 2 parts:
    # part1 = string combining values (semicolon-separated) from the 5 biological fields
    # part2 = string combining values (semicolon-separated) from the 1 experiment field plus the 5 biological fields
    # the string for each part is *digested* to make a unique string of 11 ASCII characters
    # the two parts are separate by underscore	
	for index, row in samples.iterrows():
		time_stamp = row['Timestamp']
		part1 = h11(make_string_from_fields(row, biological_fields))
		part2 = h11(make_string_from_fields(row, experiment_fields + technical_fields))
		sample_id = '%s_%s' % (part1, part2)
		print '[info]\tSAMPLE_ID for sample with Timestamp = %s is: %s' % (time_stamp, sample_id)

# generate SAMPLE_ID values according to the selected mode
assign_sample_id(df)