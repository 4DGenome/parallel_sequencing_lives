#!/usr/bin/python



#==================================================================================================
# Created on: 2017-06-23
# Usage: ./io_metadata.py <database> <mode> ...
# Author: javier.quilez@crg.eu, @jaquol
# Goal: download and manage metadata
#==================================================================================================

# IMPORTANT:
# in the presence of fields with (some) missing values, the following error can occur:
# >sqlalchemy.exc.StatementError: (exceptions.ValueError) could not convert string to float:
# it seems python/sql thinks the expected value is a float when in reality is a string
# the solution is delete that field in the table of the database and re-run the script

# Import python modules and functions
import os, sys
import dataset
import pandas as pd 
import numpy as np
from collections import OrderedDict
import glob
import time
import datetime

# Get variables from script parameters
metadata = sys.argv[1]
mode = sys.argv[2]

# Load database
db = dataset.connect('sqlite:///%s' % metadata)


#==================================================================================================
# COMMANDS
#==================================================================================================

if mode == "download_input_metadata":
	
	# read table with downloaded input metadata
	spreadsheet = sys.argv[3]
	df = pd.read_table(spreadsheet)
	
	# create/load table and add input metadata from multiple samples
	input_metadata = db.get_table('input_metadata', primary_id = 'SAMPLE_ID', primary_type = 'String')
	for index, row in df.iterrows():
		input_metadata.upsert(row.to_dict(), ['SAMPLE_ID'])


elif mode == "quality_control_raw_reads":

	# get variables
	sample_id = sys.argv[3]
	attribute = sys.argv[4]
	value = sys.argv[5]

	# create/load table and add sample data
	# important: because this script may be used in parallel, with many instances trying to access to the data
	# do not use transactions (e.g. db.begin(), db.commit()...) as these do not allow multiple writings to the database
	quality_control_raw_reads = db.get_table('quality_control_raw_reads', primary_id = 'SAMPLE_ID', primary_type = 'String')
	new_data = {}
	new_data['SAMPLE_ID'] = sample_id
	new_data[attribute] = value
	quality_control_raw_reads.upsert(new_data, ['SAMPLE_ID'])


elif mode == "get_from_metadata":

	# get variables
	table = sys.argv[3]
	sample_id = sys.argv[4]
	attribute = sys.argv[5]
	#print metadata, mode, table, sample_id, attribute
	# load table, select sample and print attribute
	tab = db.load_table(table)
	my_sample = tab.find(SAMPLE_ID = sample_id)
	for s in my_sample:
		print s[attribute]


elif mode == "print_freeze":

	date = time.strftime("%Y-%m-%d")
	CONSEQ = sys.argv[3]
	ODIR = '%s/metadata/freezes/%s' % (CONSEQ, date)
	if not os.path.exists(ODIR):
		os.makedirs(ODIR)

	# print tables
	tables = db.tables
	for t in tables:
		print 'making freeze for table %s' % t
		otab = '%s/%s.csv' % (ODIR, t)
		result = db[t].all()
		dataset.freeze(result, format = 'csv', filename = otab)

elif mode == "add_to_metadata":

	# get variables
	table = sys.argv[3]
	sample_id = sys.argv[4]
	time_stamp = sys.argv[5]
	attribute = sys.argv[6]
	value = sys.argv[7]
	my_key = ";".join([sample_id, time_stamp])

	# important: because this script may be used in parallel, with many instances trying to access to the data
	# do not use transactions (e.g. db.begin(), db.commit()...) as these do not allow multiple writings to the database

	# hic table stores the most recent values
	tab = db.get_table(table, primary_id = 'SAMPLE_ID', primary_type = 'String')
 	new_data = {}
 	new_data['SAMPLE_ID'] = sample_id
 	new_data[attribute] = value
 	tab.upsert(new_data, ['SAMPLE_ID'])

 	# jobs table stores all values from different jobs
	tab = db.get_table('jobs', primary_id = 'JOB_ID', primary_type = 'String')
 	new_data = {}
 	new_data['JOB_ID'] = my_key
 	new_data[attribute] = value
 	tab.upsert(new_data, ['JOB_ID'])
