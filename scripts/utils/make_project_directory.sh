#!/bin/bash


#==================================================================================================
# Created on: 2017-06-25
# Usage: ./make_project_directory.sh
# Author: javier.quilez@crg.eu
# Goal: makes project directory
#==================================================================================================

PSL=/users/GR/mb/jquilez/projects/parallel_sequencing_lives

# variables
project=$1

# check variables are passed as script parameters
if [ -n "$project" ]; then
	PROJECT=$PSL/projects/$project
	if [ ! -d $PROJECT ]; then
		# make directories
		mkdir -p $PROJECT/{data,analysis}
		# make project notebook
		md=$PROJECT/project_notebook_${project}.md
		rm -f $md
		echo "# $project" >> $md
		echo `printf '%100s\n' | tr ' ' -` >> $md
		echo -e "\n**objective: ...**" >> $md
		echo -e "\n**paths are relative to $PSL**\n\n" >> $md
		echo -e "\nproject directory created at $PROJECT\n"
		echo "# Project directory stucture" >> $md
		echo `printf '%100s\n' | tr ' ' -` >> $md
		echo >> $md
		echo "- analysis: subdirectories for the different analyses, named as <date_of_analysis>_<analysis_name>" >> $md
		echo "- data: input data that are not sample-specific (which are are at /users/GR/mb/jquilez/)" >> $md
		echo "- project_notebook_<project_name>.md: this file..." >> $md
		echo >> $md
		echo >> $md
	else
		echo -e "\n$PROJECT already exists\n"
		exit
	fi
else
	echo -e "\nusage: make_analysis_directory.sh <project>\n"
	exit
fi
