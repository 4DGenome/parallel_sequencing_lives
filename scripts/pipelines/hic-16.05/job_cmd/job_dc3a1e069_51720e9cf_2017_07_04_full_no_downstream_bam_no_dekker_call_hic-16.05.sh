#!/bin/bash
#$ -N job_dc3a1e069_51720e9cf_2017_07_04_full_no_downstream_bam_no_dekker_call_hic-16.05
#$ -q long-sl7
#$ -l virtual_free=100G
#$ -l h_rt=100:00:00
#$ -o /users/GR/mb/jquilez/projects/parallel_sequencing_lives/scripts/pipelines/hic-16.05/job_out/job_dc3a1e069_51720e9cf_2017_07_04_full_no_downstream_bam_no_dekker_call_hic-16.05_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/projects/parallel_sequencing_lives/scripts/pipelines/hic-16.05/job_out/job_dc3a1e069_51720e9cf_2017_07_04_full_no_downstream_bam_no_dekker_call_hic-16.05_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 10

submitted_on=2017_07_04
pipeline_version=16.05
sample_id=dc3a1e069_51720e9cf
data_type=hic
pipeline_name=hic
pipeline_version=16.05
pipeline_run_mode=full_no_downstream_bam_no_dekker_call
io_mode=standard
CUSTOM_IN=scripts/pipelines/hic-16.05/test
CUSTOM_OUT=scripts/pipelines/hic-16.05/test
sample_to_fastqs=sample_to_fastqs.txt
submit_to_cluster=no
queue=long-sl7
memory=100G
max_time=100:00:00
slots=10
email=javier.quilez@crg.eu
integrate_metadata=yes
species=
version=
read_length=
sequencing_type=PE
seedMismatches=2
palindromeClipThreshold=30
simpleClipThreshold=12
leading=3
trailing=3
minAdapterLength=1
keepBothReads=true
minQual=3
strictness=0.999
minLength=36
restriction_enzyme=
max_molecule_length=500
max_frag_size=10000
min_frag_size=50
over_represented=0.005
re_proximity=4
reads_number_qc=100000
genomic_coverage_resolution=Mb
frag_map=True
flag_excluded=783
flag_included=0
flag_perzero=99
resolution_tad=50000
resolution_ab=100000
pis=500000
pids=250000
pnt=0.1
CUSTOM_OUT=scripts/pipelines/hic-16.05/test
PIPELINE=/users/GR/mb/jquilez/projects/parallel_sequencing_lives/scripts/pipelines/hic-16.05
config=scripts/pipelines/hic-16.05/hic.config
path_job_file=/users/GR/mb/jquilez/projects/parallel_sequencing_lives/scripts/pipelines/hic-16.05/job_cmd/job_dc3a1e069_51720e9cf_2017_07_04_full_no_downstream_bam_no_dekker_call_hic-16.05.sh

PSL=/users/GR/mb/jquilez/projects/parallel_sequencing_lives

# =================================================================================================
# CODE EXECUTION
# =================================================================================================

main() {

	echo

	#==================================================================================================
	# Configuration
	#==================================================================================================

	# additional run variables
	time_start=$(date +"%s")
	run_date=`date +"%Y-%m-%d-%H-%M"`
	job_name=$pipeline_name-$pipeline_version

	# makes that the job uses this python/tadbit
	source activate latest_tadbit
	#export PATH="/software/mb/el7.2/anaconda2/bin:$PATH"

	# pipeline scripts
	SCRIPTS=$PIPELINE/scripts

	# sequencing data
	SEQ_DATA=$PSL/data/hic/raw

	# Settings based on the io_mode parameter
	if [[ $io_mode == "custom" ]]; then
		message_info "configuration" "running $job_name with io_mode=custom, which implies:"
		message_info "configuration" "integrate_metadata=no (value given in the *.config file is overwritten)"
		message_info "configuration" "all parameter values are taken from the *.config file and used for all samples"
		integrate_metadata=no
		# output directory
		SAMPLE=$CUSTOM_OUT/$sample_id
		message_info "configuration" "data for $sample_id will be stored at $SAMPLE"
		# input FASTQ files
		ifq1_name=`grep -w $sample_id $CUSTOM_IN/sample_to_fastqs.txt |cut -f2`
	 	ifq2_name=`grep -w $sample_id $CUSTOM_IN/sample_to_fastqs.txt |cut -f3`
	 	ifq1=$CUSTOM_IN/$ifq1_name
	 	ifq2=$CUSTOM_IN/$ifq2_name
	 	message_info "configuration" "paths to input FASTQ files extracted from $CUSTOM_IN/sample_to_fastqs.txt"
		message_info "configuration" "species ($species) and assembly version ($version) extracted from the metadata"
		if [[ $species == 'homo_sapiens' ]]; then
			fasta=/users/GR/mb/jquilez/assemblies/$species/$version/ucsc/${version}_chr1-22XYM.fa
		elif [[ $species == 'mus_musculus' ]]; then
			fasta=/users/GR/mb/jquilez/assemblies/$species/$version/ucsc/${version}_chr1-19XYM.fa
		elif [[ $species == 'drosophila_melanogaster' ]]; then
			fasta=/users/GR/mb/jquilez/assemblies/$species/$version/ucsc/${version}_chr2-4XYM.fa
		elif [[ $species == 'caulobacter_crescentus' ]]; then
			fasta=/users/GR/mb/jquilez/assemblies/$species/$version/ncbi/${version}.fa			
		fi
	elif [[ $io_mode == "standard" ]]; then
		# script to in/out data from metadata
		io_metadata=$PSL/scripts/utils/io_metadata.sh
		metadata=$PSL/metadata/metadata.db
		message_info "configuration" "script to access/write from/to metadata: $io_metadata"
		message_info "configuration" "metadata database: $metadata"
		# output directory
		SAMPLE=$PSL/data/$data_type/samples/$sample_id
		message_info "configuration" "data for $sample_id will be stored at $SAMPLE"
		# input FASTQ files
		ifq1=$SEQ_DATA/*/${sample_id}*read1.fastq.gz
		ifq2=$SEQ_DATA/*/${sample_id}*read2.fastq.gz
		message_info "configuration" "paths to input FASTQ files extracted from the metadata"
		# get species and assembly version from the metadata
		if [[ $integrate_metadata == "yes" ]]; then
			species=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a 'SPECIES'`
			restriction_enzyme=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a RESTRICTION_ENZYME`
			read_length=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a READ_LENGTH`
			if [[ ${species,,} == 'homo_sapiens' ]]; then
				version=hg38_mmtv
				fasta=/users/GR/mb/jquilez/assemblies/${species,,}/$version/ucsc/${version}_chr1-22XYM.fa
			elif [[ ${species,,} == 'mus_musculus' ]]; then
				version=mm10
				fasta=/users/GR/mb/jquilez/assemblies/${species,,}/$version/ucsc/${version}_chr1-19XYM.fa
			elif [[ ${species,,} == 'drosophila_melanogaster' ]]; then
				#version=dm6
				#fasta=/users/GR/mb/jquilez/assemblies/${species,,}/$version/ucsc/${version}_chr2-4XYM.fa
				version=dm3
				fasta=/users/GR/mb/jquilez/assemblies/${species,,}/$version/flybase/${version}.fa
			fi
		fi
		message_info "configuration" "species ($species) and assembly version ($version) extracted from the metadata"
	fi

	# Logs
	LOGS=$SAMPLE/logs/$version

	# Trim reads
	TRIMMED=$SAMPLE/fastqs_processed/trimmomatic
	PAIRED=$TRIMMED/paired_end
	UNPAIRED=$TRIMMED/unpaired_reads
	ADAPTERS=/software/mb/el7.2/Trimmomatic-0.36/adapters

	# SHA cheksums
	CHECKSUMS=$SAMPLE/checksums/$version/$run_date
	checksums=$CHECKSUMS/files_checksums.sha

	# quality plots of raw reads
	QUALITY_PLOTS=$SAMPLE/plots/$version/raw_fastqs_quality_plots

	# aligned, processed and merged reads and post-mapping quality plots
	PROCESSED=$SAMPLE/results/$version/processed_reads
	POSTMAPPING_PLOTS=$SAMPLE/plots/$version/post_mapping_statistics
	COVERAGES=$SAMPLE/results/$version/genomic_coverages

	# filtered and excluded reads, post-filtering plots
	FILTERED=$SAMPLE/results/$version/filtered_reads
	DANGLING=$SAMPLE/results/$version/excluded_reads/dangling_ends
	SELF_CIRCLE=$SAMPLE/results/$version/excluded_reads/self_circle
	SUMMARY_EXCLUDED=$SAMPLE/results/$version/excluded_reads/summary_excluded_per_filter
	POSTFILTERING_PLOTS=$SAMPLE/plots/$version/post_filtering_statistics

	# downstream analyses
	DOWNSTREAM=$SAMPLE/downstream/$version

	# tools
	shasum=`which shasum`
	python=`which python`
	trimmomatic=`which trimmomatic`
	samtools=`which samtools`
	bgzip=`which bgzip`
	tabix=`which tabix`
	java=`which java`

	# get TADbit and its dependencies versions
	tadbit_and_dependencies_versions=`$python $SCRIPTS/print_tadbit_and_dependencies_version.py`

	# store general parameters into the metadata
	if [[ $integrate_metadata == "yes" ]]; then
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a PIPELINE_RUN_MODE -v $pipeline_run_mode
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a QUEUE -v $queue
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a MEMORY -v $memory
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a MAX_TIME -v $max_time
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a SLOTS -v $slots
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a ASSEMBLY_VERSION -v $version
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a JOB_NAME -v $job_name		
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a PATH_JOB_FILE -v $path_job_file		
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a TADBIT_AND_DEPENDENCIES_VERSIONS -v $tadbit_and_dependencies_versions
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a LAST_RUN_DATE -v $run_date		
	fi

	echo

	#==================================================================================================
	# Execute modules
	#==================================================================================================

	if [[ $pipeline_run_mode == 'full' ]]; then
		preliminary_checks
		raw_fastqs_quality_plots
		trim_reads_trimmomatic
		align_and_merge
		post_mapping_statistics
		reads_filtering
		post_filtering_statistics
		map_to_bam
		downstream_bam
		dekker_call
		clean_up
	elif [[ $pipeline_run_mode == 'preliminary_checks' ]]; then preliminary_checks
	elif [[ $pipeline_run_mode == 'raw_fastqs_quality_plots' ]]; then raw_fastqs_quality_plots
	elif [[ $pipeline_run_mode == 'trim_reads_trimmomatic' ]]; then trim_reads_trimmomatic
	elif [[ $pipeline_run_mode == 'align_and_merge' ]]; then align_and_merge
	elif [[ $pipeline_run_mode == 'post_mapping_statistics' ]]; then post_mapping_statistics
	elif [[ $pipeline_run_mode == 'reads_filtering' ]]; then reads_filtering
	elif [[ $pipeline_run_mode == 'post_filtering_statistics' ]]; then post_filtering_statistics
	elif [[ $pipeline_run_mode == 'map_to_bam' ]]; then map_to_bam
	elif [[ $pipeline_run_mode == 'downstream_bam' ]]; then downstream_bam
	elif [[ $pipeline_run_mode == 'dekker_call' ]]; then dekker_call	
	elif [[ $pipeline_run_mode == 'clean_up' ]]; then clean_up
	elif [[ $pipeline_run_mode == 'full_no_clean_up' ]]; then
		preliminary_checks
		raw_fastqs_quality_plots
		trim_reads_trimmomatic
		align_and_merge
		post_mapping_statistics
		reads_filtering
		post_filtering_statistics
		map_to_bam
		downstream_bam
		dekker_call
	elif [[ $pipeline_run_mode == 'full_no_downstream_bam_no_dekker_call' ]]; then
		preliminary_checks
		raw_fastqs_quality_plots
		trim_reads_trimmomatic
		align_and_merge
		post_mapping_statistics
		reads_filtering
		post_filtering_statistics
		map_to_bam
		clean_up
	fi
	echo

	# Final message
	message_info "pipeline" "completed successfully"
	message_time_pipeline 

}



#==================================================================================================
# FUNCTIONS DEFINITIONS
#==================================================================================================


# =================================================================================================
# Pipeline progress functions
# =================================================================================================

# Outputs a message about the task being done
message_info() {
	step_name=$1
	message=$2
	echo -e "INFO \t`date +"%Y-%m-%d %T"` \t[$step_name] \t$message"
}

# Outputs a message about the error found and exits
message_error() {
	step_name=$1
	message=$2
	echo -e "ERROR \t`date +"%Y-%m-%d %T"` \t[$step_name] \t$message"
	exit	
}

# Outputs a warning message about the task being done
message_warn() {
	step_name=$1
	message=$2
	echo -e "WARN \t`date +"%Y-%m-%d %T"` \t[$step_name] \t$message"
}

# Outputs a message with the time in seconds spent in a given step
message_time_step() {
	step_name=$1
	field_name="TIME_${step_name^^}"
	time0=$2
	time1=$(date +"%s")
	length=$(($time1-$time0))
	echo -e "TIME \t`date +"%Y-%m-%d %T"` \t[$step_name] \tstep time for completion (seconds) = $length"
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a $field_name -v $length
	fi
	echo
}

# Outputs the total time in seconds for the pipeline to run
message_time_pipeline() {
	field_name="TIME_COMPLETE_PIPELINE"
	time0=$time_start
	time1=$(date +"%s")
	length=$(($time1-$time0))
	echo -e "TIME \t`date +"%Y-%m-%d %T"` \t[pipeline] \ttotal time for completion (seconds) = $length"
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a $field_name -v $length
	fi
	echo
}


# ========================================================================================
# Checks that the selected sample(s) is in the database
# Gets analysis parameters (for Trimmomatic and TADbit) from the configuration file
# Gets FASTQ files and additional metadata
# ========================================================================================

preliminary_checks() {

	step="preliminary_checks"
	time0=$(date +"%s")

	# Check that a sample with the provided SAMPLE_ID exists
	if [[ $integrate_metadata == "yes" ]]; then
		sample_check=`$io_metadata -m check_sample -s $sample_id`
		if [[ $sample_check == "no_sample" ]]; then
			message_error $step "$sample_id not found in $metadata. Exiting..."
		elif [[ $sample_check == "multiple_samples" ]]; then
			message_error $step "$sample_id has multiple entries in $metadata. Exiting..."
		else
			message_info $step "$sample_id found in $metadata"
		fi
	fi

	# check FASTQ files exist
	if [ -f $ifq1 ] && [ -f $ifq2 ]; then
		mkdir -p $SAMPLE
		mkdir -p $CHECKSUMS
		# save a SHA checksums of the FASTQ files
		# because different compressions would then have different checksums, checksums are generated on the uncompressed FASTQ
		# after checksums are generated, the uncompressed file is deleted for the sake of space
		#tfq1=`echo $ifq1 | sed "s/\.gz//g"`
		#tfq2=`echo $ifq2 | sed "s/\.gz//g"`
		zcat $ifq1 |$shasum | awk -v file=`echo $ifq1 | sed 's/.gz//g'` '{print $1,"",file}' >> $checksums
		zcat $ifq2 |$shasum | awk -v file=`echo $ifq2 | sed 's/.gz//g'` '{print $1,"",file}' >> $checksums
		#zcat $ifq1 > $tfq1
		#zcat $ifq2 > $tfq2
		#shasum $ifq1 >> $checksums
		#shasum $tfq1 >> $checksums
		#shasum $ifq2 >> $checksums
		#shasum $tfq2 >> $checksums
		#rm -f $tfq1 $tfq2
		# Get sequencing information from the header of the FASTQ reads (some fields should be shared across all reads)
		fq_header=`zcat $ifq1 | head -n 1 | sed s'/ /:/g'`
		sequencing_instrument_name=`echo $fq_header | cut -f1 -d':'`
		sequencing_run_id=`echo $fq_header | cut -f2 -d':'`
		sequencing_flowcell_id=`echo $fq_header | cut -f3 -d':'`
		sequencing_flowcell_lane=`echo $fq_header | cut -f4 -d':'`
		sequencing_index_fq=`echo $fq_header | cut -f11 -d':'`
		# update metadata
		if [[ $integrate_metadata == "yes" ]]; then
			$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a PATH_FASTQ_READ1 -v $ifq1
			$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a PATH_FASTQ_READ2 -v $ifq2
			$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a 'SEQUENCING_INSTRUMENT_NAME' -v $sequencing_instrument_name
			$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a 'SEQUENCING_RUN_ID' -v $sequencing_run_id
			$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a 'SEQUENCING_FLOWCELL_ID' -v $sequencing_flowcell_id
			$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a 'SEQUENCING_FLOWCELL_LANE' -v $sequencing_flowcell_lane
			message_info $step "paths to read1 and read2 saved to metadata database"
			# Check that the sequencing index introduced as part of the metadata agrees with that found in the FASTQ reads
			sequencing_index=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a SEQUENCING_INDEX`
			if [[ $sequencing_index == $sequencing_index_fq ]]; then
		 		message_info $step "Sequencing index added as part of the metadata agrees with that observed in the FASTQ file"
			else
				message_warn $step "Sequencing index added as part of the metadata does not agree with that observed in the FASTQ file"
			fi
			# Determine read length and update metadata if necessary
			# Read length from the first read of the FASTQ files (it assumes all reads in the file have the same lenght!)
			first_read=`zcat $ifq1 | head -2 | tail -n 1`	
			read_length_fq=${#first_read}
			# Read length from the metadata
			read_length_metadata=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a 'READ_LENGTH'`
			if [[ "$read_length_fq" != "$read_length_metadata" ]]; then
				message_warn $step "read length obtained from the metadata ($read_length_metadata bp) differs from that seen in the FASTQ ($read_length_fq bp)"
				message_warn $step "the latter will be kept"
				$io_metadata -m update_input_metadata -s $sample_id -a 'READ_LENGTH' -v $read_length_fq
			else
				message_info $step "read length added as part of the metadata agrees with that observed in the FASTQ file"
			fi
		fi
		# Check that FASTQ for the reference genome sequence exists
		if ! [[ -f "$fasta" ]]; then
			message_error $step "FASTA file $fasta does not exist! Exiting..."
		else
			message_info $step "genome reference FASTA file found at $fasta"
			$shasum $fasta >> $checksums
			if [[ $integrate_metadata == "yes" ]]; then
				$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a PATH_REFERENCE_FASTA -v $fasta
				message_info $step "paths to reference FASTA saved to metadata database"
			fi
		fi
	else
		message_error $step "$ifq1 and/or $ifq2 not found. Exiting..."
	fi

	message_time_step $step $time0

}


# ========================================================================================
# Quality plots of the raw reads
# ========================================================================================

raw_fastqs_quality_plots() {

	step="raw_fastqs_quality_plots"
	time0=$(date +"%s")

	message_info $step "restriction enzyme = $restriction_enzyme"
	message_info $step "a subset of $reads_number_qc reds will be used to generate the quality plots"
	
	# Make plots
	mkdir -p $QUALITY_PLOTS
	message_info $step "making quality plots of the raw FASTQ files... plots saved at $QUALITY_PLOTS"
	returned_values=`$python $SCRIPTS/fastqs_quality_plots.py $QUALITY_PLOTS \
  											$ifq1 \
  											$ifq2 \
  											$reads_number_qc \
  											$restriction_enzyme`

 	# Extract the percentage of dangling-ends and ligated sites for read1 and read2
  	percentage_dangling_ends_read1=`echo $returned_values |cut -f1 -d';'`
  	percentage_ligated_sites_read1=`echo $returned_values |cut -f2 -d';'`
  	percentage_dangling_ends_read2=`echo $returned_values |cut -f3 -d';'`
  	percentage_ligated_sites_read2=`echo $returned_values |cut -f4 -d';'`

  	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a PERCENTAGE_DANGLING_ENDS_READ1 -v $percentage_dangling_ends_read1
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a PERCENTAGE_LIGATED_SITES_READ1 -v $percentage_ligated_sites_read1
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a PERCENTAGE_DANGLING_ENDS_READ2 -v $percentage_dangling_ends_read2
		$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a PERCENTAGE_LIGATED_SITES_READ2 -v $percentage_ligated_sites_read2
		message_info $step "percentages of dangling ends and ligated sites in read1 and read2 saved to metadata database"
	fi
 
	message_time_step $step $time0

}


# =================================================================================================
# Trim adapter and low-quality ends
# =================================================================================================

trim_reads_trimmomatic() {

	step="trim_reads_trimmomatic"
	time0=$(date +"%s")

	# paths
	mkdir -p $PAIRED
	mkdir -p $UNPAIRED
	mkdir -p $LOGS
	mkdir -p $CHECKSUMS
	step_log=$SAMPLE/logs/${sample_id}_${step}_paired_end.log
	paired1=$PAIRED/${sample_id}_read1.fastq.gz
	paired2=$PAIRED/${sample_id}_read2.fastq.gz
	unpaired1=$UNPAIRED/${sample_id}_read1.fastq.gz
	unpaired2=$UNPAIRED/${sample_id}_read2.fastq.gz
	params="$ifq1 $ifq2 $paired1 $unpaired1 $paired2 $unpaired2"
	ODIR=$PAIRED

	# Get metadata
	#if [[ $integrate_metadata == "yes" ]]; then
	#	read_length=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a READ_LENGTH`
	#fi

	# adapter trimming: the trimmomatic program directory contains a folder with the adapter sequences for
	# the Illumina sequencers in use. 'TruSeq3-PE.fa' is used, which contains the adapter sequences for the HiSeq
	message_info $step "sequencing type = $sequencing_type" 
	message_info $step "trimming TruSeq3 adapter sequences for HiSeq, NextSeq or HiSeq"
	message_info $step "trimming low-quality reads ends using trimmomatic's recommended practices"
	seqs=$ADAPTERS/TruSeq3-$sequencing_type.fa
	targetLength=$read_length
	$java -jar $trimmomatic $sequencing_type \
 					$params \
 					ILLUMINACLIP:$seqs:$seedMismatches:$palindromeClipThreshold:$simpleClipThreshold:$minAdapterLength:$keepBothReads \
 					LEADING:$leading \
 					TRAILING:$trailing \
 					MAXINFO:$targetLength:$strictness \
 					MINLEN:$minLength >$step_log 2>&1

	# parse step log to extract generated metadata
	message_info $step "parse step log to extract generated metadata"
 	n_reads_trimmed=`grep Surviving $step_log | cut -f3 -d':' | cut -f1 -d'(' | sed "s/ //g"`
	message_info $step "reads after trimming = $n_reads_trimmed"

	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a ADAPTERS_SEQS -v $seqs
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a SEED_MISMATCHES -v $seedMismatches
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a PALINDROME_CLIP_THRESHOLD -v $palindromeClipThreshold
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a SIMPLE_CLIP_THRESHOLD -v $simpleClipThreshold
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a MIN_ADAPTER_LENGTH -v $minAdapterLength
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a KEEP_BOTH_READS -v $keepBothReads
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a LEADING -v $leading
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a TRAILING -v $trailing
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a TARGET_LENGTH -v $targetLength
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a STRICTNESS -v $strictness
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a MIN_LENGTH -v $minLength
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a N_READS_TRIMMED -v $n_reads_trimmed
		message_info $step "trimmomatic parameters and numbe of trimmed reads added to metadata"
	fi

	# delete intermediate files
	message_info $step "trimmed reads are in $ODIR"
	if [[ $sequencing_type == "PE" ]]; then
		message_info $step "unpaired reads are deleted"
		rm -fr $UNPAIRED
	fi

	message_time_step $step $time0

}


# ========================================================================================
# align with GEM, process mapped reads according to restriction enzyme fragments and Merge mapped "read1" and "read2"
# ========================================================================================

align_and_merge() {

	step="align_and_merge"
	time0=$(date +"%s")

	# paths
	mkdir -p $CHECKSUMS

	# Get metadata
	#if [[ $integrate_metadata == "yes" ]]; then
	#	restriction_enzyme=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a RESTRICTION_ENZYME`
	#	species=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a SPECIES`
	#	read_length=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a READ_LENGTH`
	#fi
	message_info $step "species = $species"
	message_info $step "assembly version = $version"
	message_info $step "restriction enzyme = $restriction_enzyme"
	message_info $step "read length = $read_length"

	# genome reference FASTA
	#fasta=/users/GR/mb/jquilez/assemblies/${species,,}/$version/ucsc/$version.fa
	if ! [[ -f "$fasta" ]]; then
		message_error $step "FASTA file $fasta does not exist! Exiting..."
	else
		message_info $step "genome reference FASTA file found at $fasta"
	fi

	# processed FASTQs (i.e. after `trim_reads_trimmomatic`)
	paired1=$PAIRED/$sample_id*read1.fastq.gz
	paired2=$PAIRED/$sample_id*read2.fastq.gz

	# Added to overcome en error when running TADbit 'get_intersection' function
	export LC_ALL=en_US.UTF-8
	export LANG=en_US.UTF-8

	step_log=$LOGS/${sample_id}_${step}_paired_end.log
	# # Prevent overriding existing processed reads
	# if [ -f $step_log ]; then
    #		message_error $step "processed reads already mapped (see $step_log). Exiting..."
	# fi

	# Mapping
	message_info $step "mapping, processing reads according to restriction enzyme fragments and merging aligments for read1 and read2..."
	gem_index=`echo $fasta | sed "s/\.fa/\.gem/g"`
	$python $SCRIPTS/map_process_merge.py $gem_index $SAMPLE $species $read_length $paired1 $paired2 $restriction_enzyme $fasta $slots $frag_map $version > $step_log
	rm -fr $SAMPLE/mapped_reads/tmp_dir*
	rm -fr $SAMPLE/results/processed_reads/tmp*
	message_info $step "output saved in $step_log"

	# data integrity
	$shasum $PROCESSED/${sample_id}*both_map.tsv >> $checksums

	message_time_step $step $time0

}


# ========================================================================================
# Post-mapping descriptive statistics
# ========================================================================================

post_mapping_statistics() {

	step="post_mapping_statistics"
	time0=$(date +"%s")

	# paths
	maps1=$PROCESSED/$sample_id*read1*map.tsv
	maps2=$PROCESSED/$sample_id*read2*map.tsv
	mkdir -p $POSTMAPPING_PLOTS
	mkdir -p $COVERAGES

	# get metadata
	step_log_trim_reads_trimmomatic=$SAMPLE/logs/${sample_id}_trim_reads_trimmomatic_paired_end.log
 	n_reads_trimmed=`grep Surviving $step_log_trim_reads_trimmomatic | cut -f3 -d':' | cut -f1 -d'(' | sed "s/ //g"`

	# Generate plots
	message_info $step "generating post mapping descriptive statistics plots..."
	returned_values=`$python $SCRIPTS/mappings_descriptive_statistics.py $PROCESSED \
														$POSTMAPPING_PLOTS \
 														$COVERAGES \
 														$maps1 \
 														$maps2 \
 														$n_reads_trimmed \
 														$genomic_coverage_resolution`

	# parse metadata
	fraction_mapped_read1=`echo $returned_values | cut -f1 -d';' | cut -f1 -d','`
	fraction_mapped_read2=`echo $returned_values | cut -f1 -d';' | cut -f2 -d','`
	counts_to_distance_slope=`echo $returned_values | cut -f2 -d';'`
	both_reads_mapped=`cat $PROCESSED/$sample_id*both_map.tsv | grep -v "# CRM" | wc -l`
	message_info $step "fraction of reads mapped read1 = $fraction_mapped_read1"
	message_info $step "fraction of reads mapped read2 = $fraction_mapped_read2"
	message_info $step "number of pairs in which both read1 and read2 are mapped = $both_reads_mapped"
	message_info $step "counts-to-distance slope = $counts_to_distance_slope"

	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a FRACTION_MAPPED_READ1 -v $fraction_mapped_read1
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a FRACTION_MAPPED_READ2 -v $fraction_mapped_read2
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a COUNTS_TO_DISTANCE_SLOPE -v $counts_to_distance_slope
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a BOTH_READS_MAPPED -v $both_reads_mapped
	fi

	message_time_step $step $time0

}


# ========================================================================================
# Filtering of reads
# ========================================================================================

reads_filtering() {

	step="reads_filtering"
	time0=$(date +"%s")

	# Paths
	both_reads_mapped=`cat $PROCESSED/$sample_id*both_map.tsv | grep -v "#" | wc -l`
	mkdir -p $FILTERED
	mkdir -p $DANGLING 
	mkdir -p $SELF_CIRCLE
	mkdir -p $SUMMARY_EXCLUDED
	mkdir -p $CHECKSUMS

	# Filter reads based on multiple quality parameters
	message_info $step "filtering reads based on multiple quality parameters..."
	message_info $step "filtered reads to be used for subsequent analyses are at $FILTERED"
	message_info $step "excluded reads due to dangling ends are at $DANGLING"
	message_info $step "excluded reads due to self circle are at $SELF_CIRCLE"
	message_info $step "summary of excluded reads at $SUMMARY_EXCLUDED"
	n_filtered_reads=`$python $SCRIPTS/filter_reads.py $PROCESSED \
									$FILTERED \
									$DANGLING \
									$SELF_CIRCLE \
									$SUMMARY_EXCLUDED \
									$max_molecule_length \
									$over_represented \
									$min_frag_size \
									$max_frag_size \
									$re_proximity \
									$both_reads_mapped`

	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a N_FILTERED_READS -v $n_filtered_reads
		infile=$SUMMARY_EXCLUDED/*_summary_excluded_per_filter.txt
		for filter in `cut -f1 -d$'\t' $infile | grep -v exclusion`; do
			reads_number=`cat $infile | grep -w $filter | cut -f5`
			reads_fraction=`cat $infile | grep -w $filter | cut -f4`
			filter=`echo ${filter^^} | sed 's/-/_/g'`
			filter="EXCLUDED_$filter"
		 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a $filter -v "$reads_number;$reads_fraction"
		done
	 	message_info $step "numbers and fractions of filtered and excluded reads added to the metadata"	
	fi

	# data integrity
	filtered_reads=$FILTERED/*filtered_map.tsv
	$shasum $filtered_reads >> $checksums

	message_time_step $step $time0

}


# ========================================================================================
# Post-filtering descriptive statistics
# ========================================================================================

post_filtering_statistics() {

	step="post_filtering_statistics"
	time0=$(date +"%s")

	# Paths
	filtered_reads=$FILTERED/*filtered_map.tsv
	dangling_ends=$DANGLING/*dangling_ends_map.tsv
	self_circle=$SELF_CIRCLE/*self_circle_map.tsv
 	mkdir -p $POSTFILTERING_PLOTS

	# Make plots
	message_info $step "generating post filtering descriptive statistics plots... saved at $POSTFILTERING_PLOTS:"
	message_info $step "- filtered reads: sequencing coverage along chromosomes, coverage values and interaction matrix"
	message_info $step "- dangling ends: sequencing coverage along chromosomes and coverage values"
	message_info $step "- self-circle ends: sequencing coverage along chromosomes and coverage values"
	$python $SCRIPTS/filtered_descriptive_statistics.py $filtered_reads \
													$dangling_ends \
													$self_circle \
													$POSTFILTERING_PLOTS \
													$COVERAGES \
													$genomic_coverage_resolution

	message_time_step $step $time0

}


# ========================================================================================
# Convert MAP to BAM
# ========================================================================================

map_to_bam() {

	step="map_to_bam"
	time0=$(date +"%s")

	mkdir -p $CHECKSUMS
	message_info $step "converting MAP to BAM"
	imap=$PROCESSED/*both_map.tsv

	if [[ $io_mode == "custom" ]]; then
		ODIR=$PROCESSED
		obam=$ODIR/$(basename $imap .tsv).bam
		$python $SCRIPTS/tadbit_map2sam_stdout_modified.py $imap | $samtools view -Su - | $samtools sort - -o $obam
		$samtools index $obam
	elif [[ $io_mode == "standard" ]]; then
		ODIR=/users/project/4DGenome/data/hic/samples/$sample_id/results/$version/processed_reads
		mkdir -p $ODIR
		obam=$ODIR/$(basename $imap .tsv).bam
		$python $SCRIPTS/tadbit_map2sam_stdout_modified.py $imap | $samtools view -Su - | $samtools sort - -o $obam
		$samtools index $obam
		ln -sf $obam $PROCESSED/$(basename $obam)
		ln -sf $obam.bai $PROCESSED/$(basename $obam).bai
	fi

	# save a SHA checksums of the alignment file BAM
	# because different compressions would then have different checksums, checksums are generated on the uncompressed BAM (i.e. SAM)
	# after checksums are generated, the uncompressed file is deleted for the sake of space
	$samtools view -h $obam |$shasum | awk -v file=$obam.sam '{print $1,"",file}' >> $checksums
	#$samtools view -h $obam.bam > $obam.sam
	#$shasum $obam.bam >> $checksums
	#$shasum $obam.sam >> $checksums
	#rm -f $obam.map

	message_time_step $step $time0

}


# ========================================================================================
# Downstream BAM
# ========================================================================================

downstream_bam() {

	step="downstream_bam"
	time0=$(date +"%s")

	# paths
	ibam=$PROCESSED/*both_map.bam
	mkdir -p $DOWNSTREAM
	step_log=$LOGS/${sample_id}_${step}_paired_end.log

	# perform several downstream analyses
	message_info $step "perform several downstream analyses"
	$python $SCRIPTS/tadbit_after_bam_v2.py $ibam $flag_excluded $flag_included $flag_perzero $DOWNSTREAM/${sample_id}_ $slots $resolution_ab $resolution_tad &> $step_log

	# arrange, merge chromosomes, compress and index TADs
	obed=$DOWNSTREAM/$(basename $(ls -1 $DOWNSTREAM/${sample_id}_tads*tsv | head -n 1 |  sed 's,_chr.*$,,1'))_allchr.bed.gz
	grep -v "^#" $DOWNSTREAM/${sample_id}_tads_*.tsv | awk -v res=$resolution_tad -v OFS="\t" '{gsub(/^.*tads_/, "", $1); gsub(/_.*$/, "", $1); $2 = ($2 - 1) * res; $3 = ($3 - 1) * res; print}' | sort -k1,1 -k2,2n | $bgzip -fc > $obed
	$tabix -fp bed $obed
	rm $DOWNSTREAM/${sample_id}_tads_chr*
	n_tads=`zcat $obed | wc -l`
	message_info $step "number of TADs = $n_tads"

	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a FLAG_EXCLUDED -v $flag_excluded
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a FLAG_INCLUDED -v $flag_included
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a FLAG_PERZERO -v $flag_perzero
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a RESOLUTION_TAD -v $resolution_tad
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a RESOLUTION_AB -v $resolution_ab
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a N_TADS -v $n_tads
	fi
	
	message_time_step $step $time0

}


# ========================================================================================
# Dekker call
# ========================================================================================

dekker_call() {

	step="dekker_call"
	time0=$(date +"%s")

	# paths
	ibam=$PROCESSED/*both_map.bam
	mkdir -p $DOWNSTREAM
	step_log=$LOGS/${sample_id}_${step}_paired_end.log

	# call TADs with the Dekker method
	message_info $step "call TADs using Dekker's method"
	resolution_nice=`$python $SCRIPTS/nice.py $resolution_tad`
	MY_TMP=$DOWNSTREAM/my_tmp 
	mkdir -p $MY_TMP
	$SCRIPTS/dekker_call.r $DOWNSTREAM/${sample_id}_normalized_${resolution_nice}.tsv.gz $ibam $resolution_tad $slots $DOWNSTREAM/${sample_id} $pis $pids $pnt $MY_TMP &>>$step_log
	rm -fr $MY_TMP

	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a PIS -v $pis
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a PIDS -v $pids
	 	$io_metadata -m add_to_metadata -t 'hic' -s $sample_id -u $run_date -a PNT -v $pnt
	fi
	
	message_time_step $step $time0

}


# ========================================================================================
# Deletes intermediate files
# ========================================================================================

clean_up() {

	step="clean_up"
	time0=$(date +"%s")

	message_info $step "deleting the following intermediate files/directories:"
	message_info $step "$SAMPLE/fastqs_processed"
	message_info $step "$SAMPLE/mapped_reads"
	rm -fR $SAMPLE/fastqs_processed
	rm -fR $SAMPLE/mapped_reads
	rm -fR $SAMPLE/results/*/processed_reads/${sample_id}_read*_map.tsv
	message_time_step $step $time0

}


# execute main function
main
