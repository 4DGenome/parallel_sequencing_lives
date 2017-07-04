#!/bin/bash
#$ -N quality_control_b7fa2d8db_bfac48760
#$ -q short-sl7
#$ -l virtual_free=15G
#$ -l h_rt=06:00:00
#$ -o /users/GR/mb/jquilez/projects/parallel_sequencing_lives/projects/jquilez/analysis/2017-06-28_process_hic_samples/job_out/quality_control_b7fa2d8db_bfac48760_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/projects/parallel_sequencing_lives/projects/jquilez/analysis/2017-06-28_process_hic_samples/job_out/quality_control_b7fa2d8db_bfac48760_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 1
/software/mb/bin/fastqc --extract /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2016-02-11/b7fa2d8db_bfac48760_read1.fastq.gz -o /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2016-02-11/fastqc; rm -f /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2016-02-11/fastqc/b7fa2d8db_bfac48760*read1_fastqc.zip
/software/mb/bin/fastqc --extract /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2016-02-11/b7fa2d8db_bfac48760_read2.fastq.gz -o /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2016-02-11/fastqc; rm -f /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2016-02-11/fastqc/b7fa2d8db_bfac48760*read2_fastqc.zip
/users/GR/mb/jquilez/projects/parallel_sequencing_lives/scripts/utils/io_metadata.sh -m quality_control_raw_reads -s b7fa2d8db_bfac48760 -p PE
