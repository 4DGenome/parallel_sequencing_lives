#!/bin/bash
#$ -N quality_control_dc3a1e069_51720e9cf
#$ -q short-sl7
#$ -l virtual_free=15G
#$ -l h_rt=06:00:00
#$ -o /users/GR/mb/jquilez/projects/didactic_dataset/projects/jquilez/analysis/2017-06-28_process_hic_samples/job_out/quality_control_dc3a1e069_51720e9cf_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/projects/didactic_dataset/projects/jquilez/analysis/2017-06-28_process_hic_samples/job_out/quality_control_dc3a1e069_51720e9cf_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 1
/software/mb/bin/fastqc --extract /users/GR/mb/jquilez/projects/didactic_dataset/data/hic/raw/8/1/2015/dc3a1e069_51720e9cf*read1.fastq.gz -o /users/GR/mb/jquilez/projects/didactic_dataset/data/hic/raw/8/1/2015/fastqc; rm -f /users/GR/mb/jquilez/projects/didactic_dataset/data/hic/raw/8/1/2015/fastqc/dc3a1e069_51720e9cf*read1_fastqc.zip
/software/mb/bin/fastqc --extract /users/GR/mb/jquilez/projects/didactic_dataset/data/hic/raw/8/1/2015/dc3a1e069_51720e9cf*read2.fastq.gz -o /users/GR/mb/jquilez/projects/didactic_dataset/data/hic/raw/8/1/2015/fastqc; rm -f /users/GR/mb/jquilez/projects/didactic_dataset/data/hic/raw/8/1/2015/fastqc/dc3a1e069_51720e9cf*read2_fastqc.zip
/users/GR/mb/jquilez/projects/didactic_dataset/scripts/utils/io_metadata.sh -m quality_control_raw_reads -s dc3a1e069_51720e9cf -p PE
