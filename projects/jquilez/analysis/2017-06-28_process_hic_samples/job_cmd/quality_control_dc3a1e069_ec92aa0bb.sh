#!/bin/bash
#$ -N quality_control_dc3a1e069_ec92aa0bb
#$ -q short-sl7
#$ -l virtual_free=15G
#$ -l h_rt=06:00:00
#$ -o /users/GR/mb/jquilez/projects/parallel_sequencing_lives/projects/jquilez/analysis/2017-06-28_process_hic_samples/job_out/quality_control_dc3a1e069_ec92aa0bb_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/projects/parallel_sequencing_lives/projects/jquilez/analysis/2017-06-28_process_hic_samples/job_out/quality_control_dc3a1e069_ec92aa0bb_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 1
/software/mb/bin/fastqc --extract /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/dc3a1e069_ec92aa0bb_read1.fastq.gz -o /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/fastqc; rm -f /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/fastqc/dc3a1e069_ec92aa0bb*read1_fastqc.zip
/software/mb/bin/fastqc --extract /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/dc3a1e069_ec92aa0bb_read2.fastq.gz -o /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/fastqc; rm -f /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/fastqc/dc3a1e069_ec92aa0bb*read2_fastqc.zip
/users/GR/mb/jquilez/projects/parallel_sequencing_lives/scripts/utils/io_metadata.sh -m quality_control_raw_reads -s dc3a1e069_ec92aa0bb -p PE
