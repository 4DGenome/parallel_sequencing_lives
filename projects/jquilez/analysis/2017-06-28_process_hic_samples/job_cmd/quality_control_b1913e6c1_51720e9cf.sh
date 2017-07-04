#!/bin/bash
#$ -N quality_control_b1913e6c1_51720e9cf
#$ -q short-sl7
#$ -l virtual_free=15G
#$ -l h_rt=06:00:00
#$ -o /users/GR/mb/jquilez/projects/parallel_sequencing_lives/projects/jquilez/analysis/2017-06-28_process_hic_samples/job_out/quality_control_b1913e6c1_51720e9cf_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/projects/parallel_sequencing_lives/projects/jquilez/analysis/2017-06-28_process_hic_samples/job_out/quality_control_b1913e6c1_51720e9cf_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 1
/software/mb/bin/fastqc --extract /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/b1913e6c1_51720e9cf_read1.fastq.gz -o /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/fastqc; rm -f /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/fastqc/b1913e6c1_51720e9cf*read1_fastqc.zip
/software/mb/bin/fastqc --extract /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/b1913e6c1_51720e9cf_read2.fastq.gz -o /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/fastqc; rm -f /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/fastqc/b1913e6c1_51720e9cf*read2_fastqc.zip
/users/GR/mb/jquilez/projects/parallel_sequencing_lives/scripts/utils/io_metadata.sh -m quality_control_raw_reads -s b1913e6c1_51720e9cf -p PE
