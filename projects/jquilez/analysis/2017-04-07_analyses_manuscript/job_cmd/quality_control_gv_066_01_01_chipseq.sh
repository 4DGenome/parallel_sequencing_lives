#!/bin/bash
#$ -N quality_control_gv_066_01_01_chipseq
#$ -q 
#$ -l virtual_free=
#$ -l h_rt=
#$ -o /users/GR/mb/jquilez/projects/conseq/projects/jquilez/analysis/2017-04-07_analyses_manuscript/job_out/quality_control_gv_066_01_01_chipseq_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/projects/conseq/projects/jquilez/analysis/2017-04-07_analyses_manuscript/job_out/quality_control_gv_066_01_01_chipseq_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 
/software/mb/bin/fastqc --extract /gv_066_01_01_chipseq*read1.fastq.gz -o ; rm -f /gv_066_01_01_chipseq*read1_fastqc.zip
/users/GR/mb/jquilez/utils/io_metadata.sh -m quality_control_raw_reads -s gv_066_01_01_chipseq -p
