# 2017-06-29_process_hic_samples
----------------------------------------------------------------------------------------------------

**objective: process Hi-C samples**

**paths are relative to /users/GR/mb/jquilez/projects/parallel_sequencing_lives**



## [2017-06-29] Run FastQC
----------------------------------------------------------------------------------------------------

```
# variables
#samples="dc3a1e069_51720e9cf b1913e6c1_51720e9cf dc3a1e069_ec92aa0bb b7fa2d8db_bfac48760"
#process=quality_control
#project=jquilez
#analysis=2017-06-28_process_hic_samples
#submit_to_cluster=no
#integrate_metadata="yes"
scripts/utils/quality_control.sh
```



## [2017-06-29] Run Hi-C pipeline
----------------------------------------------------------------------------------------------------

Configuration file used:
```
; This configuration file follows the INI file format (https://en.wikipedia.org/wiki/INI_file)

[data_type]
data_type			= hic

[samples]
samples				=dc3a1e069_51720e9cf b1913e6c1_51720e9cf dc3a1e069_ec92aa0bb b7fa2d8db_bfac48760				; e.g.: `samples=s1 s2 s3`, use 'test1 test2' for testing purposes

[pipeline]
pipeline_name		= hic
pipeline_version	= 16.05
pipeline_run_mode	= full_no_downstream_bam_no_dekker_call

[IO mode]
io_mode				= standard									; standard = ******, custom = in and out dir specified
CUSTOM_IN			= scripts/pipelines/hic-16.05/test 		; only used if pipeline_io_mode=custom
CUSTOM_OUT			= scripts/pipelines/hic-16.05/test		; only used if pipeline_io_mode=custom
sample_to_fastqs	= sample_to_fastqs.txt				; file with paths, relative to CUSTOM_IN, to read1 (and read2) FASTQs, only used if pipeline_io_mode=custom

[cluster options]
submit_to_cluster	= no					; the following are only applied if submit_to_cluster=yes
queue				= long-sl7				; for real data = long-sl65, for test = short-sl65
memory				= 100G					; for real data = 100G, for test = 20G
max_time			= 100:00:00 			; for real data = 100:00:00, for test = 6:00:00
slots				= 10 					; for real data = 10, for test = 1
email				= javier.quilez@crg.eu	; email to which start/end/error emails are sent

[metadata]
integrate_metadata	= yes					; yes: metadata is stored into database

[genome]
species					= 					; required if integrate_metadata=no, otherwise, ignored and retrieved from the metadata
version					= 					; required if integrate_metadata=no, otherwise, ignored: hg38_mmtv (homo_sapiens), mm10 (mus_musculus) and dm6(droshophila_melanogaster)
read_length				=  					; required if integrate_metadata=no, otherwise, ignored and retrieved from the metadata

[trimmomatic]
; for recommended values see http://www.broadinstitute.org/cancer/software/genepattern/modules/docs/Trimmomatic/
; and those used in the supplementary data of the Trimmomatic paper (Bolger et al. 2014)
sequencing_type			= PE					; PE=paired-end, SE=single-end
seedMismatches			= 2
palindromeClipThreshold	= 30
simpleClipThreshold		= 12
leading					= 3
trailing				= 3
minAdapterLength		= 1
keepBothReads			= true
minQual					= 3
strictness				= 0.999
minLength				= 36

[tadbit]
restriction_enzyme		= 				; required if integrate_metadata=no, otherwise, retrieved from the metadata
max_molecule_length		= 500
max_frag_size			= 10000
min_frag_size			= 50
over_represented		= 0.005
re_proximity			= 4
reads_number_qc			= 100000
genomic_coverage_resolution	= Mb
frag_map				= True

[downstream]
flag_excluded			= 783
flag_included			= 0
flag_perzero			= 99
resolution_tad			= 50000 				; in bp
resolution_ab			= 100000				; in bp
pis						= 500000				; Dekker/Crane parameter
pids					= 250000				; Dekker/Crane parameter
pnt						= 0.1					; Dekker/Crane parameter
```

Execute pipeline:
```

```


