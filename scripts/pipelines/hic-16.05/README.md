# README

**Pipeline to process HiC data**


<br>

## Modules

The pipeline is broken down into modules:

1. `preliminary_checks`
	- check that a sample with the provided SAMPLE_ID exists in the metadata
	- check FASTQ files exist
	- save a SHA checksums of the FASTQ files
	- retrieve information from the sequencer, and store it in the metadata; a warning will be produced if the sequencing index introduced in the metadata does not agree with that seen in the first read of the FASTQ
	- compare the read length as seen in the FASTQ with that introduced either in the metadata or the configuration file 
	- check the FASTQ for the reference genome sequence exists
2. `raw_fastqs_quality_plots`
	- generate quality plots of the raw reads (read1 and read2 separately)
	- extract the percentage of dangling-ends and ligated sites for read1 and read2
3. `trim_reads_trimmomatic`
	- trim sequencing adapters (by default, Illumina's TruSeq3, which are typically used for the HiSeq and NextSeq sequencers) and low-quality ends from the reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
	- extract the number of trimmed reads
4. `align_and_merge`
	- align reads (read1 and read2 separately) to the genome reference sequence (using [GEM](http://www.nature.com/nmeth/journal/v9/n12/full/nmeth.2221.html)), process reads according to the restriction enzyme and merge alignments from read1 and read2
	- save a SHA checksums of the alignments (MAP format) in which both reads align
5. `post_mapping_statistics`
	- generate plots from the alignments: decay of interaction counts with genomic distance, distribution of dangling ends lengths, genomic covarage and interaction matrix of mapped reads in 1-Mb bins and proportion of mapped reads
	- extract metrics: fraction of read1 and read2 alignments, fraction of reads in which both read1 and read2 are mapped and the counts-to-distance slope
6. `reads_filtering`
	- filter reads based on multiple quality parameters; as of 2016-05-17 reads that meet these filters are excluded: self-circle, dangling-end, error, extra dangling-ends, duplicated and random breaks (for more info see [TADbit](http://3dgenomes.github.io/TADbit/tutorial/tutorial_0_mapping.html))
	- save a SHA checksums of the filtered alignments (*tsv file)
7. `post_filtering_statistics`
	- generate plots after filtering reads: genomic coverage and interaction matrix of filtered reads in 1-Mb; also, 1-Mb genomic coverage of dangling ends and self circle reads
8. `map_to_bam`:
	- convert MAP file with the alignments to BAM format
9. `downstream_bam`
	- perform several downstream analysis using the generated BAM (e.g. add SAM-like flags, find TADs and A/B compartments)
10. `dekker_call`
	- calculated insulation index and TAD border using Dekker's method
11. `clean_up`
	- delete (relatively large) intermediate files


The modules can be executed altogether or individually (see [Configuration file](#configuration-file)). The diagram below shows the order in which modules are sequentially executed (numbers), when the full pipeline is run, and the dependencies between modules in case they want to be run individually (e.g. all modules require that `trim_reads_trimmomatic` has been executed):

![hic-16.05](https://github.com/4DGenome/conseq/blob/master/docs/figures_github_repo/hic-16.05/hic-16.05.001.png)


<br>

## Scripts

- `hic.sh`: script with the code of the pipeline
- `hic_submit`: wrapper script that both:
	- retrieves configuration variables and parameter values from the `hic.config` file
	- (if applies) submits jobs (one per sample) to execute the pipeline in a Univa Grid Engine HPC cluster 
- `hic.config`:  configuration file with the list of samples and the hard-coded parameter values (see [Configuration file](#configuration-file))


<br>

## Pipeline execution

```
hic_submit.sh hic.config
```


<br>

## Configuration file

```
; This configuration file follows the INI file format (https://en.wikipedia.org/wiki/INI_file)

[data_type]
data_type			= hic

[samples]
samples				=66950b082_c478f1d09 ad1a9f5b0_c478f1d09 0f24b004c_95a8cd511 0f24b004c_0e31e17c7				; e.g.: `samples=s1 s2 s3`, use 'test1 test2' for testing purposes

[pipeline]
pipeline_name		= hic
pipeline_version	= 16.05
pipeline_run_mode	= full

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

In `pipeline_run_mode`, specifying the name of the module runs that module. Also:
- `full`: runs all modules sequentially
- `full_no_clean_up`: intermediate files are not deleted (for testing purposes; note that large files will be left on disk)
- `full_no_downstream_bam_no_dekker_call`: same as `full` but `downstream_bam` and `dekker_call` modules are not executed (there is no point in doing so for a small number of reads, e.g. 1,000, as it is the case of this didactic dataset)


<br>

## `io_mode` and `integrate_metadata`

When `io_mode = standard`, FASTQ input file(s) and output directory are pre-defined. This behaviour can be changed with `io_mode = custom`, where the `hic.config` file provides the path to the input FASTQs and a file listing them as well as the path to the output directory. As this is a non-standard usage of the pipeline, `integrate_metadata` is internally set to **no** so values for all the variables in the `hic.config` are required. To try the custom mode one can use the data in the `test` directory.


<br>

## Dependencies and edits


- `export PATH="/software/mb/el7.2/anaconda2/bin:$PATH"` makes that a specific version of Python and TADbit are used so it should be changed accordingly
- edit `fasta` so that it points to the genome assembly reference sequence
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- edit `$ADAPTERS` so that it points to the `adapters` subdirectory of the downloaded version of Trimmomatic
- `shasum`
- `bgzip`
- `tabix`
- [Python](https://www.python.org/)
- [SAMtools](http://samtools.sourceforge.net/)

