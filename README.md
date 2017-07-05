# Didactic dataset

In ["Parallel sequencing lives, or what makes large sequencing projects successful"]() ~~update the name of the manuscript and the link~~ we used the life of `b1913e6c1_51720e9cf`, a Hi-C sample, to (i) illustrate the challanges posed by the management and analysis of high-throughput sequencing (HTS) samples and (ii) propose best-practices for doing it successfully. Here we further illustrate these recommendations with a didactic dataset of 4 Hi-C samples (including `b1913e6c1_51720e9cf`).


<br>

## Table of Contents
- [Installation and usage](#installation-and-usage)
- [Dependencies](#dependencies)
- [Metadata](#metadata)
	- [Metadata collection](#metadata-collection)
	- [Metadata SQL database](#metadata-sql-database)
	- [Extract metadata](#extract-metadata)
	- [Print freeze](#print-freeze)
- [Sequencing index concordance](#sequencing-index-concordance)
- [Sample identification](#sample-identification)
- [Structured and hierarchical data organisation](#structured-and-hierarchical-data-organisation)
	- [(1) Raw data](#1-raw-data)
	- [(2) Processed data](#2-processed-data)
	- [(3) Analysis results](#3-analysis-results)
- [Automation of analysis pipelines](#automation-of-analysis-pipelines)
	- [Scalability](#scalability)
	- [Parallelization](#parallelization)
	- [Automatic configuration](#automatic-configuration)
	- [Modularity](#modularity)
- [Interactive web application](#interactive-web-application)
- [Documentation](#documentation)

<br>

## Installation and usage

Download the entire repository with:
```bash
git clone https://github.com/4DGenome/parallel_sequencing_lives.git
```

In many of the scripts in [scripts](https://github.com/4DGenome/parallel_sequencing_lives/tree/master/scripts/) paths are relative to the `$PSL` Unix variable defined at the beginning of the script, which is set to `/users/GR/mb/jquilez/projects/parallel_sequencing_lives` (the absolute path to the repository directory in the machine where it was developed). As an example see:
```bash
head -n 12 scripts/utils/check_sequencing_index_concordance.sh
```

The scripts are written so that they can be executed from the directory where the repository is cloned by conveniently changing the `$PSL` value, which can be achieved for all scripts with:
```bash
TARGET_DIR=my_home_directory
for s in `find scripts -name "*.sh"`; do
	IDIR='\/users\/GR\/mb\/jquilez\/projects\/parallel_sequencing_lives'
	sed -i "s/$IDIR/$TARGET_DIR/g" $s
done
```


<br>

## Dependencies

- bgzip
- FastQC
- java
- Python and its packages:
	- os
	- sys
	- dataset
	- pandas
	- numpy
	- collections
	- glob
	- time
	- datetime
	- hashlib
- samtools
- shasum
- tabix
- trimmomatic
- unzip
- wget


<br>

## Metadata

Sequencing reads (FASTQ files) should not be the only data making a HTS experiment. Metadata describe and provide information about the sequenced DNA sample and are thus required at different steps of the analysis of HTS data. Despite their importance very often metadata are scattered, inaccurate, insufficient or even missing, and that there is a decay in the availability of metadata for older sequencing samples. Factors contributing to this situation include (i) disconnection between the experiment and the analysis (in other words, the experimentalist may not be aware of the information needed for the analysis), (ii) short-term view of the experiment (performed just for a specific ongoing project without considering its potential future use), (iii) the initial investment required for establishing a metadata collection system as well as the subsequent inertia of filling the form, and (iv) high turnover of people. Altogether, this results in a poor description of sequencing samples and can affect performance.

### Metadata collection

In our projects Metadata are collected samplewise via an online [Google Form](https://www.google.com/forms/about/) - a snapshot of the form is shown below.

![google_form](https://github.com/4DGenome/parallel_sequencing_lives/blob/master/figures/google_form.png)

Once the form is completed for one sample, its metadata as well as of those of other submitted samples are available for approved users in the associated online Google spreadsheet. This spreadsheet can be [published as a text file](https://support.google.com/docs/answer/37579?co=GENIE.Platform%3PSLesktop&hl=en), which we use to dump the metadata into a local SQL database. 

In [this table](https://github.com/4DGenome/parallel_sequencing_lives/blob/master/tables/table_features_metadata_collection_system.pdf) we propose desired features for a metadata collection system, and in [this other one](https://github.com/4DGenome/parallel_sequencing_lives/blob/master/tables/table_metadata_fields.xlsx) we describe the metadata fields that are collected for each Hi-C experiment.

Below we provide a script to download the metadata in an [example online text file](https://zenodo.org/record/817549/files/metadata.tsv) and dump them into a SQL database:

```bash
scripts/utils/io_metadata.sh -m download_input_metadata
```
The `-m` command selects the mode. When `download_input_metadata` is passed, the metadata is downloaded from the corresponding URL and added the database. Note that 2 files are generated:
```bash
metadata/metadata.db
metadata/metadata.tsv
```
- `metadata.db`; SQL database with several tables storing the metadata collected for each HTS experiment as well as information resulting from running the analysis pipeline (see below)
- `metadata.tsv`: metadata as a tab-separated values file


### Metadata SQL database

Below is the scheme of the metadata SQL database, each table with its name (top), a subset of its fields and its primary unique key highlighted in orange. The metadata collected by the user are dumped into the `input_metadata` table, which can be associated through the primary key `SAMPLE_ID` to other information generated from the analysis of the sequencing data. For instance, `quality_control_raw_reads` stores information related to the quality of the raw reads reported by [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). As many tables as analysis pipelines can be added to collect parameters used and metrics generated by these (e.g. hic, rnaseq). Because the latter will store information only for the last time a pipeline is executed for a given sample, including a table like `jobs` to keep track of different executions may be useful (e.g. benchmarking of pipelines or parameters).

![metadata_sql_database](https://github.com/4DGenome/parallel_sequencing_lives/blob/master/figures/metadata_sql_database.png)


### Extract metadata

Once the input metadata are stored in the SQL database, they can be accessed with:
```bash
scripts/utils/io_metadata.sh -m get_from_metadata -s b1913e6c1_51720e9cf -t input_metadata -a CELL_TYPE
scripts/utils/io_metadata.sh -m get_from_metadata -s b1913e6c1_51720e9cf -t input_metadata -a SEQUENCING_INDEX
```
- sample (`-s`)
- table in the metadata database (`-t`)
- attribute (`-a`) that is printed out

From which we can confirm that `b1913e6c1_51720e9cf` derives from the T47D cell type and the sequencing index used in the DNA library preparation.


### Print freeze

Although this didactic dataset only includes 4 samples, for projects with a continuous inclusion of samples it may be useful to have a _freeze_ of the metadata at a given moment:
```
scripts/utils/io_metadata.sh -m print_freeze
```
Prints the content of the SQL metadata database (one `*.csv` file per table) into a dated directory:
```
ls -l metadata/freezes/
```


<br>

## Sequencing index concordance

Because of the high yield of current sequencing technologies (in the order of tens to hundreds of gigabase pairs), a common practice is loading multiple biological samples into the sequencing instrument; the DNA from each sample is labeled with a different sequence index so that the reads generated can be assigned to the biological sample from which they originate.

For instance, from the metadata we can get the sequencing index for the 4 samples of the didactic dataset:
```bash
cut -f2,20  metadata/metadata.tsv |column -t
```

Ideally, the sequencing index is also contained in the first line of each sequence read (a [FASTQ file](https://en.wikipedia.org/wiki/FASTQ_format) normally uses four lines per sequence), as shown for `b1913e6c1_51720e9cf` with:
```
zcat data/hic/raw/2015-04-28/b1913e6c1_51720e9cf_read1.fastq.gz |head
```

Therefore, it is a good practice to compare for each sample the sequencing index in the metadata and in the associated FASTQ files, which can be done with:
```bash
scripts/utils/check_sequencing_index_concordance.sh b1913e6c1_51720e9cf
```

This sanity check can be incorporated into the analysis pipelines so that a warning is issued unless both sequencing indexes match (see below).


<br>

## Sample identification

As happened to `T47D_rep2` in our story, often sequencing samples and the associated files (e.g. FASTQ files) are identified with names that _describe_ the HTS experiment and/or are easy to remember for the person who performed it. However, this practice has several undesired consequences. Identical or similar identifiers referring to different HTS experiments as well as vague sample names, especially if there is no associated metadata, can preclude any analysis or lead to errors. Moreover, such unsystematic sample naming undermines the capability to process samples programmatically (e.g. search for data, parallelize scripts), which impairs automation and may also lead to errors.

Therefore, we established a system to generate unique sample identifiers (ID) in an automated manner based on a selection of fields from the metadata that uniquely points to a sequencing experiment. As illustrated below, two sets of either biological or technical fields that unequivocally defined a sequencing sample were identified. Then, for a given sample the values of the biological fields treated as text are concatenated and computationally digested into a 9-mer, and the same procedure is applied to the technical fields. The two 9-mers are combined to form the sample identifier (ID), as happened for b1913e6c1_51720e9cf. 

![sample_id_scheme](https://github.com/4DGenome/parallel_sequencing_lives/blob/master/figures/sample_id_scheme.png)

For example, the sample ID values for the samples of the didactic dataset were generated with:
```bash
scripts/utils/sample_id_generator.sh
```
which prints:
```
[info]	SAMPLE_ID for sample with Timestamp = 8/10/2015 14:13:44 is: b1913e6c1_51720e9cf
[info]	SAMPLE_ID for sample with Timestamp = 8/10/2015 14:35:30 is: b1913e6c1_51720e9cf
[info]	SAMPLE_ID for sample with Timestamp = 8/10/2015 14:38:24 is: dc3a1e069_ec92aa0bb
[info]	SAMPLE_ID for sample with Timestamp = 2/22/2016 12:35:00 is: b7fa2d8db_bfac48760
```

Despite the apparent non-informativeness of this sample ID approach, it easily allows identifying:
- biological replicates, which will share the first 9-mer (e.g. `b1913e6c1_51720e9cf` and `dc3a1e069_ec92aa0bb`)
- experiments generated in the same batch, which will share the second 9-mer (e.g. `b1913e6c1_51720e9cf` and `b1913e6c1_51720e9cf`; `dc3a1e069_ec92aa0bb` was _sequenced_ in the same run, 2015-08-10, but in different sequencing libraries)

While the specific fields used to generate the sample ID can vary, it is important that they unambiguously define a sequencing sample (otherwise duplicated identifiers can emerge) and that they are always combined in the same order to ensure reproducibility. Indeed, another advantage of this naming scheme is that the integrity of the metadata can be checked, as altered metadata values will lead to a different sample ID.


<br>

## Structured and hierarchical data organisation

Collecting metadata and labelling HTS experiments efficiently is useless if data cannot be located. Unfortunately, we have observed that the raw and processed sequencing data as well as the results derived from them tend to be stored in a rather untidy fashion. Such situations may happen if files are stored without anticipating additional sequencing runs and analyses, processed and analysed on the fly or managed by several people with different or missing perceptions of organizing data (the latter is a frequent case when people leave and must be replaced). Difficulties to find data are aggravated by the existence of duplicated sample names and lack of recorded metadata.

Alternatively, we suggest a structured and hierarchical organisation that reflects the way in which sequencing data are generated and analyzed. First, **raw data** (i.e. FASTQ files) are processed sample-wise with relatively standard but tunable analysis pipelines (“core analysis pipelines” in the figure below) which generate a variety of directories and files. In a second step, **processed data** from one or more samples are combined to perform downstream analyses and produce the **analysis results**.

![stages_hts_data](https://github.com/4DGenome/parallel_sequencing_lives/blob/master/figures/stages_hts_data.png)

Considering this, we propose the following considerations for storing the raw data, processed data and analysis results.


### (1) Raw data

In general, experiments are sequenced in different multi-sample runs separated in time, so it is convenient storing the raw data from the sequencing samples grouped by the run in which they were generated:

```bash
tree -C -L 2 data/*/raw
# -C colors directories and files differently
# -L 2 only shows 2 levels of the tree for clarity
```

Dated run directories contain:
- Compressed FASTQ files with the sequencing reads (`*fastq.gz`); note that FASTQ files contain the unique SAMPLE_ID and what read of the paired-end sequencing they contain (for single-end data, we find useful to still use `read1` for consistency)
- information about their quality of the raw reads (e.g. [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports)


### (2) Processed data

Allocate one directory for each sequencing sample:
```
ll data/hic/samples
```

Include subdirectories for the output of the different steps of the pipeline:
```
ll data/hic/samples/b1913e6c1_51720e9cf/plots/hg38_mmtv
```

Include subdirectories for the output of the logs of the programs used:
```
data/hic/samples/b1913e6c1_51720e9cf/logs/b1913e6c1_51720e9cf_trim_reads_trimmomatic_paired_end.log
data/hic/samples/b1913e6c1_51720e9cf/logs/hg38_mmtv/b1913e6c1_51720e9cf_align_and_merge_paired_end.log
```

Include subdirectories for the file integrity verifications:
```
data/hic/samples/b1913e6c1_51720e9cf/checksums/hg38_mmtv/2017-06-29-15-24/files_checksums.sha
```
Everytime the Hi-C pipeline is run on `b1913e6c1_51720e9cf` a time stamped directory (`2017-06-29-15-24`) with a [SHA](https://en.wikipedia.org/wiki/Sha1sum) file is generated. 
```
9ea81ab473cde0bdaa660e06c1c1a07f0b40d2fb  /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/b1913e6c1_51720e9cf_read1.fastq
bb631437d44e57643cd8d931ebdddde7a834db66  /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/raw/2015-04-28/b1913e6c1_51720e9cf_read2.fastq
cff9077b44948ed7be700f5c1d9959f6ee390ceb  /users/GR/mb/jquilez/assemblies/homo_sapiens/hg38_mmtv/ucsc/hg38_mmtv_chr1-22XYM.fa
c30f48fcdb6cd84d9c8da0d7515e5237b18a9e40  /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/samples/b1913e6c1_51720e9cf/results/hg38_mmtv/processed_reads/b1913e6c1_51720e9cf_both_map.tsv
de317ed0e4ffe8ecb082b20600abb72ad5450c64  /users/GR/mb/jquilez/projects/parallel_sequencing_lives/data/hic/samples/b1913e6c1_51720e9cf/results/hg38_mmtv/filtered_reads/b1913e6c1_51720e9cf_filtered_map.tsv
2c0ab11c9ef7caafd1b3cc309cbe778dd100167d  /users/project/4DGenome/data/hic/samples/b1913e6c1_51720e9cf/results/hg38_mmtv/processed_reads/b1913e6c1_51720e9cf_both_map.bam.sam
```

As shown above, the `*.sha` file contains alphanumerical text strings associated to key files used in the pipeline (e.g. input FASTQs, genome reference sequece), which can be used to assess the integtity of such files.

Note that many of the paths above include `hg38_mmtv`. This not only is informative about the version of the human genome used to process the data but also allows to accommodate variations in the analysis pipelines without over-writting existing data. With the corresponding changes in the Hi-C pipeline, `b1913e6c1_51720e9cf` can be re-processed on `hg19` and the output of the pipeline will be stored in _parallel_ paths. This can be extended to other variables such as the aligner used, etc.


### (3) Analysis results

Finding effortlessly what, when and how data are processed is crucial.

When one has to perform an analysis some basic questions arise: for what/who is it? what, when and who is it done?

Firstly, we found convenient allocating a directory for each of the users who requests an analysis, especially when there are many of them. While saving the downstream analyses grouped by projects is also a sensible option, this may be straightforward only for analyses under the umbrella of a well-defined broad project. Moreover, very often the name initially given to a project when its directory is created may become imprecise as the project evolves, analyses are unrelated to any existing project or a given analysis is used in many projects.

For instance, all the analyses generated for the [manuscript]()~~link to manuscript~~ and this didactic dataset are allocated in the `projects/jquilez/analysis` directory, which was generated with:
```bash
scripts/utils/make_project_directory.sh jquilez
```

Each of the analyses is named with a timestamp (which allows chronological sorting of the analyses) plus a descriptive tag from a controlled vocabulary (e.g. "process_hic_samples")
```
projects/jquilez/analysis/2017-04-07_analyses_manuscript
projects/jquilez/analysis/2017-06-29_process_hic_samples
```

And each of the analysis directories includes well-defined subdirectories for data, figures, scripts, etc.

Importantly, this structured and hierarchical organisation of the data facilitates both humand and computer searches. As an example, the command below allows quickly checking that the trimming step of the pipeline completed successfully for all Hi-C samples:
```bash
cat data/hic/*/*/logs/*trim_reads_trimmomatic_paired_end.log |grep "Completed successfully"
``` 


<br>

## Automation of analysis pipelines

Analysing HTS data is hardly ever a one-time task. Samples are often sequenced at different time points so core analysis pipelines have to be executed for every new sequencing batch. Also, samples may need to be re-processed when analysis pipelines are modified substantially (for instance, by including new programs or changing key parameter values) to ensure that data from different samples are comparable in the downstream analysis. At the downstream level, repeating analyses with different datasets or variables, just to name a few variations, is a common task. We therefore identified four desired features for the code used in the analysis.

![automation_pipelines](https://github.com/4DGenome/parallel_sequencing_lives/blob/master/figures/automation_pipelines.png)

### Scalability

Firstly, it needs to be scalable, that is, effortlessly executed for a single sample or for hundreds. For instance, the single command:
```bash
scripts/pipelines/hic-16.05/hic_submit.sh scripts/pipelines/hic-16.05/hic.config
```
launches the [`hic-16.05`](https://github.com/4DGenome/parallel_sequencing_lives/tree/master/scripts/pipelines/hic-16.05) Hi-C pipeline in all the samples in the [`hic.config`](https://github.com/4DGenome/parallel_sequencing_lives/blob/master/scripts/pipelines/hic-16.05/hic.config) file.

### Parallelization

Processing hundreds of samples sequentially is impractical so, secondly, code has to be parallelizable to exploit multi-core computing architectures to process multiple samples simultaneously and speed up the individual steps within the analysis.

[`hic.config`](https://github.com/4DGenome/parallel_sequencing_lives/blob/master/scripts/pipelines/hic-16.05/hic.config) allows configuring whether the samples pipeline scripts:
```
ll scripts/pipelines/hic-16.05/job_cmd
```
are executed sequentially or submitted to the queuing system of a computing cluster as well as the cluster options (memory and time allocated, number of nodes, etc).

_Note that parallelization of [`hic-16.05`](https://github.com/4DGenome/parallel_sequencing_lives/tree/master/scripts/pipelines/hic-16.05) is implemented for [this type] of cluster._

Moreover, the pipeline code in [`hic.sh`](https://github.com/4DGenome/parallel_sequencing_lives/blob/master/scripts/pipelines/hic-16.05/hic.sh) is adapted to use the specified number of nodes, if possible.

### Automatic configuration

Automatic configuration of the variable values such as species or read length is necessary so that these need not to be set for each sample. Each pipeline script retrieves the required variable values from `metadata/metadata.db` through the `scripts/utils/io_metadata.sh` (see above).

### Modularity

The pipeline code in [`hic.sh`](https://github.com/4DGenome/parallel_sequencing_lives/blob/master/scripts/pipelines/hic-16.05/hic.sh) 
is grouped into modules that can be executed all sequentially or individually by specifying it in [`hic.config`](https://github.com/4DGenome/parallel_sequencing_lives/blob/master/scripts/pipelines/hic-16.05/hic.config). See the [README](https://github.com/4DGenome/parallel_sequencing_lives/tree/master/scripts/pipelines/hic-16.05) of `hic-16.05` for more information on the modules available.


<br>

## Interactive web application

HTS analysis workflows generate a great number of files; for instance, Processing `b1913e6c1_51720e9cf` generates many files:
```bash
tree data/hic/samples/b1913e6c1_51720e9cf
```
However, such files can be useless for some users, if for instance they do not have access to the computer where the data are stored, files are too big to be opened with text editors and/or users lack the skills to manipulate them with the right tools (e.g. Unix, [BEDtools](http://bedtools.readthedocs.io/en/latest/), [SAMtools](http://samtools.sourceforge.net/)). Even if this is not the case, as the number of samples increases, better ways to visualize the data than inspecting files individually are essential.

Therefore, we suggest to implement interactive web applications that display the processed data and allow performing specific analyses in a user-friendly manner. As an example, we take advantage of our structured and hierarchical data organisation as well as the available metadata to deploy a web application to visualise processed Hi-C data using [Shiny](https://shiny.rstudio.com/).

![hisee](https://github.com/4DGenome/parallel_sequencing_lives/blob/master/figures/hisee.png)

We recommend defining the specific features of such web applications with their potential users, because implementing them requires effort and attempting to comprehend unnecessary functions may lead to a loss of time.


<br>

## Documentation

From the moment HTS data are generated, they go through several procedures (e.g. compression, alignment, statistical analysis) that will eventually generate results, typically in the form of text, tables and figures. Very often the details of how these are generated are absent or vaguely documented, which may result in little understanding of the results, irreproducibility and hampers the identification of errors.

On the contrary, we recommend to document as much as possible all the parts involved in the analysis and here provide some tips for doing so.

1. Write in README files how and when software and accessory files (e.g. genome reference sequence, annotation) are obtained. As an example:

```
# 2016-01-14: Download hg38 full dataset from UCSC Genome Browser
# -------------------------------------------------------------------------------------

# Download from UCSC's FTP
INDIR="//hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/"
OUTDIR="$HOME/assemblies/homo_sapiens/hg38/ucsc"
mkdir -p $OUTDIR
rsync -a -P rsync:$INDIR $OUTDIR
# Untar and uncompress FASTA files
FASTA=$OUTDIR/chromFa
mkdir -p $FASTA
tar -zxvf $OUTDIR/hg38.chromFa.tar.gz
mv ./chroms/*.fa $FASTA
rm -rf ./chroms
# Concatenate chromosome FASTA files into a single one (only for autosomes plus chrX)
chroms="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"
ofasta=$OUTDIR/hg38.fa
rm -f $ofasta
for c in $chroms; do
        cat $FASTA/$c.fa >> $ofasta
done
# Check the resulting file
cat $ofasta | grep ">"
```

2. Allocate a directory for virtually any task, as shown in:
```
projects/jquilez/analysis/2017-04-07_analyses_manuscript
projects/jquilez/analysis/2017-06-29_process_hic_samples
```

3. Code core analysis pipeline to log the output of the programs and verify files integrity (see [Structured and hierarchical data organisation](#structured-and-hierarchical-data-organisation))e.g. hg19).

4. Document procedures using [Markdown](https://daringfireball.net/projects/markdown/), [Jupyter Notebooks](http://jupyter.org/), [RStudio](https://www.rstudio.com/) or alike.

5. Specify the non-default variable values that are used. For instance, the file documenting the how Hi-C samples were processed:
```
projects/jquilez/analysis/2017-06-29_process_hic_samples/2017-06-29_process_hic_samples.md
```
The file above contains:
- the variable values passed to the script that performed the quality control of the raw reads (`scripts/utils/quality_control.sh`)
- a copy of the configuration file used to run the `hic-16.05` Hi-C pipeline.

In this way there is an exact record of the parameter values used in the analysis.


<br>



