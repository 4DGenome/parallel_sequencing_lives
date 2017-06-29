# Didactic dataset

In ["Parallel sequencing lives, or what makes large sequencing projects successful"]() ~~update the name of the manuscript and the link~~ we used the life of `dc3a1e069_51720e9cf`, a Hi-C sample, to (i) illustrate the challanges posed by the management and analysis of high-throughput sequencing (HTS) samples and (ii) propose best-practices for doing it successfully. Here we further illustrate these recommendations with a didactic dataset with 4 Hi-C samples (including `dc3a1e069_51720e9cf`).


<br>

## Table of Contents
- [Installation and usage](#installation-and-usage)
- [Metadata](#metadata)
- [Sequencing index concordance](#sequencing-index-concordance)
- [Sample identification](#sample-identification)
- [Structured and hierarchical data organisation](#structured-and-hierarchical-data-organisation)
- [Automation of analysis pipelines](#automation-of-analysis-pipelines)

<br>

## Installation and usage

Download the entire repository with:
```
git clone https://github.com/4DGenome/didactic_dataset.git
```

In many of the scripts in [scripts](https://github.com/4DGenome/didactic_dataset/tree/master/scripts/) paths are relative to the `$DD` Unix variable defined at the beginning of the script, which is set to `/users/GR/mb/jquilez/projects/didactic_dataset` (the absolute path to the repository directory in the machine where it was developed). As an example see:
```
head -n 12 scripts/utils/check_sequencing_index_concordance.sh
```

The scripts are written so that they can be executed from the directory where the repository is cloned by conveniently changing the `$DD` value, which can be achieved for all scripts with:
```
TARGET_DIR=my_home_directory
for s in scripts/*/*.sh; do
	IDIR='\/users\/GR\/mb\/jquilez\/projects\/didactic_dataset'
	sed -i "s/$IDIR/$TARGET_DIR/g" $s
done
```


<br>

## Metadata

Sequencing reads (FASTQ files) should not be the only data making a HTS experiment. Metadata describe and provide information about the sequenced DNA sample and are thus required at different steps of the analysis of HTS data. Despite their importance very often metadata are scattered, inaccurate, insufficient or even missing, and that there is a decay in the availability of metadata for older sequencing samples. Factors contributing to this situation include (i) disconnection between the experiment and the analysis (in other words, the experimentalist may not be aware of the information needed for the analysis), (ii) short-term view of the experiment (performed just for a specific ongoing project without considering its potential future use), (iii) the initial investment required for establishing a metadata collection system as well as the subsequent inertia of filling the form, and (iv) high turnover of people. Altogether, this results in a poor description of sequencing samples and can affect performance.

### Metadata collection

In our projects Metadata are collected via an online [Google Form](https://www.google.com/forms/about/) and stored both online (associated Google Spreadsheet). The latter is also accessible as a text file (as explained [here](https://support.google.com/docs/answer/37579?co=GENIE.Platform%3DDesktop&hl=en)), which we use to dump the metadata into a local SQL database. 

In [this table](https://github.com/4DGenome/didactic_dataset/blob/master/tables/table_features_metadata_collection_system.pdf) we propose desired features for a metadata collection system, and in [this other one](https://github.com/4DGenome/didactic_dataset/blob/master/tables/table_metadata_fields.xlsx) we describe the metadata fields that are collected for each Hi-C experiment.

Below we provide a script to download the metadata in an [example online text file](https://zenodo.org/record/817549/files/metadata.tsv) and dump them into a SQL database:

```
scripts/utils/io_metadata.sh -m download_input_metadata
```
The `-m` command selects the mode. When `download_input_metadata` is passed, the metadata is downloaded from the corresponding URL and added the database. Note that 2 files are generated:
```
metadata/metadata.db
metadata/metadata.tsv
```
- `metadata.db`; SQL database with several tables storing the metadata collected for each HTS experiment as well as information resulting from running the analysis pipeline (see below)
- `metadata.tsv`: metadata as a tab-separated values file


### Metadata SQL database

Below is the scheme of the metadata SQL database, each table with its name (top), a subset of its fields and its primary unique key highlighted in orange. The metadata collected by the user are dumped into the `input_metadata` table, which can be associated through the primary key `SAMPLE_ID` to other information generated from the analysis of the sequencing data. For instance, `quality_control_raw_reads` stores information related to the quality of the raw reads reported by [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). As many tables as analysis pipelines can be added to collect parameters used and metrics generated by these (e.g. hic, rnaseq). Because the latter will store information only for the last time a pipeline is executed for a given sample, including a table like `jobs` to keep track of different executions may be useful (e.g. benchmarking of pipelines or parameters).

![metadata_sql_database](https://github.com/4DGenome/didactic_dataset/blob/master/figures/metadata_sql_database.png)


### Extract metadata

Once the input metadata are stored in the SQL database, they can be accessed with:
```
scripts/utils/io_metadata.sh -m get_from_metadata -s dc3a1e069_51720e9cf -t input_metadata -a CELL_TYPE
scripts/utils/io_metadata.sh -m get_from_metadata -s dc3a1e069_51720e9cf -t input_metadata -a SEQUENCING_INDEX
```
- sample (`-s`)
- table in the metadata database (`-t`)
- attribute (`-a`) that is printed out

From which we can confirm that `dc3a1e069_51720e9cf` derives from the T47D cell type and the sequencing index used in the DNA library preparation.


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

Ideally, the sequencing index is also contained in the first line of each sequence read (a [FASTQ file](https://en.wikipedia.org/wiki/FASTQ_format) normally uses four lines per sequence), as shown for `dc3a1e069_51720e9cf` with:
```
zcat data/hic/raw/2015-04-28/dc3a1e069_51720e9cf_read1.fastq.gz |head
```

Therefore, it is a good practice to compare for each sample the sequencing index in the metadata and in the associated FASTQ files, which can be done with:
```bash
scripts/utils/check_sequencing_index_concordance.sh dc3a1e069_51720e9cf
```

This sanity check can be incorporated into the analysis pipelines so that a warning is issued unless both sequencing indexes match (see below).


<br>

## Sample identification

As happened to `T47D_rep1` in our story, often sequencing samples and the associated files (e.g. FASTQ files) are identified with names that _describe_ the HTS experiment and/or are easy to remember for the person who performed it. However, this practice has several undesired consequences. Identical or similar identifiers referring to different HTS experiments as well as vague sample names, especially if there is no associated metadata, can preclude any analysis or lead to errors. Moreover, such unsystematic sample naming undermines the capability to process samples programmatically (e.g. search for data, parallelize scripts), which impairs automation and may also lead to errors.

Therefore, we established a system to generate unique sample identifiers (ID) in an automated manner based on a selection of fields from the metadata that uniquely points to a sequencing experiment. As illustrated below, two sets of either biological or technical fields that unequivocally defined a sequencing sample were identified. Then, for a given sample the values of the biological fields treated as text are concatenated and computationally digested into a 9-mer, and the same procedure is applied to the technical fields. The two 9-mers are combined to form the sample identifier (ID), as happened for dc3a1e069_51720e9cf. 

![sample_id_scheme](https://github.com/4DGenome/didactic_dataset/blob/master/figures/sample_id_scheme.png)

For example, the sample ID values for the samples of the didactic dataset were generated with:
```
scripts/utils/sample_id_generator.sh
```
which prints:
```
[info]	SAMPLE_ID for sample with Timestamp = 8/10/2015 14:13:44 is: dc3a1e069_51720e9cf
[info]	SAMPLE_ID for sample with Timestamp = 8/10/2015 14:35:30 is: b1913e6c1_51720e9cf
[info]	SAMPLE_ID for sample with Timestamp = 8/10/2015 14:38:24 is: dc3a1e069_ec92aa0bb
[info]	SAMPLE_ID for sample with Timestamp = 2/22/2016 12:35:00 is: b7fa2d8db_bfac48760
```

Despite the apparent non-informativeness of this sample ID approach, it easily allows identifying:
- biological replicates, which will share the first 9-mer (e.g. `dc3a1e069_51720e9cf` and `dc3a1e069_ec92aa0bb`)
- experiments generated in the same batch, which will share the second 9-mer (e.g. `dc3a1e069_51720e9cf` and `b1913e6c1_51720e9cf`; `dc3a1e069_ec92aa0bb` was _sequenced_ in the same run, 2015-08-10, but in different sequencing libraries)

While the specific fields used to generate the sample ID can vary, it is important that they unambiguously define a sequencing sample (otherwise duplicated identifiers can emerge) and that they are always combined in the same order to ensure reproducibility. Indeed, another advantage of this naming scheme is that the integrity of the metadata can be checked, as altered metadata values will lead to a different sample ID.


<br>

## Structured and hierarchical data organisation

Collecting metadata and labelling HTS experiments efficiently is useless if data cannot be located. Unfortunately, we have observed that the raw and processed sequencing data as well as the results derived from them tend to be stored in a rather untidy fashion. Such situations may happen if files are stored without anticipating additional sequencing runs and analyses, processed and analysed on the fly or managed by several people with different or missing perceptions of organizing data (the latter is a frequent case when people leave and must be replaced). Difficulties to find data are aggravated by the existence of duplicated sample names and lack of recorded metadata.

Alternatively, we suggest a structured and hierarchical organisation that reflects the way in which sequencing data are generated and analyzed. 

![stages_hts_data](https://github.com/4DGenome/didactic_dataset/blob/master/figures/stages_hts_data.png)


### Raw data

In general, experiments are sequenced in different multi-sample runs separated in time, so it is convenient storing the raw data from the sequencing samples grouped by the run in which they were generated. 


Sequencing run directories can contain not only the FASTQ files with the sequencing reads but also information about their quality (e.g. FastQC [11] reports). Conversely, we discourage storing herein modified, subsetted or merged FASTQ files to ensure that analyses start off from the very same set of reads. 


### Analysis

**Make project/user directory**

```
scripts/utils/make_project_directory.sh jquilez
```
which generates the directory:
```
ll projects/jquilez/
```

**Make analysis directory***

```
scripts/utils/make_analysis_directory.sh jquilez process_hic_samples
```
which generates the directory:
```
ll projects/jquilez/analysis/2017-06-28_process_hic_samples
```

**Quality control of raw reads**

scripts/utils/quality_control.sh




## Automation of analysis pipelines

Analysing HTS data is hardly ever a one-time task. Samples are often sequenced at different time points (Fig. 1b) so core analysis pipelines have to be executed for every new sequencing batch. Also, samples may need to be re-processed when analysis pipelines are modified substantially (for instance, by including new programs or changing key parameter values) to ensure that data from different samples are comparable in the downstream analysis. At the downstream level, repeating analyses with different datasets or variables, just to name a few variations, is a common task. We therefore identified four desired features for the code used in the analysis.
