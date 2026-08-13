# 01_alignment

## 0 series: Preparing input for alignment

### Reference genome 

Here, the index genome files are downloaded and GM28040 is removed, after which cellranger/STARsolo are used to construct the required reference genome files for cellranger/STARsolo respectively

> ./0a_download_index_genome_files.sh

> ./0b_create_cellranger_index_genome.sh

> ./0c_create_STAR_index_genome.sh

__Note: Do not run 0b and 0c in parallel! This will lead to errors.__

### Harmonising fastq file names required for input to alignment scripts

The alignment scripts work on fastq.gz files following a standardized naming convention, which necessitates manual steps.

More details on this in the dev-notes.md file. 

After downloading the data and wrangling it into the shape and structure outlined in the dev-notes, run the following script on each of the datasets, e.g., 

> ./0d_bcl2_formatting.sh SC2-13

## 1 series: Running alignment

## Processing Differences

Difference chemistries will be processed differently.

### Technology Type
1. 10X protocols
	- 10X have their own software for alignment called CellRanger
2. DropSeq protocols
	- STARsolo is used
3. SmartSeq protocols
	- STARsolo is also used but using a different 'mode'

### Suspension Type

Be very careful in how you process single cell vs single nuclear datasets.
Single nuclear datasets will capture pre-mRNA also as such introns must be included.
Single cell datasets can safely exclude introns (and is encouraged).


### 1a series: 10x, cellranger

specify whether it is sn or sc RNAseq; and substitute "auto" with "multiome" if you have samples with multiome chemistry, e.g., SC2-22

single samples can be aligned directly, providing positional arguments: 
1) dataset_dir
2) sample name
3) sequencing type: sn | sc
4) chemistry type: auto | multiome

> ./1a_10x.sh \
	SC2-23 \
	SC2-23_SAMN34591640 \
	sn \
	auto

or with slurm:

(first set your email in ../../config.toml)

> ./1aa_submit_one_sample_10x.sh \
	SC2-23 \
	SC2-23_SAMN34591640 \
	sn \
	auto


you can also submit slurm jobs for all samples of a dataset:

> ./1ab_submit_all_samples_10x.sh \
	SC2-23 \
	sn \
	auto

### 1b series: smartseq

not automated yet

### 1c series: dropseq

not automated yet

## 2 series: alignment QC metrics 

> 2_run_combine_alignment_qc_metrics_to_csv.sh 

combines qc metrics from each sample into a combined csv file from which the estimated number of cells is used in the next module as input for cellbender.

