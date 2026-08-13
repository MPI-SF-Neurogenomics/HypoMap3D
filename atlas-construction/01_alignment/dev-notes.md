# 01_alignment

The goal of this module is download raw fastq data, format them ready for processing, and finally running the alignment.

# 1. Downloading data

First, download the raw data into a structure sufficing <DATADIR>/1a_source/<dataset_name>.
Here, <DATADIR> is the path to the data directory specified in config.toml, 
and <dataset_name> is the name of the dataset according to our metadata file, e.g.,

<DATADIR>/1a_source/SC1-01

The goal is to have a flat structure of fastq.gz files.

To achieve this, some manual work will be required.
There are some useful helper scripts outlined in the following section.

## Noteworthy datasets warranting special care

### SC1-07

For this dataset there was no raw data available. 

Hence, for this dataset we skipped the alignment and cellbender steps, and instead downloaded it into the 2a_cellbender/wo_empty_droplets folder, 
and manually changed the file names to match the file name pattern SC1-07_<sample_name>.h5

### SC2-04

This dataset has to be downloaded in the form of 1 fastq file per cell. For alignment, a manifest file with all fastqs listed was passed to STARsolo, treating all cells to come from a single sample.
Make sure that the structure of the output follows DATADIR/1c_aligned/<dataset_name>/<sample_name>/ as it might break other scripts otherwise. We used sample name = dataset name = SC2-04.

# 2. Manual massaging of your downloaded files

The steps you have to take here depend on the dataset in question.

After massaging your downloaded data your data directory should look like this:

dataset_id/fastq_gz_of_sample.fastq.gz

e.g.

SC2-22/SRR31177245_1.fastq.gz # SRR pattern
SC2-22/Gabrb3_P65_M_mut_L003_R1_001.fastq.gz # BCL2-like pattern

The names of the files can be either in SRR format or can contain a bcl2-like pattern.

In the following, a few useful helpers will be shown.
Some example datasets where the helpers hat to be applied are shown further down.

### untar balling

If the samples of the dataset comes in .tar format.

> ./hlpr_untar_fastq_files.sh SC2-13

will untar all samples in the SC2-13 dir, assuming flat file structure.

### gzipping .fastq files into .fastq.gz (multiprocessing)

Some files come in unzipped format. 

Gzip them using 

> ./hlpr_gzip_files.sh SC2-25

logs to SC2-25/logs in that case.

Again, assuming flat file structure.

### bam to fastq

Some files have to be transformed from bam to fastq format.

> ./hlpr_run_bamtofastq.sh SC2-24 --jobs 4

will do so using multiprocessing and 4 workers. 
Omitting the jobs argument will use all available workers.

### Moving files of subdirectories into parent directory

This creates flat file structures in case your the downloaded data comes in the form of e.g.

SC2-13/<sample_name>/some_file.file_ending 

> ./hlpr_move_subdir_files_into_target_dir.sh SC2-13

will create a flat file structure for dataset SC2-13.

there are some additional arguments:

--dry-run: allows you to check what would happen to the filenames during the move. omit that flag to actually perform the operation.

--prepend-subdir: will prepend the subdir names to the filenames found in the subdirs, .e.g, would move <sample_name>_some_file.file_ending into the SC2-13 dir in the above example

--delete-empty-subdirs will remove all empty subdirs, e.g. move the files out of the <sample_name> dir, then delete that dir

example:

> ./hlpr_move_subdir_files_into_target_dir.sh SC2-13 --dry-run --prepend-subdir

### Removing trailing .1 from fastq.gz.1 files

Some datasets come with a ".1" appended to their fastq.gz file names.

This will crash some scripts later.

The following script gets rid of the ".1". 

Renames the files INPLACE.

> ./hlpr_remove_trailing_.1_from_fastqgz.sh SC2-23

You can do a dry run to see what would happen to the files without actually renaming anything by passing the --dry-run or -n flag

> ./remove_trailing_.1_from_fastqgz.sh /path/to/directory [--dry-run|-n]

# 3. Harmonising names to a standardized bcl2 format 

For this there is the 0d_bcl2_formatting.sh script.

It will create symlinks of harmonised file names with the biosample/sample_id mapped, and the bcl2 pattern determined/harmonised into <DATADIR>/1b_bcl_FASTQ

The result is a structure where we have symlinks in the format <DATADIR>/1b_bcl_FASTQ/<dataset_name>/<dataset_name>_<BioSampleID>_<bcl_string>.fastq.gz 

Here, dataset_name can be a string that has a "-" in it, such as sringA-stringB.
BioSampleID is everything up to the bcl string. BioSampleID can be a collection of strings divided by "_". afterwards, there is the bcl string. 

Example: SC2-13_M2TX_201001_203_D01_S01_L003_I1_001.fastq.gz

dataset_name: SC2-13

BioSampleID: M2TX_201001_203_D01

bcl_string: S01_L003_I1_001

The script processes data depending on their names:
- files with SRRXXXXXX_X.fastq.gz format
- files with bcl2-like format (e.g., ..._LXX_RX_XXX.fastq.gz)

If the data is in the bcl2-like format you may have to provide some metadata.

## SRRXXXXXX_X.fastq.gz format

If the data comes in the SRRXXXXXX_X.fastq.gz format, the bcl2 pattern will be determined automatically.

## bcl2-like format

Files with bcl2-like format such as Cnga2_P65_M_ctrl_1_L002_R2_001.fastq or AW20221001_RF_0_1_H3TGHBGXN_bamtofastq_S1_L001_R2_001.fastq.gz will be detected and harmonised, and mapped to their BioSample/sample_id afterwards.

For such file names you have to ensure that there is a column in your metadata file, i.e., either SraRunTableMetaData.csv or nemo manifest, that can be used to map the file to the BioSample/sample_id.

You may have to add this information manually if it is not provided in some column already. 

If you are indeed unlucky, add a column into the metadata file which has the information of the file name up to the bcl2 part, e.g., AW20221001_RF_0_1_H3TGHBGXN_bamtofastq should be in the same row as its corresponding BioSample (if SRARunTableMetaData) or sample_id (if NeMO Manifest) for the file AW20221001_RF_0_1_H3TGHBGXN_bamtofastq_S1_L001_R2_001.fastq.gz. 

The name of the column does not really matter as the metadata file present in the datadir will be searched for a column that contains the filename (or a substring of it) and use that column to map the filename to the sample_id/biosample depending on whether the metadata file is a SraRunTableMetaData file or a Nemo manifest.  

Further notes:
- bcl2 patterns will be harmonised, e.g., replacing "." with "_" between the S L R etc.
- the bcl2 pattern will also be detected if the S information from the bcl2 string is missing (e.g. dataset SC2-22: Cnga2_P65_M_ctrl_1_L002_R2_001.fastq). in those cases, S will be set to S01.

Script execution:

From within repo:
> ./hlpr_bcl2_formatting.sh SC2-23

for a single dataset or 

> ./0d_bcl2_formatting_all_datasets.sh 

for all of them at the same time.

Either way, will create logs at 1b_bcl_FASTQ/logs.


# Examples

## History of downloading and processing some of the datasets

### SC2-13

Files had to be untared

> ./hlpr_untar_fastq_files.sh SC2-13

and then moved from within their subdirectories into SC2-13 folder

> ./hlpr_move_subdir_files_into_target_dir.sh SC2-13 --delete-empty-subdirs

### SC2-22

Trailing .1 had to be removed, and a column for mapping the bcl2-like dataset names to BioSamples had to be added.

Also, SRR31177245 and SRR31177284 had to be gzipped first.

SRR31177245 = SAMN44516421 apparently lacks the actual R2 file when downloading SRR file as cell ranger throws errors that reads are too short.
hence, it must be downloaded using [SRA_tools](https://github.com/ncbi/sra-tools/wiki). 


### SC2-24

The data for the libraries GSM8524882 - GSM85248825 ("batch 2") were downloaded as bam files, e.g., from [here](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR30725023&display=data-access).

The other sample GSM8524881 ("batch 1") was available as fastq.gz

bam files had to be transferred from bam to fastq

> ./hlpr_run_bamtofastq.sh SC2-24 --jobs 4

Then, all created fastq files in the bam subdirectories had to be moved one up while prepending the subdir name

> ./hlpr_move_subdir_files_into_target_dir.sh SC2-24/Batch2_RF_possorted_genome_bam_fastq --prepend-subdir --delete-empty-subdirs 

after which all files were collected in the SC2-24 dir without prepending the subdir names

> ./hlpr_move_subdir_files_into_target_dir.sh SC2-24 --delete-empty-subdirs

To become able to map BioSample information, a new column "bamtofastq_filename" was added to SC2-24_SRARunTableMetaData that could be used to map filenames emerging from the above, such as AW20221001_H2B_MCHERRY_cDNA_0_1_H3TGHBGXN_bamtofastq to their BioSample, e.g., SAMN43828887.

Downloaded files that were already fastq files were from one BioSample, and came in 4 sets of 3 files, each set for a different lane.
Lane information could not be inferred, and was thus added to the "lane" column in SC2-24_SRARunTableMetaData.
Lane information was obtained by manually checking lane information from the place where data was downloaded, i.e., [here](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR30725026&display=data-access).
For other samples, lane information was calculated on the fly by ordering Run Ids in ascending order and setting lanes accordingly.




