# **Poresippr** (under development)
![ALT Poresippr](https://github.com/madhubioinfo/poresippr_code/blob/main/newimage.png)

Poresippr is a bioinformatics pipeline it automates the process of running Guppy for basecalling and performs downstream analysis such as mapping with Minimap2 and Samtools. The script takes a CSV input file containing necessary parameters and runs in a loop, continuously processing new data ( every 30 minutes).

## Requirements

Python 3.x <br>
Guppy basecaller ( This uses nvidia gpu nodes) <br>
Minimap2 <br>
samtools <br>
Singularity (optional, when running with singularity container) <br>

## Installation requires ubuntu > 22.0

- **Python 3.0** 
    sudo apt install python3.12

- **Minimap2**
    Install minimap2 and installation instructions from here https://github.com/lh3/minimap2

- **samtools**
    Install samtools from here https://www.htslib.org/download/

- **Guppy gpu installation**
    Installation instructions are found here https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/linux-guppy

 - **poresippr singularity container**
     If you don't want to go through the above installation instruction , you can pull the singularity container by singularity pull **  for this you need to have singularity istalled from here https://docs.sylabs.io/guides/3.0/user-guide/installation.html

## Author

Mathu Malar C
Mathu.Malar@inspection.gc.ca

##Introduction

Poresippr is a python workflow which is run while the nanopore sequencing is in progress. This program is designed to do gpu guppy fast basecalling , barcoding and do the mapping with minimap2 with stx database and calculating the coverage and detecting stx genes. This program is used for rapid detection of virulence genes. The database we use is curated stx genes for classifying the stec types. Also this program run every 30 minutes, until the new nanopore data is written. The poresippr database is available here xxx. 

- **guppy basecalling**: Currently this program uses Nanopore Flowcell version R10.4.1 and  Barcode kit was used SQK-RBK114.24. For running guppy nvidia Gpu nodes are required. It also uses the dna_r10.4.1_e8.2_260bps_fast.cfg settings which usually comes with the guppy installation.

- **Read Mapping**: Minimap2 was used for mapping nanopore reads with poresippr database XXX  and detects the gene targets using samtools. It ouputs the csv file for every iteration of 30 minutes containing gene targets and genome coverage for the sample. 

##Input requirements

An input CSV file containing the following information is needed. It should include:
- The location of the Poresippr database (reference)
- The path for nanopore FAST5 reads
- The path where output files, such as CSV and barcoded reads, will be saved
- The appropriate model/config required for basecalling
- The barcode kit used during sequencing
- The barcode numbers

`input.csv`:

```csv
reference,fast5_dir,output_dir,config,barcode,barcode_values
/home/olcbio/PoreSippRDB_240509.fasta,/home/olcbio/240125_MC26299/test,/home/olcbio/240125_MC26299/test_out,dna_r10.4.1_e8.2_260bps_fast.cfg,SQK-RBK114-24,06,07,08,09,10,11
```

##Running the code




