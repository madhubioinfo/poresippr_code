# **Poresippr** 
![ALT Poresippr](https://github.com/madhubioinfo/poresippr_code/blob/main/newimage.png)

Poresippr is a bioinformatics pipeline that automates the process of running Guppy for basecalling, barcoding and performs downstream analysis such as mapping with Minimap2 and Samtools. The script takes a CSV input file containing the necessary parameters and runs in a loop, continuously processing new data every 30 minutes.

## Requirements

Python 3.x <br>
Guppy basecaller ( This uses nvidia gpu nodes) 6.5.7 <br>
Minimap2 2.24-r1122 <br>
samtools 1.13 <br>
Singularity (optional, when running with singularity container) 4.1.0 <br>

## Installation requires ubuntu > 22.0

- **Python 3.0** 
    sudo apt install python 3.12

- **Minimap2**
    Install minimap2 and installation instructions from here https://github.com/lh3/minimap2

- **samtools**
    Install samtools from here https://www.htslib.org/download/

- **Guppy gpu installation**
    Installation instructions are found here https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/linux-guppy

 - **poresippr singularity container**
     If you don't want to go through the above installation instructions, you can download the Singularity container from here https://zenodo.org/records/12628468. For this, you need to have Singularity installed. You can install it from [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).


## Author

Mathu Malar C
Mathu.Malar@inspection.gc.ca

## Introduction

Poresippr is a Python workflow designed to run while nanopore sequencing is in progress. This program performs GPU-accelerated Guppy basecalling, barcoding, and mapping with Minimap2 using the stx database. It calculates coverage and detects stx gene targets, facilitating the rapid detection of virulence genes. The database used is curated for stx genes to classify STEC types. Additionally, this program runs every 30 minutes, processing new nanopore data as it becomes available. The Poresippr database is not available for public usage yet.

- **guppy basecalling**: Currently, this program uses the Nanopore Flowcell version R10.4.1 and the barcode kit SQK-RBK114.24. Running Guppy requires NVIDIA GPU nodes. The program also uses the `dna_r10.4.1_e8.2_260bps_fast.cfg` settings, which typically come with the Guppy installation.

- **Read Mapping**: Minimap2 is used for mapping nanopore reads with the Poresippr database and detects gene targets using samtools. It outputs a CSV file every 30 minutes containing gene targets and genome coverage for each of the sample.
. 

## Input requirements

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
/home/olcbio/PoreSippRDB_240509.fasta,/home/olcbio/240125_MC26299/test,/home/olcbio/240125_MC26299/test_out,dna_r10.4.1_e8.2_260bps_fast.cfg,SQK-RBK114-24,"06,07,08,09,10,11"
```

> [!NOTE]
>Ensure that the barcode numbers are enclosed in double quotes (") as shown in the barcode_values field. This is necessary to correctly parse multiple barcode values.

## __Instructions to Run the Code__

To run the Python script, use the following command:

```sh
git clone https://github.com/madhubioinfo/poresippr_code.git
cd poresippr_code
python poresippr_basecall_scheduler.py input.csv

To build the singularity container and run from here use the below command [ need to have sudo access and GPU]
sudo singularity build poresippr_container.sif cuda_singularity.def
singularity run --nv poresippr_container.sif input.csv

To download and run the script as a Singularity container, use the following command still need GPU:
click and download the zipped folder containing code https://zenodo.org/records/12628468
singularity run --nv poresippr_container.sif input.csv
```

## Output

**The pipeline will output the following files:**

```
test_out/
├── fail
└── pass
    ├── barcode01
    ├── barcode02
    ├── barcode..
    └── unclassified
├── barcode01_iteration1.csv
├── barcode01_iteration2.csv
└── barcode..XXX.csv
```


## citations

If you use the code for your research please use https://doi.org/10.5281/zenodo.12628468 for citation
