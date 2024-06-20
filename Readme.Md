# **Poresippr** (under development)
![ALT Poresippr](https://github.com/madhubioinfo/poresippr_code/blob/main/newimage.png)

Poresippr is a bioinformatics pipeline it automates the process of running Guppy for basecalling and performs downstream analysis such as mapping with Minimap2 and Samtools. The script takes a CSV input file containing necessary parameters and runs in a loop, continuously processing new data ( every 30 minutes).

## Requirements

Python 3.x <br>
Guppy basecaller ( This uses nvidia gpu nodes) <br>
Minimap2 <br>
Samtools <br>
Singularity (optional, when running with singularity container) <br>

## Author

Mathu Malar C
Mathu.Malar@inspection.gc.ca

####Usage

####Singularity
