# **Poresippr** (under development)
![ALT Poresippr](https://github.com/madhubioinfo/poresippr_code/blob/main/gitimage.png)

Poresippr is a bioinformatics pipeline it automates the process of running Guppy for basecalling and performs downstream analysis such as mapping with Minimap2 and Samtools. The script takes a CSV input file containing necessary parameters and runs in a loop, continuously processing new data ( every 30 minutes).

## Requirements

Python 3.x
Guppy basecaller ( This uses nvidia gpu nodes)
Minimap2
Samtools
Singularity (optional, when running with singularity container)

##Author

Mathu Malar C
Mathu.Malar@inspection.gc.ca

####Usage

####Singularity
