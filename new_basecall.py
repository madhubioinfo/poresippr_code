import os
import subprocess
import csv
import sys

#####Author Mathu Malar C Mathu.Malar@inspection.gc.ca######
def run_command(command):
    try:
        process = subprocess.run(command, shell=True, check=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {command}")
        print(f"Error message: {e}")

def concatenate_fastq_files(barcode_path, output_file):
    # Listing all FASTQ files in the barcode directory
    fastq_files = [f for f in os.listdir(barcode_path) if f.endswith('.fastq')]
    # Sorting the list to ensure consistent order
    fastq_files.sort()
    # Opening the output file for writing
    with open(output_file, 'w') as outfile:
        # Concatenate all FASTQ files into the output file
        for fastq_file in fastq_files:
            with open(os.path.join(barcode_path, fastq_file), 'r') as infile:
                outfile.write(infile.read())

# Reading the inputs from CSV file
if len(sys.argv) != 2:
    print("Usage: python script.py <input.csv>")
    sys.exit(1)

csv_file_path = sys.argv[1]  # it looks for input.csv file always, change it if you want
with open(csv_file_path, 'r') as csvfile:
    csv_reader = csv.DictReader(csvfile)
    for row in csv_reader:
        reference = row['reference']
        fast5_dir = row['fast5_dir']
        output_dir = row['output_dir']
        config = row['config']
        barcode = row['barcode']
        barcode_values = [int(x) for x in row['barcode_values'].split(',')]
       
        # Creating the output directory if it doesn't exist
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            print(f"Error creating output directory: {output_dir}")
            print(f"Error message: {e}")

        # Runing guppy fast basecalling
        guppy_command = f"guppy_basecaller --input_path {fast5_dir} --save_path {output_dir} --config {config} --barcode_kits {barcode} -x auto -r"
        print(f"Running guppy command please sit back: {guppy_command}")
        run_command(guppy_command)

        # Processing the each provided barcode number by user as an argument in csv file , enclosed in double quotes
        for barcode_num in barcode_values:
            barcode_dir = f"barcode{barcode_num:02d}"  # Forming and looking for the barcode directory name
            barcode_path = os.path.join(output_dir, "pass", barcode_dir)
            if os.path.isdir(barcode_path):
                # Concatenating FASTQ files into a single fastq file
                concatenated_fastq_file = os.path.join(output_dir, f"{barcode_dir}.fastq")
                concatenate_fastq_files(barcode_path, concatenated_fastq_file)

                # Aligning using Minimap2 and create sorted BAM file
                bam_file = os.path.join(output_dir, f"{barcode_dir}_sorted.bam")
                minimap2_command = f"minimap2 -ax map-ont {reference} {concatenated_fastq_file} | samtools view -@ 5 -bS - | samtools sort -o {bam_file} -"
                print(f"Running Minimap2 command: {minimap2_command}")
                run_command(minimap2_command)

                # Calculating coverage and sorting the output from samtools out
                csv_file = os.path.join(output_dir, f"{barcode_dir}_coverage.csv")
                samtools_command = f"samtools coverage {bam_file} | cut -f 1,4 | awk '$2 > 0' | sort -rnk 2,2 > {csv_file}"
                print(f"Running Samtools command sit tight: {samtools_command}")
                run_command(samtools_command)

                # Adding header to the output CSV file
                with open(csv_file, 'r+') as f:
                    content = f.read()
                    f.seek(0, 0)
                    f.write("Stx_type,number_of_reads_mapped\n" + content)

                # Deleting the temporary concatenated FASTQ files else it will throw memory error
                os.remove(concatenated_fastq_file)
            else:
                print(f"Oops no barcode directory found for {barcode_dir}")
