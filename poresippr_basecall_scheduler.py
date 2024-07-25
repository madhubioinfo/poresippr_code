import os
import subprocess
import csv
import sys
import time
import argparse
import signal
from multiprocessing import Process, Value

##### Author Mathu Malar C Mathu.Malar@inspection.gc.ca ######

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

def print_usage():
    if "SINGULARITY_NAME" in os.environ:
        print("Usage: singularity run --nv mycontainer.sif <input.csv> <metadata.csv>")
    else:
        print("Usage: python poresippr_basecall_scheduler.py <input.csv> <metadata.csv>")

def signal_handler(signum, frame):
    global terminate
    print("Signal received, terminating the script...")
    terminate = True

# Set up argument parsing
parser = argparse.ArgumentParser(description='This script runs Guppy basecaller and downstream analysis such as mapping with minimap2 and samtools.')
parser.add_argument('csv_file_path', help='Path to the input CSV file.')
parser.add_argument('metadata_csv_path', help='Path to the metadata CSV file.')
args = parser.parse_args()

# Reading the inputs from CSV file
if len(sys.argv) != 3:
    print_usage()
    sys.exit(1)

csv_file_path = args.csv_file_path
metadata_csv_path = args.metadata_csv_path

# Read Metadata CSV and create a mapping of barcode to seqID and olnid
barcode_to_seqid = {}
barcode_to_olnid = {}
with open(metadata_csv_path, 'r') as metafile:
    meta_reader = csv.DictReader(metafile)
    header = meta_reader.fieldnames
    print(f"Metadata CSV columns: {header}")  # Debug: Print column names
    for row in meta_reader:
        print(f"Row: {row}")  # Debug: Print each row
        barcode_to_seqid[int(row['Barcode'])] = row['SEQID']
        barcode_to_olnid[int(row['Barcode'])] = row['OLNID']

def main_loop(complete, iteration):
    while not terminate:  # Loop to keep running Guppy and downstream analysis
        with open(csv_file_path, 'r') as csvfile:
            csv_reader = csv.DictReader(csvfile)
            for row in csv_reader:
                reference = row['reference']
                fast5_dir = row['fast5_dir']
                output_dir = row['output_dir']
                config = row['config']
                barcode = row['barcode']
                barcode_values = [int(x.strip().replace('"', '')) for x in row['barcode_values'].split(',')]

                # Creating the output directory if it doesn't exist
                try:
                    os.makedirs(output_dir, exist_ok=True)
                except OSError as e:
                    print(f"Error creating output directory: {output_dir}")
                    print(f"Error message: {e}")

                # Running guppy fast basecalling
                guppy_command = f"guppy_basecaller --input_path {fast5_dir} --save_path {output_dir} --config {config} --barcode_kits {barcode}  -x auto -r"
                print(f"Running guppy command: {guppy_command}")
                run_command(guppy_command)

                # Processing each provided barcode number by user as an argument in CSV file, enclosed in double quotes
                for barcode_num in barcode_values:
                    barcode_dir = f"barcode{barcode_num:02d}"  # Forming and looking for the barcode directory name
                    barcode_path = os.path.join(output_dir, "pass", barcode_dir)
                    if os.path.isdir(barcode_path):
                        seqid = barcode_to_seqid.get(barcode_num, f"barcode{barcode_num:02d}")
                        olnid = barcode_to_olnid.get(barcode_num, f"OLN{barcode_num:02d}")
                        concatenated_fastq_file = os.path.join(output_dir, f"{seqid}.fastq")
                        concatenate_fastq_files(barcode_path, concatenated_fastq_file)

                        # Aligning using Minimap2 and create sorted BAM file
                        bam_file = os.path.join(output_dir, f"{seqid}_sorted.bam")
                        minimap2_command = f"minimap2 -ax map-ont {reference} {concatenated_fastq_file} | samtools view -@ 5 -bS - | samtools sort -o {bam_file} -"
                        print(f"Running Minimap2 command: {minimap2_command}")
                        run_command(minimap2_command)

                        # Calculating coverage and sorting the output from samtools out
                        csv_file = os.path.join(output_dir, f"{seqid}_iteration{iteration.value}.csv")
                        samtools_command = f"samtools coverage {bam_file} | cut -f 1,4 | awk '$2 > 0' | sort -rnk 2,2 | sed 's/\\t/,/g' > {csv_file}"
                        print(f"Running Samtools command: {samtools_command}")
                        run_command(samtools_command)

                        # Adding header to the output CSV file
                        with open(csv_file, 'r+') as f:
                            content = f.read()
                            f.seek(0, 0)
                            f.write("gene_name,number_of_reads_mapped\n" + content)

                        file_size = os.path.getsize(concatenated_fastq_file)
                        if file_size > 0:
                            genome_coverage_value = file_size / 5000000
                            output_string = ["genome_coverage", f"{genome_coverage_value}X"]
                            with open(csv_file, mode='a', newline='') as file:
                                writer = csv.writer(file)
                                writer.writerow(output_string)
                        else:
                            print(f"Error: Unable to determine file size for {concatenated_fastq_file}")

                        # Deleting the temporary concatenated FASTQ files else it will throw a memory error
                        os.remove(concatenated_fastq_file)
                        os.remove(bam_file)
                    else:
                        print(f"No barcode directory found for {barcode_dir}")

        iteration.value += 1
        if terminate:
            break
        print("Waiting for 30 minutes before running Guppy again...")
        time.sleep(1800)
        complete.value = True

if __name__ == "__main__":
    terminate = False
    complete = Value('b', False)
    iteration = Value('i', 1)
    p = Process(target=main_loop, args=(complete, iteration))

    # Register the signal handler for SIGINT (Ctrl+C) and SIGTERM
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    # Start the subprocess
    p.start()

    # Wait for the subprocess to finish
    while not complete.value:
        try:
            p.join(timeout=1)
        except SystemExit:
            print('Terminating subprocess...')
            p.terminate()
            p.join()
        if complete.value:
            break

    if terminate:
        print('Terminating subprocess...')
        p.terminate()
        p.join()
