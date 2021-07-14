# -------------------------------
# Title: MetaGEN_CoverM.py
# Author: Silver A. Wolf
# Last Modified: Fri, 11.06.2021
# Version: 0.0.4
# -------------------------------

# Imports
from glob import glob
import argparse
import os

# Settings
#name_project = "MetaSUB/"
name_project = "Horses_Gut/"
output_folder = "output/" + name_project
version = "0.0.1"

def run_coverm(coverm_input, coverm_output, threads):
	print("Step 11/11 - AMR Scan [CoverM]:\n")
	
	read_list = sorted([name for name in os.listdir(coverm_input) if fnmatch(name, "*_R1.fastq.gz")])
	c = 0
	
	os.system("mkdir -p " + coverm_output)
	
	# PE + SE
	for read in read_list:
		sample_name = read.split("_R1")[0]
		read_1 = coverm_input + read
		read_2 = coverm_input + sample_name + "_R2.fastq.gz"
		read_3 = coverm_input + sample_name + "_R3.fastq.gz"
		print("CoverM: Mapping " + sample_name + " to MegaRes database.")
		os.system("coverm contig " +
				  "-1 " + read_1 + " " +
				  "-2 " + read_2 + " " +
				  "--single " + read_3 + " " +
				  "-r megares_full_database_v2.00.fasta " + 
				  "-p bwa-mem " + 
				  "-m tpm " +
				  "-o " + coverm_output + sample_name + ".txt "
				  "-t " + threads
				 )

	print("CoverM: " + str(c) + " files successfully analyzed.")
	print("CoverM: Finished.\n")
	
# Main
def main():
	parser = argparse.ArgumentParser(description = "")
	parser.add_argument("-t",
						"--threads",
						type = int,
						default = "64",
						required = False,
						help = "Amount of threads used for running MetaGEN"
					   )
	args = parser.parse_args()
	
	print("Running MetaGEN-CoverM Pipeline Version " + version + "\n")
	
	run_coverm(output_folder + "fastp/",
				 output_folder + "coverm/",
				 str(args.threads)
				)
	
	print("MetaGEN-CoverM: Finished.")

if __name__ == "__main__":
	main()