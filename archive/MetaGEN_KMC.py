# -------------------------------
# Title: MetaGEN_KMC.py
# Author: Silver A. Wolf
# Last Modified: Fri, 23.07.2021
# Version: 0.0.2
# -------------------------------

# Imports
from fnmatch import fnmatch, fnmatchcase
from glob import glob
import argparse
import os

# Settings
#name_project = "MetaSUB/"
name_project = "Horses_Gut/"
output_folder = "output/" + name_project
version = "0.0.2"

def run_coverm(kmc_input, kmc_output, threads):
	print("Step 12/12 - Calculating k-mer statistics [KMC]:\n")
	
	read_list = sorted([name for name in os.listdir(kmc_input) if fnmatch(name, "*_R1.fastq.gz")])
	os.system("mkdir -p " + kmc_output)
	c = 0
	
	# PE + SE
	for read in read_list:
		sample_name = read.split("_R1")[0]
		read_1 = kmc_input + read
		read_2 = kmc_input + sample_name + "_R2.fastq.gz"
		read_3 = kmc_input + sample_name + "_R3.fastq.gz"
		
		f = open("tmp/sample_list.txt", "w")
		f.write(read_1 + "\n" + read_2 + "\n" + read_3)
		f.close()
		
		print("KMC: Assessing k-mers for " + sample_name + ".")
		os.system("kmc " +
				  "@tmp/sample_list.txt " +
				  kmc_output + sample_name + " "
				  "tmp/ " +
				  "-m100 " +
				  "-sm " +
				  "-fq " +
				  "-ci0 " +
				  "-cs999 " +
				  "-t" + threads
				 )
		os.system("rm tmp/sample_list.txt")
		c = c + 1

	print("KMC: " + str(c) + " files successfully analyzed.")
	print("KMC: Finished.\n")
	
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
	
	print("Running MetaGEN-KMC Pipeline Version " + version + "\n")
	
	run_coverm(output_folder + "fastp/",
			   output_folder + "kmc/",
			   str(args.threads)
			  )
	
	print("MetaGEN-KMC: Finished.")

if __name__ == "__main__":
	main()