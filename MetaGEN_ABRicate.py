# -------------------------------
# Title: MetaGEN_ABRicate.py
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
name_project = "Horses/"
output_folder = "output/" + name_project
version = "0.0.4"

def run_abricate(metaspades_input, abricate_output, threads):
	print("Step 10/10 - Scanning for AMR Genes [ABRicate]:\n")
	os.system("mkdir -p " + abricate_output)
	
	file_list = glob(metaspades_input + "*.gz")

	c = 0
	for assembly in file_list:
		sample_name = assembly.split("/")[-1].split(".fa")[0]
		print("ABRicate: Analyzing " + sample_name + ".")
		os.system("abricate" +
				  " --db megares" +
				  " --threads " + threads +
				  " --nopath" + 
				  " " + assembly +
				  " > " + abricate_output + sample_name + ".tab"
				 )
		c = c + 1
		print("")
	
	print("ABRicate: Summarizing results.\n")
	os.system("abricate --summary --nopath " + abricate_output + "*.tab > " + abricate_output + "summary.tab")
	
	print("ABRicate: " + str(c) + " files successfully analyzed.")
	print("ABRicate: Finished.\n")
	
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
	
	print("Running MetaGEN-ABRicate Pipeline Version " + version + "\n")
	
	run_abricate(output_folder + "metaspades/",
				 output_folder + "abricate/",
				 str(args.threads)
				)
	
	print("MetaGEN-ABRicate: Finished.")

if __name__ == "__main__":
	main()