# -------------------------------
# Title: MetaGEN_CoverM.py
# Author: Silver A. Wolf
# Last Modified: Fri, 16.07.2021
# Version: 0.0.3
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
version = "0.0.3"

def run_coverm(coverm_input, coverm_output, threads):
	print("Step 11/11 - AMR Scan [CoverM]:\n")
	
	#print("USEARCH: Cluster MegaRES database (90% identity).")
	#os.system("./usearch11.0.667_i86linux32 " +
	#		  "-cluster_fast amr/megares_full_database_v2.00.fasta " +
	#		  "-id 0.90 " +
	#		  "-centroids amr/megares.fasta " +
	#		  "-uc amr/megares.uc"
	#		 )
	
	read_list = sorted([name for name in os.listdir(coverm_input) if fnmatch(name, "*_R1.fastq.gz")])
	os.system("mkdir -p " + coverm_output)
	c = 0
	
	# PE + SE
	for read in read_list:
		sample_name = read.split("_R1")[0]
		read_1 = coverm_input + read
		read_2 = coverm_input + sample_name + "_R2.fastq.gz"
		read_3 = coverm_input + sample_name + "_R3.fastq.gz"
		print("CoverM: Mapping " + sample_name + " to MegaRES database.")
		os.system("coverm contig " +
				  "-1 " + read_1 + " " +
				  "-2 " + read_2 + " " +
				  "--single " + read_3 + " " +
				  "-r amr/megares.fasta " + 
				  "-p minimap2-sr " + 
				  "-m tpm " +
				  "-o " + coverm_output + sample_name + ".txt "
				  "-t " + threads
				 )
		c = c + 1
	print("CoverM: " + str(c) + " files successfully analyzed.")

	print("CoverM: Merging results.")
	coverm_results = sorted([name for name in os.listdir(coverm_output) if fnmatch(name, "*.txt")])
	ref_names = []
	ref_seqs = "MegaRes"
	
	for c in coverm_results:
		ref_seqs = ref_seqs + "\t" + c.split(".txt")[0]
	
	with open(coverm_output + coverm_results[0], "r") as ref_file:
		for line in ref_file:
			ref_names.append(line.split("\t")[0])
	
	ref_names.remove("Contig")
	f = open(coverm_output + "summary.tsv", "w")
	f.write(ref_seqs + "\n")
	
	for r in ref_names:
		val = []
		for l in coverm_results:
			with open(coverm_output + l, "r") as cover_file:
				for line in cover_file:
					if line.split("\t")[0] == r:
						val.append(line.split("\t")[1])
						break
		f.write(r + "\t" + "\t".join(val) + "\n")
	
	f.close()
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