# -------------------------------
# Title: MetaGEN_Main.py
# Author: Silver A. Wolf
# Last Modified: Thu, 10.06.2021
# Version: 0.3.1
# -------------------------------

# Imports
from fnmatch import fnmatch, fnmatchcase
from glob import glob
import argparse
import csv
import os

# Settings

# MetaSUB
#name_project = "MetaSUB/"
#read_1_identifier = "read1.fastq.gz"
#read_2_identifier = "read2.fastq.gz"

# Horses
name_project = "Horses/"
read_1_identifier = "R1_001.fastq.gz"
read_2_identifier = "R2_001.fastq.gz"

# General
input_folder = "input/" + name_project
input_metadata = input_folder + "metadata/"
input_sequences = input_folder + "sequences/"
output_folder = "output/" + name_project
version = "0.3.1"

def download_metasub(city):
	print("Step 1/8 - Fetching Data [MetaSUB]:\n")
	os.system("mkdir -p " + input_metadata)
	os.system("mkdir -p " + input_sequences)
	
	print("MetaSUB: Downloading Metadata Table.")
	os.system("wget -N -P " + input_metadata + " " +
			  "https://github.com/MetaSUB/MetaSUB-metadata/raw/master/complete_metadata.csv")
	
	with open(input_metadata + "complete_metadata.csv") as metadata_file:
		for line in csv.reader(metadata_file, delimiter = ","):
			if len(line) > 0:
				if city in line[4]:
					sample_name = line[0].strip()
					sample_list = sample_name.split("_")
					print("MetaSUB: Downloading " + sample_name + ".")
					command = ("wget -N -P " + input_sequences + " " +
							   "http://s3.wasabisys.com/metasub/human_filtered_data/hudson_alpha_library/" +
							   sample_list[0] + "/" + sample_list[1] + "/" + sample_name +
							   ".filter_human_dna.nonhuman_read"
							  )
					os.system(command + "1.fastq.gz")
					os.system(command + "2.fastq.gz")

	# Move known SE sequences
	print("Unix: Cleaning up output.")
	os.system("mv " + input_sequences + "haib17CEM4890_H7KYMCCXY_SL273100* " + input_metadata)
	
	print("MetaSUB: Finished.\n")	
	
def run_fastp(memory, threads):
	print("Step 2/8 - Quality Control [fastp]:\n")
	os.system("mkdir -p " + output_folder + "fastp/reports/")
	
	read_list_pe = sorted([name for name in os.listdir(input_sequences) if fnmatch(name, "*" + read_1_identifier)])
	fastp_dir = output_folder + "fastp/"
	c = 0
	
	for read in read_list_pe:
		# Define output naming convention
		if name_project == "MetaSUB/":
			sample_name = read.split(".filter")[0].split("_")[-1]
		else:
			sample_name = read.split("-L")[0]

		print("fastp: Analyzing " + sample_name + " (PE).")

		bbsplit_R12_scaf = output_folder + "bbmap/" + sample_name + "_scaffolds_12.txt"
		bbsplit_R12_stat = output_folder + "bbmap/" + sample_name + "_stats_12.txt"
		bbsplit_R3_scaf = output_folder + "bbmap/" + sample_name + "_scaffolds_3.txt"
		bbsplit_R3_stat = output_folder + "bbmap/" + sample_name + "_stats_3.txt"
		
		read1_in = input_sequences + read
		read2_in = input_sequences + read.split(read_1_identifier)[0] + read_2_identifier
		
		read1_out = fastp_dir + sample_name + "_R1.fastq.gz"
		read2_out = fastp_dir + sample_name + "_R2.fastq.gz"
		read3_out = fastp_dir + sample_name + "_R3.fastq.gz"
		
		read1_tmp = "tmp/" + sample_name + "_R1.tmp.fastq.gz"
		read2_tmp = "tmp/" + sample_name + "_R2.tmp.fastq.gz"
		read3_tmp = "tmp/" + sample_name + "_R3.tmp.fastq.gz"
		
		html_out = fastp_dir + "reports/" + sample_name + ".fastp.html"
		json_out = fastp_dir + "reports/" + sample_name + ".fastp.json"
		
		print("fastp: Performing QC.")
		
		os.system("fastp" +
				  " --in1 " + read1_in +
				  " --out1 " + read1_tmp +
				  " --in2 " + read2_in +
				  " --out2 " + read2_tmp +
				  " --unpaired1 " + read3_tmp +
				  " --unpaired2 " + read3_tmp +
				  " --html " + html_out +
				  " --json " + json_out +
				  " --thread " + threads +
				  " --compression 9" +
				  " --correction" +
				  " --cut_tail" +
				  " --cut_window_size 4" +
				  " --cut_mean_quality 20" +
				  " --detect_adapter_for_pe" +
				  " --overrepresentation_analysis" +
				  " --overrepresentation_sampling 20"
				 )
		
		#print("BBMap: Generating Host Index.")
		#os.system("rm -r references/ref/")
		#os.system("bbsplit.sh" +
		#		  " ref=references/" +
		#		  " path=references/" +
		#		  " threads=" + threads +
		#		  " -Xmx" + memory + "g"
		#		 )
		
		print("BBMap: Removing Host Genome Contamination.")
		
		os.system("bbsplit.sh" +
				  " in1=" + read1_tmp +
				  " in2=" + read2_tmp +
				  " outu1=" + read1_out +
				  " outu2=" + read2_out +
				  " path=references/" +
				  " threads=" + threads +
				  " refstats=" + bbsplit_R12_stat +
				  " scafstats=" + bbsplit_R12_scaf +
				  " -Xmx" + memory + "g"
				 )
		
		os.system("bbsplit.sh" +
				  " in=" + read3_tmp +
				  " outu=" + read3_out +
				  " path=references/" +
				  " threads=" + threads +
				  " refstats=" + bbsplit_R3_stat +
				  " scafstats=" + bbsplit_R3_scaf +
				  " -Xmx" + memory + "g"
				 )
		
		print("Unix: Cleaning up output.")
		os.system("rm -r tmp/*")
			
		c = c + 1
		print("")
		
	print("fastp: " + str(c) + " files successfully analyzed.")
	print("fastp: Finished.\n")

def run_kraken2_assembly(database, confidence_score, read_length, bracken_threshold, spades_dir, bracken_dir, kraken_dir, threads):
	print("Step 3/8 - Taxonomic Classification [kraken2]:\n")
	os.system("mkdir -p " + bracken_dir)
	os.system("mkdir -p " + kraken_dir)
	
	print("kraken2: Computing Taxonomic Classification.")
	file_list = glob(spades_dir + "*.gz")
	levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
	c = 0
	
	for assembly in file_list:
		sample_name = assembly.split("/")[-1].split(".gz")[0]
		print("kraken2: Analyzing " + sample_name + ".")
		os.system("kraken2" +
				  " --db " + database +
				  " --threads " + threads +
				  " --confidence " + confidence_score +
				  " --report " + kraken_dir + sample_name + ".report" +
				  " " + assembly +
				  " > " + kraken_dir + sample_name + ".stdout"
				 )
		for level in levels:
			print("Bracken: Analyzing " + sample_name + " (" + level + ").")
			os.system("bracken" +
					  " -d " + database +
					  " -i " + kraken_dir + sample_name + ".report" +
					  " -o " + bracken_dir + sample_name + "." + level.lower() + ".bracken" +
					  " -w " + bracken_dir + sample_name + "." + level.lower() + ".report" +
					  " -r " + read_length +
					  " -l " + level[0] +
					  " -t " + bracken_threshold
					 )
		c = c + 1
		print("")

	for level in levels:
		print("Bracken: Merging " + str(c) + " results (" + level + ").")
		os.system("combine_bracken_outputs.py" +
				  " --files " + bracken_dir + "*." + level.lower() + ".bracken" +
				  " -o " + bracken_dir + "bracken_" + level.lower() + ".bracken"
				 )
		
		print("kraken-biom: Exporting bracken reports (" + level + ").")
		os.system("kraken-biom " + bracken_dir + "*." + level.lower() + ".report" +
				  " -o " + bracken_dir + "bracken_" + level.lower() + ".biom" +
				  " --fmt " + "json" +
				  " --max " + "D" +
				  " --min " + "S"
				 )
	
		print("kraken-biom: Exporting kraken2 reports.")
		os.system("kraken-biom " + kraken_dir + "*.report" +
				  " -o " + kraken_dir + "kraken2.biom" +
				  " --fmt " + "json" +
				  " --max " + "D" +
				  " --min " + "S"
				 )
	
	print("kraken2: " + str(c) + " files successfully analyzed.")
	print("kraken2: Finished.\n")
	
def run_kraken2_reads(database, confidence_score, read_length, bracken_threshold, fastp_dir, bracken_dir, kraken_dir, threads):
	print("Step 3/8 - Taxonomic Classification [kraken2]:\n")
	os.system("mkdir -p " + bracken_dir)
	os.system("mkdir -p " + kraken_dir)
	
	print("kraken2: Computing Taxonomic Classification.")
	read_list = sorted([name for name in os.listdir(fastp_dir) if fnmatch(name, "*_R1.fastq.gz")])
	levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
	c = 0

	for read1 in read_list:
		sample_name = read1.split("_R1")[0]
		read2 = sample_name + "_R2.fastq.gz"
		print("kraken2: Analyzing " + sample_name + ".")
		os.system("kraken2" +
				  " --db " + database +
				  " --threads " + threads +
				  " --confidence " + confidence_score +
				  " --report " + kraken_dir + sample_name + ".report" +
				  " --paired " + fastp_dir + read1 + " " + fastp_dir + read2 +
				  " > " + kraken_dir + sample_name + ".stdout"
				 )
		for level in levels:
			print("Bracken: Analyzing " + sample_name + " (" + level + ").")
			os.system("bracken" +
					  " -d " + database +
					  " -i " + kraken_dir + sample_name + ".report" +
					  " -o " + bracken_dir + sample_name + "." + level.lower() + ".bracken" +
					  " -w " + bracken_dir + sample_name + "." + level.lower() + ".report" +
					  " -r " + read_length +
					  " -l " + level[0] +
					  " -t " + bracken_threshold
					 )
		c = c + 1
		print("")

	for level in levels:
		print("Bracken: Merging " + str(c) + " results (" + level + ").")
		os.system("combine_bracken_outputs.py" +
				  " --files " + bracken_dir + "*." + level.lower() + ".bracken" +
				  " -o " + bracken_dir + "bracken_" + level.lower() + ".bracken"
				 )
		
		print("kraken-biom: Exporting bracken reports (" + level + ").")
		os.system("kraken-biom " + bracken_dir + "*." + level.lower() + ".report" +
				  " -o " + bracken_dir + "bracken_" + level.lower() + ".biom" +
				  " --fmt " + "json" +
				  " --max " + "D" +
				  " --min " + "S"
				 )
	
		print("kraken-biom: Exporting kraken2 reports.")
		os.system("kraken-biom " + kraken_dir + "*.report" +
				  " -o " + kraken_dir + "kraken2.biom" +
				  " --fmt " + "json" +
				  " --max " + "D" +
				  " --min " + "S"
				 )
	
	print("kraken2: " + str(c) + " files successfully analyzed.")
	print("kraken2: Finished.\n")

def run_krona(krona_input, krona_output):
	print("Step 4/8 - Visualizing Taxonomic Classification [Krona]:\n")
	os.system("mkdir -p " + krona_output)

	#print("Krona: Updating Taxonomic Database.")
	#os.system("ktUpdateTaxonomy.sh")
	
	print("Krona: Preprocessing Taxonomic Files.")
	kraken2_list_in = sorted([name for name in os.listdir(krona_input) if fnmatch(name, "*.stdout")])
	kraken2_list_ext = []
	
	for sample in kraken2_list_in:
		new_sample = krona_input + sample + "," + sample.split(".stdout")[0]
		kraken2_list_ext.append(new_sample)
	
	kraken2_list_str = ' '.join(kraken2_list_ext)
	
	print("Krona: Visualizing Taxonomic Classification.")
	os.system("ktImportTaxonomy -q 2 -t 3 " + kraken2_list_str + " -o " + krona_output + "krona.html")	
	
	print("Krona: Finished.\n")
	
def run_multiqc(fastp_input, fastp_output, kraken2_input, kraken2_output):
	print("Step 5/8 - Quality Control [multiqc]:\n")
	
	print("MultiQC: Summarizing Results (fastp).")
	os.system("mkdir -p " + fastp_output)
	os.system("multiqc -o " + fastp_output + " " + fastp_input)
	
	print("MultiQC: Summarizing Results (kraken2).")
	os.system("mkdir -p " + kraken2_output)
	os.system("multiqc -o " + kraken2_output + " " + kraken2_input)
	
	print("multiqc: Finished.\n")	

def run_metaspades(metaspades_input, metaspades_output, threads):
	print("Step 6/8 - Metagenome Assembly [metaSPAdes]:\n")
	
	print("metaSPAdes: Creating Assemblies.")
	read_list = sorted([name for name in os.listdir(metaspades_input) if fnmatch(name, "*_R1.fastq.gz")])
	c = 0
	
	os.system("mkdir -p " + metaspades_output)
	
	# PE + SE
	for read in read_list:
		sample_name = read.split("_R1")[0]
		read_1 = metaspades_input + read
		read_2 = metaspades_input + sample_name + "_R2.fastq.gz"
		read_3 = metaspades_input + sample_name + "_R3.fastq.gz"
		print("metaSPAdes: Assembling " + sample_name + ".")
		
		os.system("spades.py" +
				  " -o " + metaspades_output + sample_name +
				  " --meta " +
				  " -1 " + read_1 +
				  " -2 " + read_2 +
				  " -s " + read_3 +
				  " -t " + threads +
				  " --tmp-dir tmp/"
				 )
		
		print("Unix: Cleaning up output.")
		os.system("mv " + metaspades_output + sample_name + "/scaffolds.fasta " + metaspades_output + sample_name + ".fa")
		os.system("rm -r " + metaspades_output + sample_name + "/")
		c = c + 1
		print("")
	
	print("Unix: Compressing Assemblies.")
	os.system("gzip " + metaspades_output + "*.fa")
	print("")
	
	print("metaSPAdes: " + str(c) + " files successfully assembled.")
	print("metaSPAdes: Finished.\n")

def run_metaquast(metaquast_input, metaquast_output, threads):
	print("Step 7/8 - Quality Control [MetaQUAST]:\n")
	
	file_list = glob(metaquast_input + "*.gz")
	c = 0
	
	for assembly in file_list:
		sample_name = assembly.split("/")[-1].split(".fa")[0]
		print("MetaQUAST: Analyzing " + sample_name + ".")
		os.system("metaquast.py " + assembly +
				  " -o " + metaquast_output + sample_name +
				  " -t " + threads +
				  " --gene-finding"
				 )
		c = c + 1
		print("")

	print("MetaQUAST: " + str(c) + " files successfully analyzed.")
	print("MetaQUAST: Finished.\n")

def run_metabat(fastp_input, metaspades_input, metabat_output, threads):
	print("Step 8/8 - Assembly Binning [MetaBAT]:\n")
	os.system("mkdir -p " + metabat_output)
	
	file_list = glob(metaspades_input + "*.gz")
	c = 0
	
	for assembly in file_list:
		sample_name = assembly.split("/")[-1].split(".fa")[0]
		print("MetaBAT: Analyzing " + sample_name + ".")
		
		bowtie_index = "tmp/" + sample_name
		bowtie_bam = bowtie_index + ".bam"
		bowtie_sam = bowtie_index + ".sam"
		bowtie_sorted = bowtie_index + "_sorted.bam"
		
		print("bowtie2: Generating Assembly Index.")
		os.system("bowtie2-build" +
				  " --quiet" +
				  " --threads " + threads +
				  " " + assembly +
				  " " + bowtie_index
				 )
		
		print("bowtie2: Mapping Reads To Assembly.")
		os.system("bowtie2" +
				  " --quiet" +
				  " -p " + threads +
				  " -x " + bowtie_index +
				  " -1 " + fastp_input + sample_name + "_R1.fastq.gz" +
				  " -2 " + fastp_input + sample_name + "_R2.fastq.gz" +
				  " -U " + fastp_input + sample_name + "_R3.fastq.gz" +
				  " -S " + bowtie_sam
				 )
		
		print("samtools: Convert SAM To BAM.")
		os.system("samtools view -bS" +
				  " -o " + bowtie_bam +
				  " " + bowtie_sam
				 )

		print("samtools: Sort BAM.")
		os.system("samtools sort " + bowtie_bam + " -o " + bowtie_sorted)
		
		print("samtools: Index BAM.")
		os.system("samtools index " + bowtie_sorted)
		
		print("MetaBAT: Binning Contigs.")
		os.system("runMetaBat.sh"
				  " -m 1500 " +
				  " -t " + threads +
				  " " + assembly +
				  " " + bowtie_sorted
				 )
		
		print("Unix: Cleaning up output.")
		os.system("mv " + sample_name + "* " + metabat_output)
		os.system("rm -r tmp/*")
		
		c = c + 1
		print("")
	
	print("MetaBAT: " + str(c) + " files successfully analyzed.")
	print("MetaBAT: Finished.\n")

# Main
def main():
	parser = argparse.ArgumentParser(description = "")
	parser.add_argument("-b",
						"--bracken_threshold",
						type = int,
						default = "10",
						required = False,
						help = "Bracken read threshold"
					   )
	parser.add_argument("-c",
						"--city_of_interest",
						type = str,
						default = "berlin",
						required = False,
						help = "MetaSUB: Name of the city of interest"
					   )
	parser.add_argument("-k",
						"--kraken2_db",
						type = str,
						default = "/scratch1/databases/kraken/20210429_kraken2_standard_new_database/",
						required = False,
						help = "Path to local kraken2 database"
					   )
	parser.add_argument("-m",
						"--memory",
						type = int,
						default = "32",
						required = False,
						help = "Amount of memory used for running MetaGEN"
					   )
	parser.add_argument("-r",
						"--read_length",
						type = int,
						default = "150",
						required = False,
						help = "Minimum expected read length"
					   )
	parser.add_argument("-s",
						"--kraken2_confidence_score",
						type = float,
						default = "0.0",
						required = False,
						help = "Kraken2 confidence score"
					   )
	parser.add_argument("-t",
						"--threads",
						type = int,
						default = "64",
						required = False,
						help = "Amount of threads used for running MetaGEN"
					   )

	args = parser.parse_args()
	
	print("Running MetaGEN Pipeline Version " + version + "\n")
	
	#download_metasub(args.city_of_interest)
	#run_fastp(str(args.memory),
	#		  str(args.threads)
	#		 )
	run_kraken2_reads(args.kraken2_db,
					  str(args.kraken2_confidence_score),
					  str(args.read_length),
					  str(args.bracken_threshold),
					  output_folder + "fastp/",
					  output_folder + "bracken/reads/",
					  output_folder + "kraken2/reads/",
					  str(args.threads)
					 )
	run_krona(output_folder + "kraken2/reads/",
			  output_folder + "krona/reads/"
			 )
	run_multiqc(output_folder + "fastp/reports/",
				output_folder + "multiqc/fastp/",
				output_folder + "kraken2/reads/",
				output_folder + "multiqc/kraken2/"
			   )
	#run_metaspades(output_folder + "fastp/",
	#			   output_folder + "metaspades/",
	#			   str(args.threads)
	#			  )
	run_kraken2_assembly(args.kraken2_db,
						 str(args.kraken2_confidence_score),
						 str(args.read_length),
						 str(args.bracken_threshold),
						 output_folder + "kraken2/reads/",
						 output_folder + "bracken/assembly/",
						 output_folder + "kraken2/assembly/",
						 str(args.threads)
						)
	#run_metaquast(output_folder + "metaspades/",
	#			  output_folder + "metaquast/",
	#			  str(args.threads)
	#			 )
	#run_metabat(output_folder + "fastp/",
	#			output_folder + "metaspades/",
	#			output_folder + "metabat/",
	#			str(args.threads)
	#		   )
	
	print("MetaGEN: Finished.")

if __name__ == "__main__":
	main()