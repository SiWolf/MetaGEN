# -------------------------------
# Title: MetaGEN_Run_Pipeline.py
# Author: Silver A. Wolf
# Last Modified: Mon, 25.01.2021
# Version: 0.2.4
# -------------------------------

# Imports
from fnmatch import fnmatch, fnmatchcase
from glob import glob
import argparse
import csv
import os

# Settings
name_project = "MetaSUB/"
input_folder = "input/" + name_project
input_pe = input_folder + "PE/"
input_se = input_folder + "SE/"
read_1_identifier = "read1.fastq.gz"
read_2_identifier = "read2.fastq.gz"
version = "0.2.4"

def set_reference(host):
	ref_seq = ""
	if host == "human":
		ref_seq = "references/GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
	elif host == "horse":
		ref_seq = "references/GCF_002863925.1_EquCab3.0_genomic.fna.gz"
	else:
		print("Error - Unknown reference species, please try a different species.")
		quit()
	return(ref_seq)

def download_metasub(city):
	print("Step 1/8 - Fetching Data [MetaSUB]:\n")
	os.system("mkdir -p " + input_pe)
	os.system("mkdir -p " + input_se)
	
	print("MetaSUB: Downloading Metadata Table.")
	os.system("wget -N -P " + input_folder + " https://github.com/MetaSUB/MetaSUB-metadata/raw/master/complete_metadata.csv")
	
	with open(input_folder + "complete_metadata.csv") as metadata_file:
		for line in csv.reader(metadata_file, delimiter = ","):
			if len(line) > 0:
				if city in line[4]:
					sample_name = line[0].strip()
					sample_list = sample_name.split("_")
					print("MetaSUB: Downloading " + sample_name + ".")
					command = "wget -N -P " + input_pe + " http://s3.wasabisys.com/metasub/human_filtered_data/hudson_alpha_library/" + sample_list[0] + "/" + sample_list[1] + "/" + sample_name + ".filter_human_dna.nonhuman_read"
					os.system(command + "1.fastq.gz")
					os.system(command + "2.fastq.gz")		

	# Move known SE sequences
	print("Unix: Cleaning up output.")
	os.system("mv " + input_pe + "haib17CEM4890_H7KYMCCXY_SL273100* " + input_se)
	
	print("MetaSUB: Finished.\n")	
	
def run_fastp(host, memory, threads):
	print("Step 2/8 - Quality Control [fastp]:\n")
	os.system("mkdir -p output/fastp/reports/")
	
	read_list_pe = sorted([name for name in os.listdir(input_pe) if fnmatch(name, "*" + read_1_identifier)])
	read_list_se = sorted([name for name in os.listdir(input_se) if fnmatch(name, "*" + read_2_identifier)])
	output_dir = "output/fastp/"
	c = 0
	
	for read in read_list_pe:
		# Output naming convention
		sample_name = read.split(".filter")[0].split("_")[-1]

		print("fastp: Analyzing " + read + " (PE).")

		bbsplit_scaf_12 = "output/bbmap/" + sample_name + "_scaffolds_12.txt"
		bbsplit_stats_12 = "output/bbmap/" + sample_name + "_stats_12.txt"
		bbsplit_scaf_3 = "output/bbmap/" + sample_name + "_scaffolds_3.txt"
		bbsplit_stats_3 = "output/bbmap/" + sample_name + "_stats_3.txt"
		
		read_1_in = input_pe + read
		read_2_in = input_pe + read.split(read_1_identifier)[0] + read_2_identifier
		read_1_out = output_dir + sample_name + "_" + read_1_identifier
		read_2_out = output_dir + sample_name + "_" + read_2_identifier
		read_3_out = output_dir + sample_name + "_read3.fastq.gz"
		read_4_out = output_dir + sample_name + "_read4.fastq.gz"
		read_html_out = output_dir + "reports/" + sample_name + ".fastp.html"
		read_json_out = output_dir + "reports/" + sample_name + ".fastp.json"
		read_merged_se_out = "tmp/" + sample_name + ".reads.se.merged.fastq.gz"
		read_unpaired_out = "tmp/" + sample_name + ".reads.unpaired.fastq.gz"
		read_tmp_1 = "tmp/" + sample_name + "_R1.tmp.fastq.gz"
		read_tmp_2 = "tmp/" + sample_name + "_R2.tmp.fastq.gz"
		
		print("fastp: Performing QC.")
		os.system("fastp --in1 " + read_1_in + " --out1 " + read_tmp_1 + " --in2 " + read_2_in + " --out2 " + read_tmp_2 + " --unpaired1 " + read_unpaired_out + " --unpaired2 " + read_unpaired_out + " --merge --merged_out " + read_merged_se_out + " --compression 9 --detect_adapter_for_pe --thread " + threads + " --html " + read_html_out + " --json " + read_json_out + " --overrepresentation_analysis --overrepresentation_sampling 20 --correction --cut_tail --cut_window_size 4 --cut_mean_quality 20")
		
		print("BBMap: Removing Host Genome Contamination.")
		os.system("bbsplit.sh in1=" + read_tmp_1 + " in2=" + read_tmp_2 + " ref=" + host + " outu1=" + read_1_out + " outu2=" + read_2_out + " threads=" + threads + " path=tmp/ -Xmx" + memory + "g scafstats=" + bbsplit_scaf_12 + " refstats=" + bbsplit_stats_12)
		
		print("Unix: Concatenating Merged and Unpaired SE Reads.")
		read_merged_final_out = "tmp/" + sample_name + ".reads.final.merged.fastq.gz" 
		os.system("cat " + read_merged_se_out + " " + read_unpaired_out + " > " + read_merged_final_out)
		
		print("BBMap: Removing Host Genome Contamination.")
		os.system("bbsplit.sh in=" + read_merged_final_out + " ref=" + host + " outu=" + read_3_out + " threads=" + threads + " path=tmp/ -Xmx" + memory + "g scafstats=" + bbsplit_scaf_3 + " refstats=" + bbsplit_stats_3)
		
		print("BBMap: Reverse Complementing Read 2.")  
		read_2_rev_out = "tmp/" + sample_name + "_read2.reverse.fastq.gz"
		os.system("reformat.sh in=" + read_2_out + " out=" + read_2_rev_out + " rcomp=t -Xmx" + memory + "g zl=9 threads=" + threads)
		
		print("BBMap: Concatenating PE Reads.")
		read_merged_pe_out = "tmp/" + sample_name + "_read1_2.merged.fastq.gz"
		os.system("fuse.sh in1=" + read_1_out + " in2=" + read_2_rev_out + " out=" + read_merged_pe_out + " fusepairs=t pad=1 -Xmx" + memory + "g ziplevel=9 threads=" + threads)
		
		print("Unix: Adjusting PE Read Separators.")
		os.system("gunzip " + read_merged_pe_out)
		read_merged_pe_out_naked = read_merged_pe_out.split(".gz")[0]
		os.system("sed -i \'s/ANA/AXA/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/ANC/AXC/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/ANG/AXG/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/ANT/AXT/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/CNA/CXA/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/CNC/CXC/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/CNG/CXG/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/CNT/CXT/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/GNA/GXA/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/GNC/GXC/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/GNT/GXT/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/GNG/GXG/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/TNA/TXA/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/TNC/TXC/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/TNG/TXG/\' " + read_merged_pe_out_naked)
		os.system("sed -i \'s/TNT/TXT/\' " + read_merged_pe_out_naked)
		os.system("gzip " + read_merged_pe_out_naked)
		
		print("Unix: Merging PE and SE Reads.")
		os.system("cat " + read_merged_pe_out + " " + read_3_out + " > " + read_4_out)
		
		print("Unix: Cleaning up output.")
		os.system("rm -r tmp/*")
				
		c = c + 1
		print("")

	for read in read_list_se:
		# Output naming convention
		sample_name = read.split(".filter")[0].split("_")[-1]

		print("fastp: Analyzing " + read + " (SE).")
		
		bbsplit_scaf = "output/bbmap/" + sample_name + "_scaffolds.txt"
		bbsplit_stats = "output/bbmap/" + sample_name + "_stats.txt"
		
		read_in = input_se + read
		read_out = output_dir + sample_name + "_read4.fastq.gz"
		read_tmp = "tmp/" + sample_name + ".reads.final.fastq.gz"
		
		read_html_out = output_dir + "reports/" + sample_name + ".fastp.html"
		read_json_out = output_dir + "reports/" + sample_name + ".fastp.json"
		
		print("fastp: Performing QC.")
		os.system("fastp -i " + read_in + " -o " + read_tmp + " --compression 9 --thread " + threads + " --html " + read_html_out + " --json " + read_json_out + " --overrepresentation_analysis --overrepresentation_sampling 20 --cut_tail --cut_window_size 4 --cut_mean_quality 20")
		
		print("BBMap: Removing Host Genome Contamination.")
		os.system("bbsplit.sh in=" + read_tmp + " ref=" + host + " outu=" + read_out + " threads=" + threads + " path=tmp/ -Xmx" + memory + "g scafstats=" + bbsplit_scaf + " refstats=" + bbsplit_stats)
		
		print("Unix: Cleaning up output.")
		os.system("rm -r tmp/*")

		c = c + 1
		print("")	
		
	print("fastp: " + str(c) + " files successfully analyzed.")
	print("fastp: Finished.\n")

def run_kraken2(database, read_length, bracken_threshold, fastp_dir, bracken_dir, kraken_dir, threads):
	print("Step 3/8 - Taxonomic Classification [kraken2]:\n")
	os.system("mkdir -p " + bracken_dir)
	os.system("mkdir -p " + kraken_dir)
	
	print("kraken2: Computing Taxonomic Classification.")
	read_list = sorted([name for name in os.listdir(fastp_dir) if fnmatch(name, "*_read4.fastq.gz")])
	c = 0
	
	levels = ["Phylum", "Genus", "Species"]
	
	for read in read_list:
		sample_name = read.split("_read4")[0]
		print("kraken2: Analyzing " + read + ".")
		os.system("kraken2" +
				  " --db " + database +
				  " --threads " + threads +
				  " --report " + kraken_dir + sample_name + ".report" +
				  " " + fastp_dir + read +
				  " > " + kraken_dir + sample_name + ".stdout"
				 )
		for level in levels:
			print("Bracken: Analyzing " + read + " (" + level + ").")
			os.system("bracken" +
					  " -d " + database +
					  " -i " + kraken_dir + sample_name + ".report" +
					  " -o " + bracken_dir + sample_name + "." + level.lower() + ".bracken" +
					  " -w " + bracken_dir + sample_name + ".report" +
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
				  " -o " + bracken_dir + "bracken_" + level.lower() + "_all.bracken"
				 )
	
	print("kraken-biom: Exporting bracken reports.")
	os.system("kraken-biom " + bracken_dir + "*.report" +
			  " -o " + bracken_dir + "bracken.biom" +
			  " --fmt " + "json"
			 )
	
	print("kraken2: " + str(c) + " files successfully analyzed.")
	print("kraken2: Finished.\n")

def run_krona(krona_input, krona_output):
	print("Step 4/8 - Visualizing Taxonomic Classification [Krona]:\n")
	os.system("mkdir -p " + krona_output)

	print("Krona: Updating Taxonomic Database.")
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

def run_megahit(megahit_input, megahit_output, threads):
	print("Step 6/8 - Metagenome Assembly [MEGAHIT]:\n")
	
	print("MEGAHIT: Creating Assemblies.")
	read_list = sorted([name for name in os.listdir(megahit_input) if fnmatch(name, "*_read1.fastq.gz")])
	c = 0
	
	# PE + SE
	for read in read_list:
		read_1 = megahit_input + read
		read_2 = megahit_input + read.split("_read1")[0] + "_read2.fastq.gz"
		read_3 = megahit_input + read.split("_read1")[0] + "_read3.fastq.gz"
		sample_name = read.split("_read1")[0]
		print("MEGAHIT: Assembling " + sample_name + ".")
		os.system("megahit" +
				  " -1 " + read_1 +
				  " -2 " + read_2 +
				  " -r " + read_3 +
				  " -m 0.5" +
				  " -t " + threads +
				  " --min-contig-len 150" +
				  " --out-prefix " + sample_name +
				  " --tmp-dir tmp/" +
				  " -o " + megahit_output + "/" + sample_name
				 )
		print("Unix: Cleaning up output.")
		os.system("mv " + megahit_output + "/" + sample_name + "/" + sample_name + ".contigs.fa" + " " +
				  megahit_output + "/" + sample_name + ".fa")
		os.system("rm -r " + megahit_output + "/" + sample_name + "/")
		c = c + 1
		print("")
	
	# SE
	read_se = megahit_input + "SL273100_read4.fastq.gz"
	sample_name = "SL273100"
	print("MEGAHIT: Assembling " + sample_name + ".")
	os.system("megahit" +
			  " -r " + read_se +
			  " -m 0.5" +
			  " -t " + threads +
			  " --min-contig-len 150" +
			  " --out-prefix " + sample_name +
			  " --tmp-dir tmp/" +
			  " -o " + megahit_output + "/" + sample_name
			 )
	print("Unix: Cleaning up output.")
	os.system("mv " + megahit_output + "/" + sample_name + "/" + sample_name + ".contigs.fa " +
			  megahit_output + "/" + sample_name + ".fa")
	os.system("rm -r " + megahit_output + "/" + sample_name + "/")
	c = c + 1
	print("")
	
	print("MEGAHIT: " + str(c) + " files successfully assembled.")
	
	print("Unix: Compressing Assemblies.")
	os.system("gzip " + megahit_output + "/*.fa")
	
	print("MEGAHIT: Finished.\n")

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

def run_metabat(fastp_input, megahit_input, metabat_output, threads):
	print("Step 8/8 - Assembly Binning [MetaBAT]:\n")
	os.system("mkdir -p " + metabat_output)
	
	file_list = glob(megahit_input + "*.gz")
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
				  " -1 " + fastp_input + sample_name + "_read1.fastq.gz" +
				  " -2 " + fastp_input + sample_name + "_read2.fastq.gz" +
				  " -U " + fastp_input + sample_name + "_read2.fastq.gz" +
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
	parser.add_argument("-b", "--bracken_threshold", type = int, default = "10", required = False, help = "Bracken read threshold")
	parser.add_argument("-c", "--city_of_interest", type = str, default = "berlin", required = False, help = "Name of the city of interest (see Metadata Table for options)")
	parser.add_argument("-k", "--kraken2_db", type = str, default = "/scratch1/databases/kraken/20200226_kraken2_standard_database/", required = False, help = "Path to local kraken2 database")
	parser.add_argument("-m", "--memory", type = int, default = "32", required = False, help = "Amount of memory used for running MetaGEN")
	parser.add_argument("-r", "--read_length", type = int, default = "150", required = False, help = "Minimum expected read length")
	parser.add_argument("-s", "--host_species", type = str, default = "human", required = False, help = "Specify the host species of the metagenomic data (human or horse)")
	parser.add_argument("-t", "--threads", type = int, default = "32", required = False, help = "Amount of threads used for running MetaGEN")

	args = parser.parse_args()
	
	print("Running MetaGEN Pipeline Version " + version + "\n")
	
	reference_genome = set_reference(args.host_species)
	#download_metasub(args.city_of_interest)
	#run_fastp(reference_genome, str(args.memory), str(args.threads))
	run_kraken2(args.kraken2_db, str(args.read_length), str(args.bracken_threshold), "output/fastp/", "output/bracken/", "output/kraken2/", str(args.threads))
	#run_krona("output/kraken2/", "output/krona/")
	#run_multiqc("output/fastp/reports/", "output/multiqc/fastp/", "output/kraken2/", "output/multiqc/kraken2/")
	#run_megahit("output/fastp/", "output/megahit/", str(args.threads))
	#run_metaquast("output/megahit/", "output/metaquast/", str(args.threads))
	#run_metabat("output/fastp/", "output/megahit/", "output/metabat/", str(args.threads))
	
	print("MetaGEN: Finished.")
		
if __name__ == "__main__":
	main()