# -------------------------------
# Title: plasmid_extraction.py
# Author: Silver A. Wolf
# Last Modified: Tue, 07.04.2023
# Version: 0.0.2
# -------------------------------

import glob
plas_threshold = 0.5

def summarize_mobilome(name, path):
	
	print("Summarizing mobilome (" + name + "):\n")

	abricate = glob.glob(path)
	
	mobile_amr_counts = open("output/08_visualization/tab_mobile_" + name + "_counts.tsv", "w")
	mobile_amr_list = open("output/08_visualization/tab_mobile_"  + name + "_list.tsv", "w")
	
	mobile_amr_counts.write("sample" + "\t" + "mobile_arg_counts" + "\n")
	mobile_amr_list.write("sample" + "\t" + "contig" + "\t" + "amr_gene" + "\t" + "mobile_likelihood" + "\t" + "taxid" + "\n")
	
	for amr_tab in abricate:
		amr_contig = ""
		amr_gene = ""
		plas_contig = ""
		plas_likelihood = 0
		kraken_contig = ""
		kraken_taxid = ""
		mobile_amr_count = 0
		sample = amr_tab.split("/")[-1].split(".")[0]

		print("Fetching results from " + sample + "...")

		with open(amr_tab) as amr_report:
			for amr_line in amr_report:
				amr_contig = amr_line.split("\t")[1].strip()
				amr_gene = amr_line.split("\t")[5].strip()
			
				plas_tab = "output/07_amr/plasclass/" + sample + ".txt"
				with open(plas_tab) as plas_report:
					for plas_line in plas_report:
						plas_contig = plas_line.split("\t")[0].strip()
						plas_likelihood = float(plas_line.split("\t")[1].strip())
					
						if plas_contig == amr_contig:
							if plas_likelihood > plas_threshold:
							
								kraken_tab = "output/04_assemblies/kraken2/" + sample + ".stdout"
								with open(kraken_tab) as kraken_report:
									for kraken_line in kraken_report:
										kraken_contig = kraken_line.split("\t")[1].strip()
										if kraken_contig == plas_contig:
											kraken_taxid = kraken_line.split("\t")[2].strip()
											break

								mobile_amr_list.write(sample + "\t" + plas_contig + "\t" + amr_gene + "\t" + str(plas_likelihood) + "\t" + kraken_taxid + "\n")
								mobile_amr_count = mobile_amr_count + 1
							break
	
		mobile_amr_counts.write(sample + "\t" + str(mobile_amr_count) + "\n")
		print("Finished fetching results from " + sample + ".\n")

	mobile_amr_counts.close()
	mobile_amr_list.close()

# Main
def main():
	summarize_mobilome("arg","output/07_amr/abricate/amr/*.tab")
	summarize_mobilome("vir","output/07_amr/abricate/vir/*.tab")

if __name__ == "__main__":
	main()