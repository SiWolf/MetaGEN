# -------------------------------
# Title: MetaGEN_AMR_Summarize.py
# Author: Silver A. Wolf
# Last Modified: Mon, 06.09.2021
# Version: 0.0.2
# -------------------------------

from fnmatch import fnmatch, fnmatchcase
import os

kraken_reports = sorted([name for name in os.listdir("output/kraken2/assembly/") if fnmatch(name, "*.stdout")])
kraken_tax_ids = []

# Generate overview of all AMR genes
with open("output/abricate/summary.tab") as areport:
	for line in areport:
		amr_list = [a.strip() for a in line.split("\t")[2:]]
		break

# Generate overview of all required taxids
for k in kraken_reports:
	sample_name = k.split(".stdout")[0]
	contigs = []
	with open("output/abricate/" + sample_name + ".tab") as kreport:
		for line in kreport:
			if line[0] != "#":
				c = line.split("\t")[1]
				if c not in contigs:
					contigs.append(c)
	for c in contigs:
		command = "grep \"" + c + "\" output/kraken2/assembly/" + k
		stdout = os.popen(command).read().split("\t")[2]
		if stdout not in kraken_tax_ids:
			kraken_tax_ids.append(stdout)

# Generate count for each individual report
summary_table = open("output/abricate/kraken2.tab", "w")#
summary_table.write("SAMPLE\t" + "\t".join(kraken_tax_ids) + "\n")
for k in kraken_reports:
	sample_name = k.split(".stdout")[0]
	contigs = []
	tax_ids = []
	tmpstr = ""
	with open("output/abricate/" + sample_name + ".tab") as kreport:
		for line in kreport:
			if line[0] != "#":
				c = line.split("\t")[1]
				contigs.append(c)
	for c in contigs:
		command = "grep \"" + c + "\" output/kraken2/assembly/" + k
		stdout = os.popen(command).read().split("\t")[2]
		tax_ids.append(stdout)
	for t in kraken_tax_ids:
		if t in tax_ids:
			count = tax_ids.count(t)
		else:
			count = 0
		tmpstr = tmpstr + "\t" + str(count)
	summary_table.write(sample_name + tmpstr + "\n")