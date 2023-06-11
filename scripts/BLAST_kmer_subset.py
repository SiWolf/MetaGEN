# -------------------------------
# Title: BLAST_kmer_subset.py
# Author: Silver A. Wolf
# Last Modified: Sun, 11.06.2023
# Version: 0.0.3
# -------------------------------

import glob
import os
import numpy as np
import multiprocessing
from multiprocessing import Pool

def kmer_blast(kmer_input):
    
    name = kmer_input.split("/")[-1].split(".")[0]
    fasta = "output/03_functional_analysis/kmc3/fasta_" + name + ".fa"
    fasta_out = open(fasta, "w")
    id = 0
    with open(kmer_input) as input_list:
        for line in input_list:
            id = id + 1
            fasta_out.write(">" + str(id) + "\n" + line.split("\t")[0].strip() + "\n")
    fasta_out.close()
    blast_out = "output/03_functional_analysis/kmc3/blast_" + name + ".txt"
    #os.system("#makeblastdb -in db/uniprot_sprot.fasta -out db/uniprot_blast_db -dbtype prot")
    os.system("blastx -query " + fasta + " -db db/uniprot_blast_db -outfmt \"6 delim=@ qseqid salltitles\" -max_target_seqs 1 -num_threads 36 > " + blast_out)

def helperfunction(kmer_list):
	for kmer_file in kmer_list:
		kmer_blast(kmer_file)

# Main
def main():
	input_list = glob.glob("output/03_functional_analysis/kmc3/unique_filter*")

	# Parallelize Script
	amount_cores = 6
	kmer_file_list = np.array_split(input_list, amount_cores)
	pool = Pool(amount_cores)
	results = pool.map(helperfunction, kmer_file_list)

if __name__ == "__main__":
	main()