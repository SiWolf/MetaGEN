# -------------------------------
# Title: BLAST_kmer_subset.py
# Author: Silver A. Wolf
# Last Modified: Thu, 30.03.2023
# Version: 0.0.1
# -------------------------------

import glob
import os

kmer_lists = (glob.glob("output/03_functional_analysis/kmc3/unique_filter*"))
threads = 6

for comparison in kmer_lists:
    name = comparison.split("/")[-1].split(".")[0]
    fasta = "output/03_functional_analysis/kmc3/" + name + ".fa"
    fasta_out = open(fasta, "w")
    id = 0
    with open(comparison) as input_list:
        for line in input_list:
            id = id + 1
            fasta_out.write(">" + str(id) + "\n" + line.split("\t")[0].strip() + "\n")
    fasta_out.close()
    blast_out = "output/03_functional_analysis/kmc3/blast_" + name + ".txt"
    #os.system("#makeblastdb -in db/uniprot_sprot.fasta -out db/uniprot_blast_db -dbtype prot")
    os.system("blastp -query " + fasta + " -db db/uniprot_blast_db -outfmt \"6 delim=@ qseqid salltitles\" -max_target_seqs 1 -num_threads " + str(threads) + " > " + blast_out)