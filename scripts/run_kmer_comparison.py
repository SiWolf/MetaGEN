# -------------------------------
# Title: run_kmer_comparsion.py
# Author: Silver A. Wolf
# Last Modified: Thu, 30.03.2023
# Version: 0.0.3
# -------------------------------

import os
threads = "216"
kmer_threshold = "254"

# Identification of unique kmers per timepoint

# t0: unique kmers to SSG
os.system("kmc_tools complex scripts/kmer_comparisons/SSG_t0.txt -t " + threads)
os.system("kmc_dump tmp/unique_SSG_t0 tmp/unique_SSG_t0.txt")
os.system("awk '{ if ($2 > " + kmer_threshold + ") { print } }' tmp/unique_SSG_t0.txt > output/03_functional_analysis/kmc3/unique_filter_SSG_t0.txt")
os.system("rm tmp/unique_*")

# t0: unique kmers to 5DG
os.system("kmc_tools complex scripts/kmer_comparisons/5DG_t0.txt -t " + threads)
os.system("kmc_dump tmp/unique_5DG_t0 tmp/unique_5DG_t0.txt")
os.system("awk '{ if ($2 > " + kmer_threshold + ") { print } }' tmp/unique_5DG_t0.txt > output/03_functional_analysis/kmc3/unique_filter_5DG_t0.txt")
os.system("rm tmp/unique_*")

# t1: unique kmers to SSG
os.system("kmc_tools complex scripts/kmer_comparisons/SSG_t1.txt -t " + threads)
os.system("kmc_dump tmp/unique_SSG_t1 tmp/unique_SSG_t1.txt")
os.system("awk '{ if ($2 > " + kmer_threshold + ") { print } }' tmp/unique_SSG_t1.txt > output/03_functional_analysis/kmc3/unique_filter_SSG_t1.txt")
os.system("rm tmp/unique_*")

# t1: unique kmers to 5DG
os.system("kmc_tools complex scripts/kmer_comparisons/5DG_t1.txt -t " + threads)
os.system("kmc_dump tmp/unique_5DG_t1 tmp/unique_5DG_t1.txt")
os.system("awk '{ if ($2 > " + kmer_threshold + ") { print } }' tmp/unique_5DG_t1.txt > output/03_functional_analysis/kmc3/unique_filter_5DG_t1.txt")
os.system("rm tmp/unique_*")

# t2: unique kmers to SSG
os.system("kmc_tools complex scripts/kmer_comparisons/SSG_t2.txt -t " + threads)
os.system("kmc_dump tmp/unique_SSG_t2 tmp/unique_SSG_t2.txt")
os.system("awk '{ if ($2 > " + kmer_threshold + ") { print } }' tmp/unique_SSG_t2.txt > output/03_functional_analysis/kmc3/unique_filter_SSG_t2.txt")
os.system("rm tmp/unique_*")

# t2: unique kmers to 5DG
os.system("kmc_tools complex scripts/kmer_comparisons/5DG_t2.txt -t " + threads)
os.system("kmc_dump tmp/unique_5DG_t2 tmp/unique_5DG_t2.txt")
os.system("awk '{ if ($2 > " + kmer_threshold + ") { print } }' tmp/unique_5DG_t2.txt > output/03_functional_analysis/kmc3/unique_filter_5DG_t2.txt")
os.system("rm tmp/unique_*")