# -------------------------------
# Title: run_kmer_comparsion.py
# Author: Silver A. Wolf
# Last Modified: Sun, 26.03.2023
# Version: 0.0.1
# -------------------------------

import os

# Identification of unique kmers per timepoint

# t0: unique kmers to SSG
os.system("kmc_tools complex scripts/kmer_comparisons/SSG_t0.txt -t 216")
#kmc_dump output/03_functional_analysis/kmc3/{wildcards.sample} {output.kmc_txt}

# t0: unique kmers to 5DG
#os.system("kmc_tools complex scripts/kmer_comparisons/5DG_t0.txt -t 216")
#kmc_dump output/03_functional_analysis/kmc3/{wildcards.sample} {output.kmc_txt}

# t1: unique kmers to SSG
#os.system("kmc_tools complex scripts/kmer_comparisons/SSG_t1.txt -t 216")
#kmc_dump output/03_functional_analysis/kmc3/{wildcards.sample} {output.kmc_txt}

# t1: unique kmers to 5DG
#os.system("kmc_tools complex scripts/kmer_comparisons/5DG_t1.txt -t 216")
#kmc_dump output/03_functional_analysis/kmc3/{wildcards.sample} {output.kmc_txt}

# t2: unique kmers to SSG
#os.system("kmc_tools complex scripts/kmer_comparisons/SSG_t2.txt -t 216")
#kmc_dump output/03_functional_analysis/kmc3/{wildcards.sample} {output.kmc_txt}

# t2: unique kmers to 5DG
#os.system("kmc_tools complex scripts/kmer_comparisons/5DG_t2.txt -t 216")
#kmc_dump output/03_functional_analysis/kmc3/{wildcards.sample} {output.kmc_txt}