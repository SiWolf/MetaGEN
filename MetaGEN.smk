# -------------------------------
# Title: MetaGEN_Main.smk
# Author: Silver A. Wolf
# Last Modified: Wed, 10.11.2021
# Version: 0.4.0
# -------------------------------

# How to run MetaGEN
#snakemake -s MetaGEN.smk -c 64 --use-conda
#snakemake -F -s MetaGEN.smk -c 64 --use-conda
#snakemake -n -s MetaGEN.smk -c 64 --use-conda
#snakemake --dag -s MetaGEN.smk -c 64 --use-conda | dot -Tsvg > MetaGEN.svg

# -------------------------------
# MetaGEN Settings
# -------------------------------

# Import packages
import os

# Specifying config file
configfile: "config/config.yml"

# Global parameters
(SAMPLES,) = glob_wildcards("input/{sample}_R1.fastq.gz")

# One rule to rule them all
rule all:
	input:
		"output/01_preprocessing/seqfu/stats.tsv",
		"output/01_preprocessing/multiqc/multiqc_report.html",
		"output/02_taxonomic_profiling/multiqc/multiqc_report.html",
		"output/02_taxonomic_profiling/kraken_biom/kraken2.biom",
		"output/02_taxonomic_profiling/bracken/species.summary",
		"output/02_taxonomic_profiling/krona/krona.html",
		expand("output/03_kmer_analysis/kmc3/{sample}.kmc_pre", sample = SAMPLES),
		expand("output/03_kmer_analysis/kmc3/{sample}.kmc_suf", sample = SAMPLES),
		#expand("output/04_assemblies/plasflow/{sample}.txt", sample = SAMPLES)
		"output/04_assemblies/metaquast/report.html",
		expand("output/04_assemblies/metabat/{sample}.fa.gz.depth.txt", sample = SAMPLES),
		expand("output/04_assemblies/metabat/{sample}.fa.gz.paired.txt", sample = SAMPLES),
		"output/05_amr/abricate/kraken2.summary",
		"output/05_amr/coverm/coverm.summary"

# -------------------------------
# V: AMR Profiling
# -------------------------------

# CoverM Summary
rule coverm_summary:
	input:
		coverm_profile = expand("output/05_amr/coverm/{sample}.txt", sample = SAMPLES)
	output:
		coverm_sum = "output/05_amr/coverm/coverm.summary"
	threads:
		1
	message:
		"[CoverM] summarizing results."
	run:
		coverm_results = sorted(input.coverm_profile)
		ref_names = []
		ref_seqs = "MegaRes_ID\tResistance_Type\tResistance_Class\tResistance_Gene\tResistance_Symbol"
		
		for c in coverm_results:
			ref_seqs = ref_seqs + "\t" + c.split(".txt")[0]
		
		with open("output/05_amr/coverm/" + coverm_results[0], "r") as ref_file:
			for line in ref_file:
				ref_names.append(line.split("\t")[0])
		
		ref_names.remove("Contig")
		f = open(output.coverm_sum, "w")
		f.write(ref_seqs + "\n")
		
		for r in ref_names:
			val = []
			for l in coverm_results:
				with open("output/05_amr/coverm/" + l, "r") as cover_file:
					for line in cover_file:
						if line.split("\t")[0] == r:
							val.append(str(int(line.split("\t")[7].replace(".", ",")) + int(line.split("\t")[17].replace(".", ","))))
							break
			f.write("\t".join(r.split("|")[0:5]) + "\t" + "\t".join(val) + "\n")
		f.close()

# CoverM
rule coverm:
	input:
		b1 = "output/01_preprocessing/bbmap/{sample}_R1.fastq.gz",
		b2 = "output/01_preprocessing/bbmap/{sample}_R2.fastq.gz",
		b3 = "output/01_preprocessing/bbmap/{sample}_R3.fastq.gz"
	output:
		coverm_profile = "output/05_amr/coverm/{sample}.txt"
	conda:
		"envs/coverm.yml"
	threads:
		16
	message:
		"[CoverM] mapping reads to the MegaRES database."
	params:
		identity = config["amr_identity"]
		coverage = config["amr_coverage"]
	shell:
		"""
		wget -N -P db/ https://megares.meglab.org/download/megares_v2.00/megares_full_database_v2.00.fasta
		coverm contig -1 {input.b1} -2 {input.b2} --single {input.b3} -r db/megares_full_database_v2.00.fasta \
		-p bwa-mem -m mean trimmed_mean covered_fraction covered_bases variance length count reads_per_base rpkm tpm \
		-o {output.coverm_profile} -t {threads} --bam-file-cache-directory tmp/ --discard-unmapped --exclude-supplementary \
		--min-read-aligned-percent {params.coverage} --min-read-percent-identity {params.identity}
		"""

# ABRicate TaxCaller
rule abricate_taxcaller:
	input:
		a_stdout = expand("output/04_assemblies/kraken2/{sample}.stdout", sample = SAMPLES),
		abricate_profile = expand("output/05_amr/abricate/{sample}.tab", sample = SAMPLES),
		abricate_sum = "output/05_amr/abricate/abricate.summary"
	output:
		kraken2_sum = "output/05_amr/abricate/kraken2.summary"
	threads:
		1
	message:
		"[ABRicate] merging results with taxonomic classifications."
	run:
		kraken_reports = sorted(input.a_stdout)
		kraken_tax_ids = []
		# Generate overview of all AMR genes
		with open(input.abricate_sum) as areport:
			for line in areport:
					amr_list = [a.strip() for a in line.split("\t")[2:]]
					break
		# Generate overview of all required taxids
		for k in kraken_reports:
			sample_name = k.split(".stdout")[0]
			contigs = []
			with open("output/05_amr/abricate/" + sample_name + ".tab") as kreport:
				for line in kreport:
					if line[0] != "#":
						c = line.split("\t")[1]
						if c not in contigs:
							contigs.append(c)
			for c in contigs:
				command = "grep \"" + c + "\" output/04_assemblies/kraken2/" + k
				stdout = os.popen(command).read().split("\t")[2]
				if stdout not in kraken_tax_ids:
					kraken_tax_ids.append(stdout)
		# Generate count for each individual report
		summary_table = open(output.kraken2_sum, "w")
		summary_table.write("SAMPLE\t" + "\t".join(kraken_tax_ids) + "\n")
		for k in kraken_reports:
			sample_name = k.split(".stdout")[0]
			contigs = []
			tax_ids = []
			tmpstr = ""
			with open("output/05_amr/abricate/" + sample_name + ".tab") as kreport:
				for line in kreport:
					if line[0] != "#":
						c = line.split("\t")[1]
						contigs.append(c)
			for c in contigs:
				command = "grep \"" + c + "\" output/04_assemblies/kraken2/" + k
				stdout = os.popen(command).read().split("\t")[2]
				tax_ids.append(stdout)
			for t in kraken_tax_ids:
				if t in tax_ids:
					count = tax_ids.count(t)
				else:
					count = 0
				tmpstr = tmpstr + "\t" + str(count)
			summary_table.write(sample_name + tmpstr + "\n")

# ABRicate Summary
rule abricate_summary:
	input:
		renamed = expand("output/05_amr/abricate/{sample}.tab", sample = SAMPLES)
	output:
		abricate_sum = "output/05_amr/abricate/abricate.summary"
	conda:
		"envs/abricate.yml"
	threads:
		1
	message:
		"[ABRicate] summarizing results."
	shell:
		"""
		abricate --summary --nopath output/05_amr/abricate/*.tab > {output.abricate_sum}
		"""

# ABRicate
rule abricate:
	input:
		renamed = "output/04_assemblies/bbmap/{sample}.fa.gz"
	output:
		abricate_profile = "output/05_amr/abricate/{sample}.tab"
	conda:
		"envs/abricate.yml"
	threads:
		16
	message:
		"[ABRicate] scanning for AMR Genes."
	params:
		identity = config["amr_identity"]
		coverage = config["amr_coverage"]
	shell:
		"""
		abricate --db megares --threads {threads} --minid {params.identity} --mincov {params.coverage} --nopath {input.renamed} > {output.abricate_profile}
		"""

# -------------------------------
# IV: Assemblies
# -------------------------------

# MetaBAT
rule metabat:
	input:
		b1 = "output/01_preprocessing/bbmap/{sample}_R1.fastq.gz",
		b2 = "output/01_preprocessing/bbmap/{sample}_R2.fastq.gz",
		b3 = "output/01_preprocessing/bbmap/{sample}_R3.fastq.gz",
		renamed = "output/04_assemblies/bbmap/{sample}.fa.gz"
	output:
		bin_depth = "output/04_assemblies/metabat/{sample}.fa.gz.depth.txt",
		bin_paired = "output/04_assemblies/metabat/{sample}.fa.gz.paired.txt"
	conda:
		"envs/metabat.yml"
	threads:
		32
	message:
		"[MetaBAT] binning assembly of {wildcards.sample}."
	shell:
		"""
		bowtie2-build --quiet --threads {threads} {input.renamed} tmp/{wildcards.sample}
		bowtie2 --quiet -p {threads} -x tmp/{wildcards.sample} -1 {input.b1} -2 {input.b2} -U {input.b3} -S tmp/{wildcards.sample}.sam
		samtools view -bS -o tmp/{wildcards.sample}.bam tmp/{wildcards.sample}.sam
		samtools sort tmp/{wildcards.sample}.bam -o tmp/{wildcards.sample}_sorted.bam
		samtools index tmp/{wildcards.sample}_sorted.bam
		runMetaBat.sh -m 1500 -t {threads} renamed {wildcards.sample}_sorted.bam
		mv {wildcards.sample}* output/04_assemblies/metabat/
		"""

# MetaQUAST
rule metaquast:
	input:
		renamed = expand("output/04_assemblies/bbmap/{sample}.fa.gz", sample = SAMPLES)
	output:
		qc_assembly = "output/04_assemblies/metaquast/report.html"
	conda:
		"envs/metaquast.yml"
	threads:
		32
	message:
		"[MetaQUAST] assessing quality of assemblies."
	shell:
		"""
		metaquast -o output/04_assemblies/metaquast/ -t {threads} {input.renamed} --glimmer --rna-finding --plots-format png --silent -m 100
		"""

#PlasFlow
#rule plasflow:
#	input:
#		renamed = "output/04_assemblies/bbmap/{sample}.fa.gz"
#	output:
#		plasmids = "output/04_assemblies/plasflow/{sample}.txt"
#	conda:
#		"envs/plasflow.yml"
#	threads:
#		1
#	message:
#		"[PlasFlow] detecting plasmid sequences in the assembly of {wildcards.sample}."
#	shell:
#		"""
#		PlasFlow.py --input {input.renamed} --output {output.plasmids}
#		"""

# kraken2
rule kraken2_assembly:
	input:
		renamed = "output/04_assemblies/bbmap/{sample}.fa.gz"
	output:
		a_report = "output/04_assemblies/kraken2/{sample}.report",
		a_stdout = "output/04_assemblies/kraken2/{sample}.stdout"
	conda:
		"envs/kraken2.yml"
	threads:
		32
	message:
		"[kraken2] assessing taxonomic content of assembly for {wildcards.sample}."
	params:
		db = config["kraken_db"]
	shell:
		"""
		kraken2 --db {params.db} --threads {threads} --report {output.a_report} --output {output.a_stdout} {input.renamed}
		"""

# bbmap rename
rule bbmap_rename:
	input:
		assembly = "output/04_assemblies/metaspades/{sample}.fa.gz"
	output:
		renamed = "output/04_assemblies/bbmap/{sample}.fa.gz"
	conda:
		"envs/bbmap.yml"
	threads:
		16
	message:
		"[bbmap] renaming contigs of {wildcards.sample}."
	params:
		min_length = config["assembly_min"]
	shell:
		"""
		rename.sh in={input.assembly} out={output.renamed} prefix={wildcards.sample} zl=9 minscaf={params.min_length} -Xmx{threads}g
		"""

# metaspades
rule metaspades:
	input:
		b1 = "output/01_preprocessing/bbmap/{sample}_R1.fastq.gz",
		b2 = "output/01_preprocessing/bbmap/{sample}_R2.fastq.gz",
		b3 = "output/01_preprocessing/bbmap/{sample}_R3.fastq.gz"
	output:
		assembly = "output/04_assemblies/metaspades/{sample}.fa.gz"
	conda:
		"envs/metaspades.yml"
	threads:
		64
	message:
		"[metaSPAdes] assembling {wildcards.sample}."
	shell:
		"""
		spades.py -o output/04_assemblies/metaspades/{wildcards.sample}/ --meta -1 {input.b1} -2 {input.b2} -s {input.b3} -t {threads}
		mv output/04_assemblies/metaspades/{wildcards.sample}/scaffolds.fasta output/04_assemblies/metaspades/{wildcards.sample}.fa
		rm -r output/04_assemblies/metaspades/{wildcards.sample}/
		gzip output/04_assemblies/metaspades/{wildcards.sample}.fa
		"""

# -------------------------------
# III: K-mer Analysis
# -------------------------------

# KMC3
rule kmc3:
	input:
		b1 = "output/01_preprocessing/bbmap/{sample}_R1.fastq.gz",
		b2 = "output/01_preprocessing/bbmap/{sample}_R2.fastq.gz",
		b3 = "output/01_preprocessing/bbmap/{sample}_R3.fastq.gz"
	output:
		kmc_pre = "output/03_kmer_analysis/kmc3/{sample}.kmc_pre",
		kmc_suf = "output/03_kmer_analysis/kmc3/{sample}.kmc_suf"
	conda:
		"envs/kmc.yml"
	threads:
		64
	message:
		"[KMC3] calculating k-mer statistics for {wildcards.sample}."
	shell:
		"""
		echo {input.b1} > tmp/sample_list.txt
		echo {input.b2} >> tmp/sample_list.txt
		echo {input.b3} >> tmp/sample_list.txt
		kmc @tmp/sample_list.txt output/03_kmer_analysis/kmc3/{wildcards.sample} tmp/ -m100 -sm -fq -ci0 -cs999 -t {threads}
		"""

# -------------------------------
# II: Taxonomic Profiling
# -------------------------------

# Krona
rule krona:
	input:
		k_stdout = expand("output/02_taxonomic_profiling/kraken2/{sample}.stdout", sample = SAMPLES)
	output:
		krona_html = "output/02_taxonomic_profiling/krona/krona.html"
	conda:
		"envs/krona.yml"
	threads:
		1
	message:
		"[Krona] visualizing taxonomic compositions."
	shell:
		"""
		ktUpdateTaxonomy.sh
		ktImportTaxonomy -q 2 -t 3 {input.k_stdout} -o {output.krona_html}
		"""

# kraken-biom
rule kraken_biom:
	input:
		b_report = expand("output/02_taxonomic_profiling/bracken/{sample}.report", sample = SAMPLES),
		k_report = expand("output/02_taxonomic_profiling/kraken2/{sample}.report", sample = SAMPLES)
	output:
		b_summary = "output/02_taxonomic_profiling/kraken_biom/bracken.biom",
		k_summary = "output/02_taxonomic_profiling/kraken_biom/kraken2.biom",
		k_taxids = "output/02_taxonomic_profiling/kraken_biom/taxids.txt"
	conda:
		"envs/kraken_biom.yml"
	threads:
		1
	message:
		"[kraken-biom] converting kraken and bracken reports to biom format."
	shell:
		"""
		kraken-biom output/02_taxonomic_profiling/bracken/*.report -o {output.b_summary} --fmt json --max D --min S
		kraken-biom output/02_taxonomic_profiling/kraken2/*.report -o {output.k_summary} --fmt json --max D --min S --otu_fp {output.k_taxids}
		"""

# bracken summary
rule bracken_summary:
	input:
		b_report = expand("output/02_taxonomic_profiling/bracken/{sample}.bracken", sample = SAMPLES)
	output:
		b_merged = "output/02_taxonomic_profiling/bracken/species.summary"
	conda:
		"envs/bracken.yml"
	threads:
		1
	message:
		"[bracken] merging bracken reports."
	shell:
		"""
		combine_bracken_outputs.py --files output/02_taxonomic_profiling/bracken/*.bracken -o {output.b_merged}
		"""

# bracken abundancies
rule bracken_abundancies:
	input:
		k_report = "output/02_taxonomic_profiling/kraken2/{sample}.report"
	output:
		b_report = "output/02_taxonomic_profiling/bracken/{sample}.report",
		b_output = "output/02_taxonomic_profiling/bracken/{sample}.bracken"
	conda:
		"envs/bracken.yml"
	threads:
		1
	message:
		"[bracken] re-estimating species abundancies for {wildcards.sample}."
	params:
		db = config["kraken_db"]
		length = config["kraken_read_length"]
	shell:
		"""
		bracken -d {params.db} -i {input.k_report} -o {output.b_output} -w {output.b_report} -r {params.length} -l S
		"""

# kraken2
rule kraken2_reads:
	input:
		b1 = "output/01_preprocessing/bbmap/{sample}_R1.fastq.gz",
		b2 = "output/01_preprocessing/bbmap/{sample}_R2.fastq.gz"
	output:
		k_report = "output/02_taxonomic_profiling/kraken2/{sample}.report",
		k_stdout = "output/02_taxonomic_profiling/kraken2/{sample}.stdout"
	conda:
		"envs/kraken2.yml"
	threads:
		32
	message:
		"[kraken2] assessing taxonomic profile for {wildcards.sample}."
	params:
		db = config["kraken_db"]
	shell:
		"""
		kraken2 --db {params.db} --threads {threads} --report {output.k_report} --output {output.k_stdout} --paired {input.b1} {input.b2}
		"""

# -------------------------------
# I: Quality Control
# -------------------------------

# MultiQC
rule multiqc:
	input:
		json = expand("output/01_preprocessing/fastp/reports/{sample}.fastp.json", sample = SAMPLES),
		k_report = expand("output/02_taxonomic_profiling/kraken2/{sample}.report", sample = SAMPLES)
	output:
		qc_fastp = "output/01_preprocessing/multiqc/multiqc_report.html",
		qc_kraken2 = "output/02_taxonomic_profiling/multiqc/multiqc_report.html"
	conda:
		"envs/multiqc.yml"
	threads:
		1
	message:
		"[multiqc] summarizing fastp and kraken2 reports."
	shell:
		"""
		multiqc -o output/01_preprocessing/multiqc/ output/01_preprocessing/fastp/reports/* -q -z
		multiqc -o output/02_taxonomic_profiling/multiqc/ output/02_taxonomic_profiling/kraken2/* -q -z
		"""

# SeqFu
rule seqfu:
	input:
		b1 = expand("output/01_preprocessing/bbmap/{sample}_R1.fastq.gz", sample = SAMPLES),
		b2 = expand("output/01_preprocessing/bbmap/{sample}_R2.fastq.gz", sample = SAMPLES),
		b3 = expand("output/01_preprocessing/bbmap/{sample}_R3.fastq.gz", sample = SAMPLES)
	output:
		stats = "output/01_preprocessing/seqfu/stats.tsv"
	conda:
		"envs/seqfu.yml"
	threads:
		1
	message:
		"[seqfu] calculating read statistics."
	shell:
		"""
		seqfu stats output/01_preprocessing/bbmap/* > {output.stats}
		"""

# bbmap split
rule bbmap_split:
	input:
		f1 = "output/01_preprocessing/fastp/{sample}_R1.fastq.gz",
		f2 = "output/01_preprocessing/fastp/{sample}_R2.fastq.gz",
		f3 = "output/01_preprocessing/fastp/{sample}_R3.fastq.gz"
	output:
		b1 = "output/01_preprocessing/bbmap/{sample}_R1.fastq.gz",
		b2 = "output/01_preprocessing/bbmap/{sample}_R2.fastq.gz",
		b3 = "output/01_preprocessing/bbmap/{sample}_R3.fastq.gz",
		ref12 = "output/01_preprocessing/bbmap/reports/{sample}_refstats_12.txt",
		ref03 = "output/01_preprocessing/bbmap/reports/{sample}_refstats_3.txt",
		scaf12 = "output/01_preprocessing/bbmap/reports/{sample}_scaffolds_12.txt",
		scaf03 = "output/01_preprocessing/bbmap/reports/{sample}_scaffolds_3.txt"
	conda:
		"envs/bbmap.yml"
	threads:
		64
	message:
		"[bbmap] removing host contamination for {wildcards.sample}."
	shell:
		"""
		bbsplit.sh in1={input.f1} in2={input.f2} outu1={output.b1} outu2={output.b2} path=references/ threads={threads} refstats={output.ref12} scafstats={output.scaf12} -Xmx{threads}g
		bbsplit.sh in={input.f3} outu={output.b3} path=references/ threads={threads} refstats={output.ref03} scafstats={output.scaf03} -Xmx{threads}g
		"""

# fastp
rule fastp:
	input:
		r1 = "input/{sample}_R1.fastq.gz",
		r2 = "input/{sample}_R2.fastq.gz"
	output:
		f1 = "output/01_preprocessing/fastp/{sample}_R1.fastq.gz",
		f2 = "output/01_preprocessing/fastp/{sample}_R2.fastq.gz",
		f3 = "output/01_preprocessing/fastp/{sample}_R3.fastq.gz",
		html = "output/01_preprocessing/fastp/reports/{sample}.fastp.html",
		json = "output/01_preprocessing/fastp/reports/{sample}.fastp.json"
	conda:
		"envs/fastp.yml"
	threads:
		16
	message:
		"[fastp] trimming reads of {wildcards.sample}."
	shell:
		"""
		fastp --in1 {input.r1} --out1 {output.f1} --in2 {input.r2} --out2 {output.f2} \
		--unpaired1 {output.f3} --unpaired2 {output.f3} --html {output.html} --json {output.json} \
		--thread {threads} --compression 9 --correction --cut_tail --cut_window_size 4 \
		--cut_mean_quality 20 --detect_adapter_for_pe --overrepresentation_analysis \
		--overrepresentation_sampling 20
		"""