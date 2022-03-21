# -------------------------------
# Title: MetaGEN_Main.smk
# Author: Silver A. Wolf
# Last Modified: Fri, 18.03.2022
# Version: 0.5.6
# -------------------------------

# How to run MetaGEN
#snakemake -s MetaGEN.smk -c 128 --use-conda
#snakemake -F -s MetaGEN.smk -c 128 --use-conda
#snakemake -n -s MetaGEN.smk -c 128 --use-conda
#snakemake --dag -s MetaGEN.smk -c 128 --use-conda | dot -Tsvg > MetaGEN.svg

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
		"output/02_taxonomic_profiling/kraken_biom/bracken.biom",
		"output/02_taxonomic_profiling/krona/krona.html",
		expand("output/03_kmer_analysis/kmc3/{sample}.kmc_pre", sample = SAMPLES),
		expand("output/03_kmer_analysis/kmc3/{sample}.kmc_suf", sample = SAMPLES),
		expand("output/04_assemblies/plasclass/{sample}.txt", sample = SAMPLES),
		expand("output/04_assemblies/metaquast/{sample}/report.html", sample = SAMPLES),
		expand("output/05_genomic_bins/checkm/{sample}.txt", sample = SAMPLES),
		"output/06_co_assembly/checkm/co_assembly.txt",
		"output/07_amr/abricate/amr/kraken2.summary",
		"output/07_amr/coverm/coverm.summary"

# -------------------------------
# VII: AMR & Virulence Profiling
# -------------------------------

# CoverM Summary
rule coverm_summary:
	input:
		coverm_profile = expand("output/07_amr/coverm/{sample}.txt", sample = SAMPLES)
	output:
		coverm_sum = "output/07_amr/coverm/coverm.summary"
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
		
		with open(coverm_results[0], "r") as ref_file:
			for line in ref_file:
				ref_names.append(line.split("\t")[0])
		
		ref_names.remove("Contig")
		f = open(output.coverm_sum, "w")
		f.write(ref_seqs + "\n")
		
		for r in ref_names:
			val = []
			for l in coverm_results:
				with open(l, "r") as cover_file:
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
		coverm_profile = "output/07_amr/coverm/{sample}.txt"
	conda:
		"envs/coverm.yml"
	threads:
		16
	message:
		"[CoverM] mapping reads to the MegaRES database."
	params:
		identity = config["amr_identity"],
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
		abricate_profile = expand("output/07_amr/abricate/amr/{sample}.tab", sample = SAMPLES),
		abricate_sum = "output/07_amr/abricate/amr/abricate.summary"
	output:
		kraken2_sum = "output/07_amr/abricate/amr/kraken2.summary"
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
			sample_name = k.split(".stdout")[0].split("/")[-1]
			contigs = []
			with open("output/07_amr/abricate/amr/" + sample_name + ".tab") as kreport:
				for line in kreport:
					if line[0] != "#":
						c = line.split("\t")[1]
						if c not in contigs:
							contigs.append(c)
			for c in contigs:
				command = "grep \"" + c + "\" " + k
				stdout = os.popen(command).read().split("\t")[2]
				if stdout not in kraken_tax_ids:
					kraken_tax_ids.append(stdout)
		# Generate count for each individual report
		summary_table = open(output.kraken2_sum, "w")
		summary_table.write("SAMPLE\t" + "\t".join(kraken_tax_ids) + "\n")
		for k in kraken_reports:
			sample_name = k.split(".stdout")[0].split("/")[-1]
			contigs = []
			tax_ids = []
			tmpstr = ""
			with open("output/07_amr/abricate/amr/" + sample_name + ".tab") as kreport:
				for line in kreport:
					if line[0] != "#":
						c = line.split("\t")[1]
						contigs.append(c)
			for c in contigs:
				command = "grep \"" + c + "\" " + k
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
		amr_profile = expand("output/07_amr/abricate/amr/{sample}.tab", sample = SAMPLES),
		vir_profile = expand("output/07_amr/abricate/vir/{sample}.tab", sample = SAMPLES)
	output:
		amr_sum = "output/07_amr/abricate/amr/abricate.summary",
		vir_sum = "output/07_amr/abricate/vir/abricate.summary"
	conda:
		"envs/abricate.yml"
	threads:
		1
	message:
		"[ABRicate] summarizing results."
	shell:
		"""
		abricate --summary --nopath output/07_amr/abricate/amr/*.tab > {output.amr_sum}
		abricate --summary --nopath output/07_amr/abricate/vir/*.tab > {output.vir_sum}
		"""

# ABRicate
rule abricate:
	input:
		renamed = "output/04_assemblies/bbmap/{sample}.fa"
	output:
		amr_profile = "output/07_amr/abricate/amr/{sample}.tab",
		vir_profile = "output/07_amr/abricate/vir/{sample}.tab"
	conda:
		"envs/abricate.yml"
	threads:
		16
	message:
		"[ABRicate] scanning for AMR and virulence genes."
	params:
		identity = config["amr_identity"],
		coverage = config["amr_coverage"]
	shell:
		"""
		abricate --db megares --threads {threads} --minid {params.identity} --mincov {params.coverage} --nopath {input.renamed} > {output.amr_profile}
		abricate --db vfdb --threads {threads} --minid {params.identity} --mincov {params.coverage} --nopath {input.renamed} > {output.vir_profile}
		"""

# -------------------------------
# VI: Co-Assembly
# -------------------------------

# Bakta
rule co_assembly_bakta:
	input:
		renamed = "output/06_co_assembly/bbmap/co_assembly.fa"
	output:
		annotation = "output/06_co_assembly/bakta/co_assembly.gff3"
	conda:
		"envs/checkm.yml"
	threads:
		64
	message:
		"[Bakta] Annotating co-assembly."
	params:
		bakta_db = config["kraken_db"]
	shell:
		"""
		bakta --db {params.bakta_db} -o output/06_co_assembly/bakta/ -p co_assembly --threads {threads} {input.renamed}
		"""

# CheckM
rule co_assembly_checkm:
	input:
		bins = "output/06_co_assembly/metabat/bin/bin.1.fa",
		co_depth = "output/06_co_assembly/metabat/depth.txt",
		co_paired = "output/06_co_assembly/metabat/paired.txt"
	output:
		co_checkm = "output/06_co_assembly/checkm/co_assembly.txt"
	conda:
		"envs/checkm.yml"
	threads:
		32
	message:
		"[CheckM] Assessing genomic bin quality of co-assembly."
	shell:
		"""
		checkm lineage_wf -t {threads} -f {output.co_checkm} -x fa output/06_co_assembly/metabat/bin/ output/06_co_assembly/checkm/
		"""

# MetaBAT
rule co_assembly_metabat:
	input:
		bam = expand("output/06_co_assembly/bowtie2/{sample}.bam", sample = SAMPLES),
		renamed = "output/06_co_assembly/bbmap/co_assembly.fa"
	output:
		bins = "output/06_co_assembly/metabat/bin/bin.1.fa",
		co_depth = "output/06_co_assembly/metabat/depth.txt",
		co_paired = "output/06_co_assembly/metabat/paired.txt"
	conda:
		"envs/metabat.yml"
	threads:
		128
	message:
		"[MetaBAT] binning assembly of co-assembly."
	params:
		min_depth = config["assembly_depth"],
		min_length = config["assembly_min"]
	shell:
		"""
		jgi_summarize_bam_contig_depths --outputDepth {output.depth} --pairedContigs {output.paired} --minContigLength {params.min_length} --minContigDepth {params.min_depth} output/06_co_assembly/bowtie2/*.bam
		metabat2 -m 1500 -a {output.depth} -i {input.renamed} -o output/06_co_assembly/metabat/bin/bin -t {threads}
		"""

# Bowtie 2
rule co_assembly_bowtie2:
	input:
		b1 = "output/01_preprocessing/bbmap/{sample}_R1.fastq.gz",
		b2 = "output/01_preprocessing/bbmap/{sample}_R2.fastq.gz",
		b3 = "output/01_preprocessing/bbmap/{sample}_R3.fastq.gz",
		index = "tmp/co_assembly.1.bt2l"
	output:
		bam = "output/06_co_assembly/bowtie2/{sample}.bam"
	conda:
		"envs/metabat.yml"
	threads:
		32
	message:
		"[Bowtie 2] Mapping reads of {wildcards.sample} to co-assembly."
	shell:
		"""
		bowtie2 --quiet --no-unal -p {threads} -x tmp/co_assembly -1 {input.b1} -2 {input.b2} -U {input.b3} -S tmp/co-{wildcards.sample}.sam
		samtools view -bS -o tmp/co-{wildcards.sample}.bam tmp/co-{wildcards.sample}.sam
		samtools sort tmp/co-{wildcards.sample}.bam -o output/06_co_assembly/bowtie2/{wildcards.sample}.bam
		samtools index output/06_co_assembly/bowtie2/{wildcards.sample}.bam
		"""

# Bowtie 2
rule co_assembly_bowtie2_index:
	input:
		renamed = "output/06_co_assembly/bbmap/co_assembly.fa"
	output:
		index = "tmp/co_assembly.1.bt2l"
	conda:
		"envs/metabat.yml"
	threads:
		64
	message:
		"[Bowtie 2] Generating index for co-assembly."
	shell:
		"""
		bowtie2-build --quiet --threads {threads} {input.renamed} tmp/co_assembly
		"""

# bbmap reformat
rule co_assembly_bbmap:
	input:
		co_assembly = "output/06_co_assembly/megahit/co_assembly.fa.gz"
	output:
		renamed = "output/06_co_assembly/bbmap/co_assembly.fa"
	conda:
		"envs/bbmap.yml"
	threads:
		16
	message:
		"[bbmap] renaming contigs of the co-assembly."
	shell:
		"""
		rename.sh in={input.co_assembly} out={output.renamed} prefix=MEGA -Xmx{threads}g
		"""

# MEGAHIT
rule co_assembly_megahit:
	input:
		b1 = expand("output/01_preprocessing/bbmap/{sample}_R1.fastq.gz", sample = SAMPLES),
		b2 = expand("output/01_preprocessing/bbmap/{sample}_R2.fastq.gz", sample = SAMPLES),
		b3 = expand("output/01_preprocessing/bbmap/{sample}_R3.fastq.gz", sample = SAMPLES)
	output:
		co_assembly = "output/06_co_assembly/megahit/co_assembly.fa.gz"
	conda:
		"envs/megahit.yml"
	threads:
		128
	message:
		"[MEGAHIT] Performing co-assembly."
	params:
		min_length = config["assembly_min"]
	shell:
		"""
		echo {input.b1} > tmp/b1.txt
		echo {input.b2} > tmp/b2.txt
		echo {input.b3} > tmp/b3.txt
		xb1=`cat tmp/b1.txt`
		xb2=`cat tmp/b2.txt`
		xb3=`cat tmp/b3.txt`
		yb1=${{xb1// /,}}
		yb2=${{xb2// /,}}
		yb3=${{xb3// /,}}
		megahit -1 "$yb1" -2 "$yb2" -r "$yb3" --kmin-1pass --k-list 27,37,47,57,67,77,87 --min-contig-len {params.min_length} -t {threads} -o tmp/co_assembly/
		mv tmp/co_assembly/final.contigs.fa output/06_co_assembly/megahit/co_assembly.fa
		gzip output/06_co_assembly/megahit/co_assembly.fa
		rm -r tmp/co_assembly/
		"""

# -------------------------------
# V: Genomic Binning
# -------------------------------

# CheckM
rule checkm:
	input:
		fasta_bins = "output/05_genomic_bins/metabat/{sample}.fa.metabat-bins32/bin.1.fa",
		stats_depth = "output/05_genomic_bins/metabat/{sample}.fa.depth.txt",
		stats_paired = "output/05_genomic_bins/metabat/{sample}.fa.paired.txt"
	output:
		checkm = "output/05_genomic_bins/checkm/{sample}.txt"
	conda:
		"envs/checkm.yml"
	threads:
		16
	message:
		"[CheckM] Assessing genomic bin quality of {wildcards.sample}."
	shell:
		"""
		checkm lineage_wf -t {threads} -f {output.checkm} -x fa output/05_genomic_bins/metabat/{wildcards.sample}.fa.metabat-bins32/ output/05_genomic_bins/checkm/{wildcards.sample}/
		"""

# MetaBAT
rule metabat:
	input:
		b1 = "output/01_preprocessing/bbmap/{sample}_R1.fastq.gz",
		b2 = "output/01_preprocessing/bbmap/{sample}_R2.fastq.gz",
		b3 = "output/01_preprocessing/bbmap/{sample}_R3.fastq.gz",
		filtered = "output/04_assemblies/bbmap/filtered/{sample}.fa"
	output:
		fasta_bins = "output/05_genomic_bins/metabat/{sample}.fa.metabat-bins32/bin.1.fa",
		stats_depth = "output/05_genomic_bins/metabat/{sample}.fa.depth.txt",
		stats_paired = "output/05_genomic_bins/metabat/{sample}.fa.paired.txt"
	conda:
		"envs/metabat.yml"
	threads:
		32
	message:
		"[MetaBAT] binning assembly of {wildcards.sample}."
	shell:
		"""
		bowtie2-build --quiet --threads {threads} {input.filtered} tmp/{wildcards.sample}
		bowtie2 --quiet --no-unal -p {threads} -x tmp/{wildcards.sample} -1 {input.b1} -2 {input.b2} -U {input.b3} -S tmp/{wildcards.sample}.sam
		samtools view -bS -o tmp/{wildcards.sample}.bam tmp/{wildcards.sample}.sam
		samtools sort tmp/{wildcards.sample}.bam -o tmp/{wildcards.sample}_sorted.bam
		samtools index tmp/{wildcards.sample}_sorted.bam
		runMetaBat.sh -m 1500 -t {threads} {input.filtered} tmp/{wildcards.sample}_sorted.bam
		mv -f {wildcards.sample}* output/05_genomic_bins/metabat/
		"""

# -------------------------------
# IV: Assemblies
# -------------------------------

# MetaQUAST
rule metaquast:
	input:
		renamed = "output/04_assemblies/bbmap/filtered/{sample}.fa"
	output:
		qc_assembly = "output/04_assemblies/metaquast/{sample}/report.html"
	conda:
		"envs/metaquast.yml"
	threads:
		128
	message:
		"[MetaQUAST] assessing quality of assemblies."
	params:
		max_refs = config["metaquast_max_refs"]
	shell:
		"""
		metaquast -o output/04_assemblies/metaquast/{wildcards.sample}/ -t {threads} {input.renamed} --plots-format png --silent --max-ref-number {params.max_refs}
		"""

#PlasClass
rule plasclass:
	input:
		filtered = "output/04_assemblies/bbmap/filtered/{sample}.fa"
	output:
		plasmids = "output/04_assemblies/plasclass/{sample}.txt"
	conda:
		"envs/plasclass.yml"
	threads:
		32
	message:
		"[PlasClass] detecting plasmid sequences in the assembly of {wildcards.sample}."
	shell:
		"""
		classify_fasta.py -f {input.filtered} -o {output.plasmids} -p {threads}
		"""

# kraken2
rule kraken2_assembly:
	input:
		renamed = "output/04_assemblies/bbmap/{sample}.fa"
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
		db = config["kraken_db"],
		confidence = config["kraken_confidence_score"]
	shell:
		"""
		kraken2 --db {params.db} --threads {threads} --confidence {params.confidence} --report {output.a_report} --output {output.a_stdout} {input.renamed}
		"""

# bbmap reformat
rule bbmap_reformat:
	input:
		assembly = "output/04_assemblies/megahit/{sample}.fa.gz"
	output:
		filtered = "output/04_assemblies/bbmap/filtered/{sample}.fa",
		renamed = "output/04_assemblies/bbmap/{sample}.fa"
	conda:
		"envs/bbmap.yml"
	threads:
		16
	message:
		"[bbmap] renaming and filtering contigs of {wildcards.sample}."
	params:
		min_length = config["assembly_min"]
	shell:
		"""
		rename.sh in={input.assembly} out={output.renamed} prefix={wildcards.sample} -Xmx{threads}g
		reformat.sh in={output.renamed} out={output.filtered} minlength={params.min_length} -Xmx{threads}g
		"""

# MEGAHIT
rule megahit:
	input:
		b1 = "output/01_preprocessing/bbmap/{sample}_R1.fastq.gz",
		b2 = "output/01_preprocessing/bbmap/{sample}_R2.fastq.gz",
		b3 = "output/01_preprocessing/bbmap/{sample}_R3.fastq.gz"
	output:
		assembly = "output/04_assemblies/megahit/{sample}.fa.gz"
	conda:
		"envs/megahit.yml"
	threads:
		64
	message:
		"[MEGAHIT] assembling {wildcards.sample}."
	shell:
		"""
		megahit -1 {input.b1} -2 {input.b2} -r {input.b3} -t {threads} -o output/04_assemblies/megahit/{wildcards.sample}/
		mv output/04_assemblies/megahit/{wildcards.sample}/final.contigs.fa output/04_assemblies/megahit/{wildcards.sample}.fa
		rm -r output/04_assemblies/megahit/{wildcards.sample}/
		gzip output/04_assemblies/megahit/{wildcards.sample}.fa
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
		"[Krona] visualizing taxonomic composition."
	shell:
		"""
		ktUpdateTaxonomy.sh
		ktImportTaxonomy -q 2 -t 3 {input.k_stdout} -o {output.krona_html}
		"""

# kraken-biom
rule kraken_biom:
	input:
		b_report = expand("output/02_taxonomic_profiling/bracken/{sample}.report", sample = SAMPLES)
	output:
		b_summary = "output/02_taxonomic_profiling/kraken_biom/bracken.biom",
		b_taxids = "output/02_taxonomic_profiling/kraken_biom/taxids.txt"
	conda:
		"envs/kraken_biom.yml"
	threads:
		1
	message:
		"[kraken-biom] converting bracken reports to biom format."
	shell:
		"""
		kraken-biom output/02_taxonomic_profiling/bracken/*.report -o {output.b_summary} --fmt json --max D --min S --otu_fp {output.b_taxids}
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
		db = config["kraken_db"],
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
		db = config["kraken_db"],
		confidence = config["kraken_confidence_score"]
	shell:
		"""
		kraken2 --db {params.db} --threads {threads} --confidence {params.confidence} --report {output.k_report} --output {output.k_stdout} --paired {input.b1} {input.b2}
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
		"[bbmap] removing host contamination in {wildcards.sample}."
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
