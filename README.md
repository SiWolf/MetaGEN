[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# MetaGEN
>An all-in-one solution for streamlined metagenomic analyses.

![MetaGEN Workflow](https://github.com/SiWolf/MetaGEN/blob/master/workflow.png)

MetaGEN is an all-in-one pipeline for the analysis of WGS shotgun metagenome data. MetaGEN takes raw paired-end NGS reads (.fastq.gz) as an input and will perform several analysis steps on these in order to produce a wide variety of outputs, including:
* QC reports
* Taxonomic tables
* Metagenomic assemblies (individual, co and bins)
* Resistome, virulence and plasmidome profiles
* Functional pathways and k-mers associated with each sample

MetaGEN requires [Snakemake](https://snakemake.readthedocs.io/en/stable/) to run, which will install all other dependencies automatically.

# Quick Start:
* Clone the MetaGEN repository
* Copy read files into the "input" folder
* Make sure the paired-end reads are named "<unique_name>_R1.fastq.gz" and "<unique_name>_R2.fastq.gz" respectively
* Edit the config.yml as needed
* Run MetaGEN

# Licence
This project is licensed under the GPLv3 Licence. See the [LICENSE](LICENSE) file for more details.