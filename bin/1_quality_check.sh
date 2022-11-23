#!/bin/bash

##### Quality check with FastQC and MultiQC


### How to install FastQC and MultiQC with conda

# Create new conda environment: conda create --name fastqc
# Activate fastqc environment: conda activate fastqc
# Install fastqc: conda install -c bioconda -c conda-forge fastqc
# Install multiqc: conda install -c bioconda -c conda-forge multiqc


### Running the script

# Activate fasqc environment: conda activate fastqc

# Create folders for results
mkdir ../output/ITS/fastqc_report
mkdir ../output/16S/fastqc_report


# Analyse ITS reads (2 runs)
for f in ../data/fastqc/ITS/bpa_10ad26d2_20221114T0127/*R1.fastq.gz ../data/ITS/bpa_10ad26d2_20221114T0127/*R2.fastq.gz;
do fastqc $f --outdir=../output/ITS/fastqc_report;
done

for f in ../data/fastqc/ITS/bpa_4f999bb9_20221114T0101/*R1.fastq.gz ../data/ITS/bpa_4f999bb9_20221114T0101/*R2.fastq.gz;
do fastqc $f --outdir=../output/ITS/fastqc_report;
done


# Analyse 16S reads (2 runs)
for f in ../data/fastqc/16S/bpa_761e3497_20221114T0152/*R1.fastq.gz ../data/16S/bpa_761e3497_20221114T0152/*R2.fastq.gz;
do fastqc $f --outdir=../output/16S/fastqc_report;
done

for f in ../data/fastqc/16S/bpa_76ac789f_20221114T0153/*R1.fastq.gz ../data/16S/bpa_76ac789f_20221114T0153/*R2.fastq.gz;
do fastqc $f --outdir=../output/16S/fastqc_report;
done


# Visualize the quality analyses for each locus and direction with MultiQC
multiqc ../output/ITS/fastqc_report/*R1_fastqc.zip -n multiqc_ITS_R1 -o ../output/ITS/fastqc_report
multiqc ../output/ITS/fastqc_report/*R2_fastqc.zip -n multiqc_ITS_R2 -o ../output/ITS/fastqc_report

multiqc ../output/16S/fastqc_report/*R1_fastqc.zip -n multiqc_16S_R1 -o ../output/16S/fastqc_report
multiqc ../output/16S/fastqc_report/*R2_fastqc.zip -n multiqc_16S_R2 -o ../output/16S/fastqc_report
