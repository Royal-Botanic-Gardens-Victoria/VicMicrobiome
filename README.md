# Soil meta-barcoding across the state of Victoria, Australia

# <span style="color:red">*** *Under development* ***</span>


# Background information

Soil samples were collected from 106 sites across the state of Victoria, at two soil depth (0–10 cm and 10–20 cm), for a total of 212 samples.


# Repository map 

In this repository you can find scripts, data and output to analyse microbial communities (fungi and bacteria) from soil samples using Illumina MiSeq. 


[`/bin`](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/bin) contains bash and R scripts to denoise Illumina MiSeq 2x300 reads, generate OTUs, assign taxonomy and trophic mode (fungi), and analyze alpha and beta diversity. 

Program used:
- [FastQC v0.11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC v1.13](https://multiqc.info/)
- [AMPtk v1.5.4](https://github.com/nextgenusfs/amptk)

[`/data`](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/data) contains the data and metadata associated to the project, apart from large fastq files (raw data).

[`/output`](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/output) contains the results and figures generated from the scripts:
- FastQC_reports: Read quality check for each locus

[`report.md`](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/blob/main/report.md) details the advances of the project, discuss the results and provide suggestions for future studies. 