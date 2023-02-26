# Soil meta-barcoding across the state of Victoria (Australia)

# Repository map 

In this repository you can find scripts, data and output to analyse microbial communities (fungi and bacteria) from soil samples using eDNA meta-barcoding. 

[`/bin`](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/bin) contains documented bash and R scripts to quality-filter Illumina MiSeq 2x300 reads, generate OTUs, filter the OTU table, assign taxonomy and trophic mode (fungi), and analyze alpha and beta diversity patterns. 

Program used:
- [FastQC v0.11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC v1.13](https://multiqc.info/)
- [Trimmomatic v0.39](http://www.usadellab.org/cms/?page=trimmomatic)
- [AMPtk v1.5.4](https://github.com/nextgenusfs/amptk)
- [R v4.0.3](https://www.r-project.org/) 

[`/data`](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/data) contains the data and metadata associated to the project, apart from large fastq files (raw data).

[`/output`](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/output) contains the final tables and figures resulting from the scripts for each locus:
- Sequence quality check (FastQC reports)
- OTU table with raw count per sample
- Taxonomy and reference sequences of OTUs
- Tables with most frequent/abundant OTUs and fungal genera
- Diversity plots

[`technical_report.md`](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/technical_report.md) discusses each step of the pipeline used to generate the OTU tables and assign functional guilds to fungal OTUs, as well as providing suggestions for future studies.

[`diversity_report.md`](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/diversity_report.md) provides some baseline diversity analyses and briefly discuss the results.


# Background information

Soil samples were collected from 106 sites across the state of Victoria, at two soil depth (0–10 cm and 10–20 cm), for a total of 212 samples. Vegetation, edaphic variables and additional meta-data associated to the project can be found [here](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/blob/main/data/VicMicrobiome_metadata.csv).

Amplicons were generated, sequenced and demultiplexed by [BioPlatform Australia](https://bioplatforms.com/projects/australian-microbiome/) for the following regions:
- ITS (fungi) with primers ITS1F/ITS4
- 16S (bacteria) with primers 27F/519R

The following control samples were included:
- Commercial mock communities: [ATCC1010MOCK](https://www.atcc.org/products/msa-1010) for ITS and [ATCC1002MOCK](https://www.atcc.org/products/msa-1002) for 16S
- Soil = Positive control as a mixture of various soil samples
- No = Negative control

Detailed protocols for sampling and bench work can be found [here](https://www.australianmicrobiome.com/protocols/).
