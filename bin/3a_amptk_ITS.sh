#!/bin/bash

##### Script to process demultiplexed Illumina MiSeq 2x300 (paired end) reads using the AMPtk pipeline
# https://amptk.readthedocs.io/en/latest/index.html


### How to install AMPtk with conda

# Create new conda environment: conda create --name amptk
# Activate amptk environment: conda activate amptk
# Install amptk: conda install -c bioconda -amptk


### Running the script

# Activate amptk environment: conda activate amptk
# Install ITS database: amptk install -i ITS --force (to update)


### Step 1. Denoise FASTQ files using vsearch (usearch available only for 32 bit)

# A mapping file is generated at this point and can be completed with the meta data associated to the project (see below under step 4).

# -i, Input folder with FASTQ files
# -o, Output file names
# -f, Forward primer name or sequence
# -r, Reverse primer name or sequence
# -l, Sequence length for trimming/padding. Default: 300
# -p, Pad reads with Ns if shorter than --trim_len (-l). Default: off
# --min_len, Minimum lenght to keep a read. Default: 100
# --reads, Paired-end or forward reads. Default: paired
# --rescue_forward, Rescue Forward Reads if PE do not merge (long amplicons). Default: on [on,off]
# --require_primer, Need to detect primer to keep a read. Default: on
# --cleanup, Remove intermediate files
# --full_length, Keep only full length sequences

amptk illumina -i ../data/ITS/fastq/bpa_10ad26d2_20221114T0127 -o T0127 -f ITS1-F -r ITS4 -l 350 -p on --min_len 150 --reads paired --rescue_forward on --require_primer off --cleanup
amptk illumina -i ../data/ITS/fastq/bpa_4f999bb9_20221114T0101 -o T0101 -f ITS1-F -r ITS4 -l 350 -p on --min_len 150 --reads paired --rescue_forward on --require_primer off --cleanup

# Move output files to new folder in /data
mkdir ../data/ITS/amptk
mv T0127.* ../data/ITS/amptk
mv T0101.* ../data/ITS/amptk

# Merge files from the two runs
cat ../data/ITS/amptk/T0127.demux.fq.gz ../data/ITS/amptk/T0101.demux.fq.gz > ../data/ITS/amptk/combined.demux.fq.gz


### Step 2. Cluster sequences at 97% similarity: This step generates an OTU table (.txt), a reference sequence file (.fasta) and a log file

# -i, Input folder (from step 1)
# -o, Output file names
# -p, Clustering radius (percent). Default: 97
# -m, Minimum size to keep an OTU (singleton filter). Default: 2
# --uchime_ref, Chimera filtering [ITS, LSU, COI, 16S, custom path]

amptk cluster -i ../data/ITS/amptk/combined.demux.fq.gz -o cluster --uchime_ref ITS
             
# Move output files to amptk folder in /data
mv cluster* ../data/ITS/amptk


### Step 3. Filter OTUs for index-bleed (based on % of total read count if no mock is available)

# -i, OTU table (from step 2)
# -f, Reference sequence file (from step 2)
# -o, Output file names
# -p,  Index bleed between samples. Default: 0.005 (0.5%)
# --min_reads_otu, Minimum size to keep an OTU. Default: 2

amptk filter -i ../data/ITS/amptk/cluster.otu_table.txt -f ../data/ITS/amptk/cluster.cluster.otus.fa -o filter -p 0.005 --min_reads_otu 10

# Move output files to amptk folder in /data
mv filter* ../data/ITS/amptk


### Step 4. Assign taxonomy to OTUs

# At this step, the mapping file generated in step 1 was modified manually:
# Samples for both run were combined
# Variables `run` (T0127, T0101) and `sample_type` (DELWP, control, NA for additional samples) were added
# All metadata provided by DELWP was added (spaces replaced by _, blank cells filled with NA)

# -i, Filtered OTU table (from step 3)
# -f, Filtered reference sequence file (from step 3)
# -o, Output folder name
# -m, Mapping file
# -d, Pre-installed database [ITS1, ITS2, ITS, 16S, LSU, COI]. Default: ITS2
# --method, Taxonomy method. Default: hybrid [utax, sintax, usearch, hybrid, rdp, blast]
# --add2db, Add FASTA files to DB on the fly.
# --tax_filter, Remove OTUs that do not match filter, i.e. Fungi to keep only fungi.

amptk taxonomy -i ../data/ITS/amptk/filter.final.txt -f ../data/ITS/amptk/filter.filtered.otus.fa -o taxonomy -m ../data/ITS/mapping_file.txt -d ITS1 --tax_filter Fungi

# Move output files to new folder "taxonomy"
mv taxonomy* ../data/ITS/amptk

# The .biom file that is generated can be uploaded in R for further diversity analyses