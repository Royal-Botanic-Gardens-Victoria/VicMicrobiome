#!/bin/bash

# Script to process demultiplexed Illumina MiSeq 2x300 reads using the AMPtk pipeline
# https://amptk.readthedocs.io

# AMPtk needs to be previously installed using conda:
# Load the AMPtk environment with `conda activate amptk`
# Install ITS database using `amptk install -i ITS`



# Step 1.
# Pre-processing FASTQ files using vsearch (usearch available only for 32 bit)
# A mapping file is generated at this point and can be completed with the meta data associated to the project (see below under step 4).

# -i, Input folder with FASTQ files
# -o, Output file names
# -f, Forward primer; scat samples: ITS5
# -r, Reverse primer
# -l, Sequence length for trimming/padding. Default: 300
# -p, Pad reads with Ns if shorter than --trim_len (-l). Default: off [on,off]
# --min_len, Minimum lenght to keep a read. Default: 100
# --full_length, Keep only full length sequences
# --reads, Paired-end or forward reads [paired, forward]. Default: paired
# --rescue_forward, Rescue Forward Reads if PE do not merge, e.g. long amplicons. Default: on [on,off]
# --require_primer [on, off]. Default=on: off for soil dataset
# --cleanup, Remove intermediate files

amptk illumina -i ../data/fastq/soil -o Soil -f ITS1-F -r ITS4 -l 300 -p on --min_len 150 --reads paired --rescue_forward on --require_primer off --cleanup
amptk illumina -i ../data/fastq/scat_demux -o Scat -f GGAAGTAAAAGTCGTAACAAGG -r ITS2 -l 300 -p on --min_len 150 --reads paired --rescue_forward on --require_primer off --cleanup



# Move output files to new folder "illumina"

mkdir ../data/illumina
mv Soil.* ../data/illumina
mv Scat.* ../data/illumina



# Merge files from soil and scat data

cat ../data/illumina/Scat.demux.fq.gz ../data/illumina/Soil.demux.fq.gz > ../data/illumina/combined.demux.fq.gz




# Step 2.
# Clustering at 97% using UPARSE
# This will generate an OTU table (.txt), a reference sequence file (.fasta) and a log file

amptk cluster -i ../data/illumina/combined.demux.fq.gz -o cluster --uchime_ref ITS -m 10

# -i, Input folder (from step 1))
# -o, Output file names
# -p, Clustering Radius (percent). Default: 97
# -m, Minimum size to keep an OTU (singleton filter). Default: 2
# --uchime_ref, Chimera filtering [ITS, LSU, COI, 16S, custom path]
             
# Move output files to new folder "cluster"

mkdir ../data/cluster
mv cluster* ../data/cluster



# Step 3.
# Filter OTUs based on index-bleed (based on % of total reads if no mock available)

amptk filter -i ../data/cluster/cluster.otu_table.txt -f ../data/cluster/cluster.cluster.otus.fa -o filter -p 0.005 --min_reads_otu 10

# -i, OTU table (from step 2)
# -f, Reference sequence file (from step 2)
# -o, Output file names
# -p,  Index bleed between samples (percent). Default: 0.005
# --min_reads_otu, Minimum size to keep an OTU. Default: 2

# Move output files to new folder "filter"

mkdir ../data/filter
mv filter* ../data/filter



# Step 4.
# Assign taxonomy to OTUs
# At this step, the meta data associated to the project can be added using the mapping file generated in step 1

amptk taxonomy -i ../data/filter/filter.final.txt -f ../data/filter/filter.filtered.otus.fa -o taxonomy -m ../data/illumina/combined.mapping_file.txt -d ITS1 --tax_filter Fungi

# -i, Filtered OTU table (from step 3)
# -f, Filtered reference sequence file (from step 3)
# -o, Output folder name
# -m, Mapping file
# -d, Pre-installed database [ITS1, ITS2, ITS, 16S, LSU, COI]. Default: ITS2
# --method, Taxonomy method. Default: hybrid [utax, sintax, usearch, hybrid, rdp, blast]
# --add2db, Add FASTA files to DB on the fly.
# --tax_filter, Remove OTUs that do not match filter, i.e. Fungi to keep only fungi.

# Move output files to new folder "taxonomy"
# The .biom file can be directly uploaded in R for further analyses

mkdir ../data/taxonomyM
mv taxonomy* ../data/taxonomyM



