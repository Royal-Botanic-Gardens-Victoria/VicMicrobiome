#!/bin/bash

##### Script to process demultiplexed Illumina MiSeq 2x300 (paired end) reads using the AMPtk pipeline
# https://amptk.readthedocs.io/en/latest/index.html


### How to install AMPtk with conda

# Create new conda environment: conda create --name amptk
# Activate amptk environment: conda activate amptk
# Install amptk: conda install -c bioconda amptk


### Running the script

# Activate amptk environment: conda activate amptk
# Install 16S database: amptk install -i 16S --force (to update)


### Step 1. Denoise FASTQ files using vsearch (usearch available only for 32 bit)

# A mapping file is generated at this point and can be completed with the meta data associated to the project (see below under step 4).

# -i, Input folder with FASTQ files
# -o, Output file names
# -f, Forward primer name or sequence (27F)
# -r, Reverse primer name or sequence (519R)
# -l, Sequence length for trimming/padding. Default: 300
# -p, Pad reads with Ns if shorter than --trim_len (-l). Default: off
# --min_len, Minimum lenght to keep a read. Default: 100
# --reads, Paired-end or forward reads. Default: paired
# --rescue_forward, Rescue Forward Reads if PE do not merge (long amplicons). Default: on [on,off]
# --require_primer, Need to detect primer to keep a read. Default: on
# --cleanup, Remove intermediate files
# --full_length, Keep only full length sequences

amptk illumina -i ../data/16S/fastq/bpa_761e3497_20221114T0152 -o T0152 -f AGAGTTTGATCMTGGCTCAG -r GWATTACCGCGGCKGCTG -l 500 -p off --min_len 150 --reads paired --rescue_forward off --require_primer off --cleanup
amptk illumina -i ../data/16S/fastq/bpa_76ac789f_20221114T0153 -o T0153 -f AGAGTTTGATCMTGGCTCAG -r GWATTACCGCGGCKGCTG -l 500 -p off --min_len 150 --reads paired --rescue_forward off --require_primer off --cleanup

# Move output files to new folder in /data
mkdir ../data/16S/amptk
mv T0152.* ../data/16S/amptk
mv T0153.* ../data/16S/amptk

# Merge files from the two runs
cat ../data/16S/amptk/T0152.demux.fq.gz ../data/16S/amptk/T0153.demux.fq.gz > ../data/16S/amptk/combined16S.demux.fq.gz


### Step 2. Cluster sequences at 97% similarity: This step generates an OTU table (.txt), a reference sequence file (.fasta) and a log file

# -i, Input folder (from step 1)
# -o, Output file names
# -p, Clustering radius (percent). Default: 97
# -m, Minimum size to keep an OTU (singleton filter). Default: 2
# --uchime_ref, Chimera filtering [16S, LSU, COI, 16S, custom path]

amptk cluster -i ../data/16S/amptk/combined16S.demux.fq.gz -o cluster16S --uchime_ref 16S
             
# Move output files to amptk folder in /data
mv cluster16S* ../data/16S/amptk


### Step 3. Filter OTUs for index-bleed (based on % of total read count if no mock is available)

# -i, OTU table (from step 2)
# -f, Reference sequence file (from step 2)
# -o, Output file names
# -p,  Index bleed between samples (percent). Default: 0.005
# --min_reads_otu, Minimum size to keep an OTU. Default: 2

amptk filter -i ../data/16S/amptk/cluster16S.otu_table.txt -f ../data/16S/amptk/cluster16S.cluster.otus.fa -o filter16S -p 0.005 --min_reads_otu 10

# Move output files to amptk folder in /data
mv filter16S* ../data/16S/amptk


### Step 4. Assign taxonomy to OTUs

# At this step, the mapping file generated in step 1 was modified manually:
# Samples for both run were combined
# Variables `run` (T0152, T0153) and `sample_type` (DELWP, control, NA for additional samples) were added
# All metadata provided by DELWP was added (spaces replaced by _, blank cells filled with NA)

# -i, Filtered OTU table (from step 3)
# -f, Filtered reference sequence file (from step 3)
# -o, Output folder name
# -m, Mapping file
# -d, Pre-installed database Select Pre-installed database [ITS1, ITS2, ITS, 16S, LSU, COI]. Default: ITS2]. Default: ITS2
# --method, Taxonomy method. Default: hybrid [utax, sintax, usearch, hybrid, rdp, blast]
# --add2db, Add FASTA files to DB on the fly.
# --tax_filter, Remove OTUs that do not match filter, i.e. Fungi to keep only fungi.

amptk taxonomy -i ../data/16S/amptk/filter16S.final.txt -f ../data/16S/amptk/filter16S.filtered.otus.fa -o taxonomy16S -m ../data/16S/mapping_file_16S.txt -d 16S

# Move output files to new folder "taxonomy"
mv taxonomy16S* ../data/16S/amptk

# The .biom file that is generated can be uploaded in R for further diversity analyses