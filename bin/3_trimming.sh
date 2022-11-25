#!/bin/bash

##### Trimming low quality reads with Trimmomatic

### Binary file for installation: http://www.usadellab.org/cms/?page=trimmomatic

### Running the script

# Since merging the reads will be done with AMPtk (script 4), here we will trim forward (R1) and reverse (R2) reads independently
# using command trimmomatic SE (single end) with
# minimum phred score of 20 over 4 bp (sliding window) and
# minimum lenght to keep a read = 100 bp


# Create folders for results
mkdir ../data/ITS/trimmed
mkdir ../data/16S/trimmed


# Trim ITS reads (2 runs)

cd ../data/ITS/fastq/bpa_10ad26d2_20221114T0127

for f in *R1.fastq.gz;
do trimmomatic SE $f ../../trimmed/$f SLIDINGWINDOW:4:20 MINLEN:100 -summary ../../trimmed/summary_$f.txt -threads 8;
done

for f in *R2.fastq.gz;
do trimmomatic SE $f ../../trimmed/$f SLIDINGWINDOW:4:20 MINLEN:100 -summary ../../trimmed/summary$f.txt -threads 8;
done

cd ../bpa_4f999bb9_20221114T0101

for f in *R1.fastq.gz;
do trimmomatic SE $f ../../trimmed/$f SLIDINGWINDOW:4:20 MINLEN:100 -summary ../../trimmed/summary$f.txt -threads 8;
done

for f in *R2.fastq.gz;
do trimmomatic SE $f ../../trimmed/$f SLIDINGWINDOW:4:20 MINLEN:100 -summary ../../trimmed/summary$f.txt -threads 8;
done


# Trim 16S reads (2 runs)

cd ../../../16S/fastq/bpa_761e3497_20221114T0152

for f in *R1.fastq.gz;
do trimmomatic SE $f ../../trimmed/$f SLIDINGWINDOW:4:20 MINLEN:100 -summary ../../trimmed/summary$f.txt -threads 8;
done

for f in *R2.fastq.gz;
do trimmomatic SE $f ../../trimmed/$f SLIDINGWINDOW:4:20 MINLEN:100 -summary ../../trimmed/summary$f.txt -threads 8;
done

cd ../bpa_76ac789f_20221114T0153

for f in *R1.fastq.gz;
do trimmomatic SE $f ../../trimmed/$f SLIDINGWINDOW:4:20 MINLEN:100 -summary ../../trimmed/summary$f.txt -threads 8;
done

for f in *R2.fastq.gz;
do trimmomatic SE $f ../../trimmed/$f SLIDINGWINDOW:4:20 MINLEN:100 -summary ../../trimmed/summary$f.txt -threads 8;
done


# The trimmed reads in folders /data/ITS/trimmed and data/16S/trimmed are used to aliment the next step (amptk script)