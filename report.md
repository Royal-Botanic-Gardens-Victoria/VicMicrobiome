# Sequence quality

The number of reads varied greatly between samples, both for ITS (fungi) and 16S (bacteria). This is likely to induce biases against rare taxa in the samples with fewer reads.

**MultiQC plot of the number of ITS forward reads per sample (fungi)**
![](output/ITS/fastqc_report/read_count_R1.png)

**MultiQC plot of the number of 16S forward reads per sample (bacteria)**
![](output/16S/fastqc_report/read_count_R1.png)

The quality of ITS reads (fungi) was below the 1% error rate (*Phred score = 20**) towards the last ca. 50 bp of the sequences. This was especially relevant for reverse (R2) reads and is likely to preclude the merging of R1 and R2 reads.


**MultiQC plot of the quality score (Phred) per base call along ITS reverse reads (fungi)**
![](output/ITS/fastqc_report/phred_score_R2.png)

The quality of 16S reads (bacteria) was below the 1% error rate for up to 2/3 of the read lenght for a proportion of the data. This is likely to induce large data loss in order to obtain robust diversity assumptions.

**MultiQC plot of the quality score (Phred) per base call along 16S reverse reads (bacteria)**
![](output/16S/fastqc_report/phred_score_R2.png)

**Phred score = 20: likelihood of finding 1 incorrect base call among 100 bases.*


# Denoising

### ITS (fungi)

The full ITS amplicon that was generated with primers ITS1F and ITS4 has an average size >600 bp, so most R1 and R2 reads could not be merged accurately (up to >90% of reads per sample). We therefore analysed the ITS1 region only, after trimming low-quality reads (Phred score <20) with [`Trimmomatic`](Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/bin/3_trimming.sh). We then used both R1 and merged reads in [`amptk`](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/bin/4a_amptk_ITS.sh). To account for the large variation in read lenght within ITS1, we set the read lenght to 350 bp and "padded" shorter reads with unknown bp (Ns).

### 16S (bacteria)
The 16S primers 27F/519R generated an amplicon of 450-490 bp per sample after merging the reads. After trimming low-quality reads (Phred score <20) with [`Trimmomatic`](Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/bin/3_trimming.sh), we set the read lenght to 500 bp (without "padding" shorter reads) and used only merged reads to improve data quality in [`amptk`](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/bin/4b_amptk_16S.sh). We avoided the recovery of unmerged R1 reads (as for ITS above) because it generated a enormous amount of dubious sequences, including in the negative control samples.


# OTU clustering

### Note
Two approaches are commonly used to identifiy microbial diversity from sequence data: 1) amplicon sequence variants (ASV) represent unique sequences, and 2) operational taxonomic units (OTU) result from clustering sequences at a similarity threshold, usually 97%. The ASV approach has been used extensively to study bacterial communities but are expected to perform poorly for fungal groups that commonly exhibit multiple variants of rRNA gene and ITS copies per genome [(Tedersoo et al. 2022)](https://onlinelibrary.wiley.com/doi/10.1111/mec.16460). ASV approaches therefore tend to overestimate richness of common fungal species (due to haplotype variation) but underestimate richness of rare species (by removing rare variants). Because community composition is  driven by abundant taxa, the results are similar betwen ASV and OTU approaches.


# Suggestions 

Improve sequence quality
- Generate amplicons with a size that corresponds to the sequence fragment (ITS1 or ITS2) so that reads can be merged acurately
- Bench work (equimolar concentrations between samples, spike-in?)

Reduce index bleed
- Do not allow any error when reading the barcode during demultiplexing (here two mismatched were allowed)