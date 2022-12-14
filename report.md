# Quality check of raw sequences

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


# Read quality filtering

### ITS (fungi)

The full ITS amplicon that was generated with primers ITS1F and ITS4 has an average size >600 bp, so most R1 and R2 reads could not be merged accurately (up to >90% of reads per sample). We therefore analysed the ITS1 region only, after trimming low-quality reads (Phred score <20) with [`Trimmomatic`](Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/bin/3_trimming.sh). We then used both R1 and merged reads in [`amptk`](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/bin/4a_amptk_ITS.sh). To account for the large variation in read lenght within ITS1, we set the read lenght to 350 bp and "padded" shorter reads with unknown bp (Ns).

### 16S (bacteria)
The 16S primers 27F/519R generated an amplicon of 450-490 bp per sample after merging the reads. After trimming low-quality reads (Phred score <20) with [`Trimmomatic`](Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/bin/3_trimming.sh), we set the read lenght to 500 bp (without "padding" shorter reads) and used only merged reads to improve data quality in [`amptk`](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/tree/main/bin/4b_amptk_16S.sh). We avoided the recovery of unmerged R1 reads (as for ITS above) because it generated a enormous amount of dubious sequences, including in the negative control samples.


# OTU clustering

### Note
Two approaches are commonly used to identifiy microbial diversity from sequence data: 1) amplicon sequence variants (ASV) represent unique sequences, and 2) operational taxonomic units (OTU) result from clustering sequences at a similarity threshold, usually 97%. The ASV approach has been used extensively to study bacterial communities but are expected to perform poorly for fungal groups that commonly exhibit multiple variants of rRNA gene and ITS copies per genome [(Tedersoo et al. 2022)](https://onlinelibrary.wiley.com/doi/10.1111/mec.16460). ASV approaches therefore tend to overestimate richness of common fungal species (due to haplotype variation) but underestimate richness of rare species (by removing rare variants). Because community composition is  driven by abundant taxa, the results are similar betwen ASV and OTU approaches.

### ITS (Fungi)

19,030,426 reads passed quality filtering and clustered into 17,991 OTUs. 2296 OTUs did not match Fungi, giving a total of 15,695 fungal OTUs generated with [amptk](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/blob/main/bin/4a_amptk_ITS.sh).

After filtering for index bleed and contaminants reads from negative control samples using [R](https://github.com/Royal-Botanic-Gardens-Victoria/VicMicrobiome/blob/main/bin/5a_filter_otu_table_ITS.R), as well as removing OTUs from other samples that were included in the runs, we obtained 2794 OTUs. 

The [mock community samples](https://www.atcc.org/products/msa-1010) indicated that index bleed was very high in both runs, with 78 to >200 OTUs generated instead of the expected 10. Our index bleed filter (0.5% of total read count per sample) reduced the OTU count of these samples to 8 or 9, which is slightly conservative comparatively to the 10 OTUs expected. Our approach nevertheless emphasizes that index bleed can overestimate overall diversity by more than x10. These low abundance OTUs therefore need to be filtered out based on the OTU count per sample.

Despite of the extensive sampling, 1570 OTUs (56%) were present in only one sample. It is therefore probable that our dataset only detected a portion of the fungal diversity occuring in Victorian soils.   

### 16S (Bacteria)

Only 232,938 reads passed quality filtering and clustered into 1,040 bacterial OTUs.



# Recommendations

Improve sequence quality
- In order to accurately merge ITS reads, it is important to generate ITS1 and ITS2 amplicons separately so that their lenght corresponds to the read size of the sequencing platform (300 bp for Illumina MiSeq 2x300).
- If the quality check of the run indicates that a large proportion of reads is below the 1% error rate, an extra quality filter (for example using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)) is recommended to filter raw reads before denoising/clustering. This greatly reduces the risk of mistaking sequencing error for genetic diversity.  

- Bench work needs to be optimized (equimolar concentrations between samples, spike-in?)

Reduce index bleed
- Do not allow any error when reading the barcode during demultiplexing (here two mismatches were allowed).
- Index bleed can overestimate overall diversity of a dataset by more than x10. We therefore recommand to filter low abundance OTUs based on the OTU count per sample.