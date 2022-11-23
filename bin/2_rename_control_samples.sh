#!/bin/bash

##### Rename control samples so that they have unique names before pooling them together


### Names of control samples
# Mock: ATCC1010MOCK (ITS) or ATCC1002MOCK (16S)
# Positive: Soil
# Negative: No


### Duplicates
# For the two ITS runs, two duplicates were used for each control, for a total of 6 control samples per run
# The first 16S run didn't have duplicate control samples (3 in total)
# The second 16S run had three duplicates control samples (9 in total)


### This script rename control samples with prefixes "Run1a/Run1b/Run1c" and "Run2a/Run2b/Run2c"

# ITS (fungi)
cd ../data/ITS/fastq/bpa_10ad26d2_20221114T0127   # run1 (2 dupplicate controls)
for f in ATCC1010MOCK_ITS_KP8LM_ACAAGGAG-GACACTGA* ; do mv -- "$f" "Run1a.$f" ; done
for f in ATCC1010MOCK_ITS_KP8LM_CGTCCGAA-CTAGAGCT* ; do mv -- "$f" "Run1b.$f" ; done
for f in No_Template_Control_ITS_KP8LM_ACAAGGAG-TAGTGTAG* ; do mv -- "$f" "Run1a.$f" ; done
for f in No_Template_Control_ITS_KP8LM_CGTCCGAA-GACACTGA* ; do mv -- "$f" "Run1b.$f" ; done
for f in Soil_DNA_ITS_KP8LM_ACAAGGAG-TGCGTACG* ; do mv -- "$f" "Run1a.$f" ; done
for f in Soil_DNA_ITS_KP8LM_CGTCCGAA-GCTCTAGT* ; do mv -- "$f" "Run1b.$f" ; done   

cd ../bpa_4f999bb9_20221114T0101   # run2 (2 dupplicate controls)
for f in ATCC1010MOCK_ITS_KP69D_AATGGAGC-GCTCTAGT* ; do mv -- "$f" "Run2a.$f" ; done
for f in ATCC1010MOCK_ITS_KP69D_ACAAGGAG-GACACTGA* ; do mv -- "$f" "Run2b.$f" ; done
for f in No_Template_Control_ITS_KP69D_AATGGAGC-TGCGTACG* ; do mv -- "$f" "Run2a.$f" ; done
for f in No_Template_Control_ITS_KP69D_ACAAGGAG-TAGTGTAG* ; do mv -- "$f" "Run2b.$f" ; done
for f in Soil_DNA_ITS_KP69D_AATGGAGC-GACACTGA* ; do mv -- "$f" "Run2a.$f" ; done
for f in Soil_DNA_ITS_KP69D_ACAAGGAG-TGCGTACG* ; do mv -- "$f" "Run2b.$f" ; done 


# 16S (bacteria)
cd ../../../16S/fastq/bpa_761e3497_20221114T0152   # run1 (no duplicate controls)
for f in ATCC1002MOCK_16S_KN48L_CGAGAAGA-CTAGTATG* ; do mv -- "$f" "Run1a.$f" ; done
for f in No_Template_Control_16S_KN48L_CGAGAAGA-TCTACACT* ; do mv -- "$f" "Run1a.$f" ; done
for f in Soil_DNA_16S_KN48L_CGAGAAGA-GATAGCGT* ; do mv -- "$f" "Run1a.$f" ; done


cd ../bpa_76ac789f_20221114T0153   # run2 (3 dupplicate controls)
for f in ATCC1002MOCK_16S_KN482_CTAGCGAA-GTCTAGTG* ; do mv -- "$f" "Run2a.$f" ; done
for f in ATCC1002MOCK_16S_KN482_TACACGAT-AAGCAGCA* ; do mv -- "$f" "Run2b.$f" ; done
for f in ATCC1002MOCK_16S_KN482_GCGGCAAT-CTAGTATG* ; do mv -- "$f" "Run2c.$f" ; done
for f in No_Template_Control_16S_KN482_CTAGCGAA-GATAGCGT* ; do mv -- "$f" "Run2a.$f" ; done
for f in No_Template_Control_16S_KN482_TACACGAT-ACGACGTG* ; do mv -- "$f" "Run2b.$f" ; done
for f in No_Template_Control_16S_KN482_GCGGCAAT-TCTACACT* ; do mv -- "$f" "Run2c.$f" ; done
for f in Soil_DNA_16S_KN482_CTAGCGAA-CTAGTATG* ; do mv -- "$f" "Run2a.$f" ; done
for f in Soil_DNA_16S_KN482_TACACGAT-ACGCGTGA* ; do mv -- "$f" "Run2b.$f" ; done
for f in Soil_DNA_16S_KN482_GCGGCAAT-GATAGCGT* ; do mv -- "$f" "Run2c.$f" ; done