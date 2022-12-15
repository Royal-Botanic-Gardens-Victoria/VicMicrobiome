### Script to filter out contaminant reads (negative controls)
### and low abundance taxa due to index bleed from the OTU table


### 1. Set working directory to source file


### 2. Load the R packages needed

packages = c("ggplot2","phyloseq","vegan","dplyr","tidyr")
lib = lapply(packages, require, character.only = TRUE)
theme_set(theme_bw())


### 3. Import data from amptk (.biom file)

amptk = import_biom("../data/ITS/amptk/taxonomyITS.biom")
amptk

head(otu_table(amptk))  ### OTU table
head(tax_table(amptk))   ### Taxonomy table
head(sample_data(amptk))   ### Sample table (meta-data)

# Rename taxonomic rank names in taxonomy table
rank_names(amptk)
colnames(tax_table(amptk)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_names(amptk)

# Clean names in taxonomy table
tax_table(amptk)
tax_table(amptk)[, colnames(tax_table(amptk))] <- gsub(tax_table(amptk)[, colnames(tax_table(amptk))], pattern = "[a-z]__", replacement = "")
tax_table(amptk)

# Rename sample table variables (too long or to avoid weird characters)
colnames(sample_data(amptk))[c(43,65,122,148,160,183,217,221)]
colnames(sample_data(amptk))[c(43,65,122,148,160,183,217,221)] = c("Ammonium_nitrogen", 
                                                                 "Clay_percent", 
                                                                 "Gravel_percent",
                                                                 "Nitrate_nitrogen",
                                                                 "Phosphorus_Colwell",
                                                                 "Silt_percent",
                                                                 "Dom_grasses_percent",
                                                                 "Dom_trees_percent")

# Correct entries in sample data (Gravel %)
sample_data(amptk)$Gravel_percent
sample_data(amptk)$Gravel_percent[sample_data(amptk)$Gravel_percent == "0"] <- "0-5"
sample_data(amptk)$Gravel_percent[sample_data(amptk)$Gravel_percent == "5"] <- "0-5"
sample_data(amptk)$Gravel_percent[sample_data(amptk)$Gravel_percent == "10-May"] <- "5-10"
sample_data(amptk)$Gravel_percent[sample_data(amptk)$Gravel_percent == "15-Oct"] <- "10-15"
sample_data(amptk)$Gravel_percent


### 4. OTU Filtering

# Function to remove reads present in negative control samples
source("R_functions/remove_nc_reads.R")

# 2 runs with two nc samples in each run
run1 = subset_samples(amptk, Run == "T0127")
run2 = subset_samples(amptk, Run == "T0101")

run1a_nc = remove_nc_reads(phyloseq=amptk, batch_id = "T0127", NC_name = "Run1a.No")
otu_table(run1) = run1a_nc   ### replace OTU table with the filtered table

run1b_nc = remove_nc_reads(phyloseq=run1, batch_id = "T0127", NC_name = "Run1b.No")
otu_table(run1) = run1b_nc

run2a_nc = remove_nc_reads(phyloseq=amptk, batch_id = "T0101", NC_name = "Run2a.No")
otu_table(run2) = run2a_nc

run2b_nc = remove_nc_reads(phyloseq=run2, batch_id = "T0101", NC_name = "Run2b.No")
otu_table(run2) = run2b_nc

amptk_filtered_nc = merge_phyloseq(run1,run2)   ### merge the two runs again
as.data.frame(otu_table(amptk_filtered_nc))$Run1a.No   ### check that all reads in nc samples have been removed


# Function to remove OTUs with low percentage count per sample (index bleed)
source("R_functions/filter_OTU_per_sample.R")

amptk_filtered_bleed = filter_OTU_per_sample(object = amptk_filtered_nc, threshold = 0.5)   ### Filter out OTU count<0.05%


# Check results: mock samples have 8-9 OTUs instead of 10 (slightly conservative)
estimate_richness(amptk, measures = "Observed")
estimate_richness(amptk_filtered_nc, measures = "Observed")
estimate_richness(amptk_filtered_bleed, measures = "Observed")


### 5. Remove control samples (and extra samples from other studies)

amptk_filtered = subset_samples(amptk_filtered_bleed, !(Sample_type %in% c("Control","NA")))

amptk_filtered = prune_taxa(taxa_sums(amptk_filtered) > 0, amptk_filtered)   ### remove OTUs that are no longer present
any(taxa_sums(amptk_filtered) == 0)


### 6. Export filtered OTU and taxonomy table

write.csv(otu_table(amptk_filtered), "../output/ITS/OTU_table_ITS.csv", row.names=TRUE)
write.csv(tax_table(amptk_filtered), "../output/ITS/Taxonomy_table_ITS.csv", row.names=TRUE)


### 7. Data transformation

# Presence/absence matrix
binary_table=decostand(otu_table(amptk_filtered),method="pa")

mc1.binary = phyloseq(otu_table(binary_table, taxa_are_rows=TRUE), 
                       sample_data(amptk_filtered), 
                       tax_table(tax_table(amptk_filtered)))

# Filter out OTUs that are present in only in 1 sample
mc2.binary = prune_taxa(taxa_sums(mc1.binary) > 1, mc1.binary)
any(taxa_sums(mc2.binary) == 1)

# Relative abundance table (as count)
mc1.rel = transform_sample_counts(amptk_filtered, function(x) 100000 * x/sum(x))
otu_table(mc1.rel) = ceiling(otu_table(mc1.rel, "matrix"))   ### transform to next integer so it looks like read count

# Relative abundance (in %)
mc1.prc = transform_sample_counts(amptk_filtered, function(x) 100 * x/sum(x))


### FINAL OTU TABLES OBTAINED FOR DOWNSTREAM ANALYSES
mc1.binary
mc2.binary
mc1.rel
mc1.prc

