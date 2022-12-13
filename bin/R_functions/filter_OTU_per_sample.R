filter_OTU_per_sample<-function(phyloseq, threshold){ 
  
  ### This function removes OTUs based on a percentage read count per sample
  # returning a phyloseq object with filtered OTU table
  
  ### Arguments:
  # phyloseq: a phyloseq object, for example imported from a .biom file; includes an OTU table, taxonomy table and sample table (metadata)
  # percentage: threshold of percentage read count per sample to filter OTUs. Note: 0.05% = 0.0005!
  
  ### Function:
  # loop to calculate the OTU count per sample and filter out OTUs with count lower than percentage threshold 
  
    for(i in 1:ncol(otu_table(phyloseq)))
    {
    otu_table(phyloseq)[,i] <- sapply(otu_table(phyloseq)[,i], function(x) ifelse(x<=sum(otu_table(phyloseq)[,i])*threshold, 0, x))
    }
  
  # remove OTUs that are no longer present in the dataset
    prune_taxa(taxa_sums(phyloseq) > 0, phyloseq) 
  
  # return results
    phyloseq
  
}