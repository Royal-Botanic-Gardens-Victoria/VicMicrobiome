filter_OTU_per_sample<-function(object, threshold){ 
  
  ### This function removes OTU counts based on a percentage read count per sample
  # calculate the total OTU count per sample and filter out OTUs with count lower than percentage threshold  
  # returning a phyloseq object with a filtered OTU table
  
  ### Arguments:
  # object: a phyloseq object, for example imported from a .biom file; includes an OTU table, taxonomy table and sample table (metadata)
  # threshold: percentage read count per sample to filter OTUs.
  
  ### Function:
  
  # transform otu_table into relative abundance (percentage) 
  prc = transform_sample_counts(object, function(x) 100 * x/sum(x))
  
  # convert into dataframe and replace OTU counts below threshold with zero
  df.perc = data.frame(otu_table(prc), check.names = FALSE)
  df.perc[df.perc < threshold] <- 0
  
  # replace original values (counts) by zeros based on percentage
  df.count=data.frame(otu_table(object),check.names = FALSE)
  df.count[df.perc==0] <- 0
   
  # re-import otu table into phyloseq
  object = phyloseq(otu_table(df.count, taxa_are_rows=TRUE), 
                                  sample_data(object),
                                  tax_table(object))
  }