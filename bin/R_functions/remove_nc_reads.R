remove_nc_reads <- function(phyloseq, batch_id, NC_name){ 
    
  ### This functions removes the reads present in the negative control (nc) sample in all samples, for a given run or PCR batch
    # returning an OTU table (matrix) w/o the nc reads 
    
  ### Arguments:
    # phyloseq: a phyloseq object, for example imported from a .biom file; includes an OTU table, taxonomy table and sample table (metadata)
    # batch_id: character vector with the name of the batch ID. Assumes it exists as a variable of "Run" in the sample table
    # NC_name: character vector with the name of the NC sample.
    
  ### Function:
    # subset by Run
    batch_id <<- batch_id # assign batch_id to the global env
    batch <- subset_samples(phyloseq, Run == batch_id)
    print("input data")
    print(batch)
    
    # get otu matrix for each Run
    batch_OTU <- as.matrix(otu_table(batch))
    
    # make NC vector
    NC <- as.vector(batch_OTU[,NC_name])
    
    # remove reads in all samples for each OTU present in NC
    batch_clean <- batch_OTU-NC
    
    # replace negative numbers with 0 since we can't have negative read counts
    batch_clean <- replace(batch_clean, batch_clean < 0, 0)
    
    ## print the results
    print(paste("NC sample", NC_name, "was used to filter",  sum(NC), "reads in total"))
    
    ## return results
    batch_clean
    
  }