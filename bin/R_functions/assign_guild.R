assign_guild<-function(object, database){ 
  
  ### This function assigns trophic modes (and other traits) to OTUs based on taxonomy
  # returning a table with taxonomy and guild/traits for each OTU 
  
  ### Arguments:
  # object: a phyloseq object, for example imported from a .biom file
  # database: Reference used, for example FungalTraits (Polme et al. 2020, Fungal Diversity 105, 1-16).
  
  ### Function:
  
  # load required packages 
  require(file2meco)
  require(microeco)

  # convert phyloseq object into a "microtable"
  meco <- phyloseq2meco(object)

  # verify that OTUs and samples information is consistent across files
  meco$tidy_dataset()

  # assign guilds
  t1 <- trans_func$new(dataset = meco)
  t1$cal_spe_func(fungi_database = database)

  # create a dataframe with taxonomy and guild/traits information
  as.data.frame(t1$res_spe_func_raw_FungalTraits)


  }

