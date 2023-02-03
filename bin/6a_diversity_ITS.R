### Diversity analyses of fungal communities in the State of Victoria


### 1. Set working directory to source file


### 2. Load data and extra packages

packages = c("MASS","cowplot","data.table")
lib = lapply(packages, require, character.only = TRUE)

source("5a_filter_otu_table_ITS.R")

mc1.binary
mc2.binary
mc1.rel
mc1.prc
guild_table



### 3. Filter OTUs by guild

unique(guild_table$primary_lifestyle)

list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")]
ecm.binary=prune_taxa(list, mc1.binary)
ecm.rel=prune_taxa(list, mc1.rel)

glo.binary=subset_taxa(mc1.binary, Phylum=="Glomeromycota")
glo.rel=subset_taxa(mc1.rel, Phylum=="Glomeromycota")

list=row.names(guild_table)[which(guild_table$primary_lifestyle %like% "saprotroph")]   ### include all sapotroph sub-cathegories
sap.binary=prune_taxa(list, mc1.binary)
sap.rel=prune_taxa(list, mc1.rel)

list=row.names(guild_table)[which(guild_table$primary_lifestyle=="plant_pathogen")]
pat.binary=prune_taxa(list, mc1.binary)
pat.rel=prune_taxa(list, mc1.rel)



### 4. Overall diversity

# List of 50 most frequent/abundant OTUS
input=mc1.rel   ### mc1.binary for frequency, mc1.rel for relative abundance
  
top50 <- names(sort(taxa_sums(input), TRUE)[1:50])
top50.table <- prune_taxa(top50, input)

df = as.data.frame(cbind(tax_table(top50.table),   ### create df with OTU taxonomy and counts
                         count = taxa_sums(top50.table)))

df$count = as.numeric(df$count)
df = df[order(df$count, decreasing = TRUE),]   ### reorder OTUs based on count

write.csv(df, "../output/ITS/50_most_abundant_OTUs.csv", row.names=TRUE)   ### export table (changing file name appropriately)


# List of 20 most frequent genera
mc1.genus=tax_glom(mc1.binary, taxrank = "Genus")   ### merge OTUs into genera

top20 <- names(sort(taxa_sums(mc1.genus), TRUE)[1:20])
top20.genus <- prune_taxa(top20, mc1.genus)
taxa_names(top20.genus) <- tax_table(top20.genus)[,"Genus"]   ### replace OTU names with genus name

df = as.data.frame(cbind(tax_table(top20.genus), count = taxa_sums(top20.genus)))  ### create df with genus name and OTU counts
df = dplyr::select(df, Genus, count)

df$count = as.numeric(df$count)  ### reorder genera based on count
df = df[order(df$count, decreasing = TRUE),]   

write.csv(df, "../output/ITS/20_most_frequent_genera.csv", row.names=FALSE)   ### export table


# Diversity plots 
data = ecm.binary  ### phyloseq object to use
factor = "Type"  ### variable to plot
tax = "Family"  ### taxonomic rank to use

p=plot_bar(data, x=factor, fill=tax, title = "") +
  geom_bar(stat="identity") +
  xlab("Vegetation type") +  ### change accordingly
  ylab("OTU frequency")
p

png(file="../Output/ITS/R_plots/ecm_per_vegtype.png")  ### change file name accordingly
p
dev.off()


### 5. Variables to test:

# Depth
# pH_solid_H2O
# Organic_carbon_percent
# Ammonium_nitrogen, too many missing data -> omitted
# Nitrate_nitrogen, too many missing data -> omitted
# Phosphorus_Colwell, too many missing data -> omitted

# Type (vegetation)
# State (vegetation)
# Clay_percent
# Gravel_percent
# Sand_percent
# Silt_percent, too many missing data -> omitted
# Dom_grasses_percent
# Dom_trees_percent


# Identifying missing data
data <- cbind(sample_data(mc1.rel))   ### Add metadata to alpha diversity table

which(data$pH_solid_H2O == "NA", arr.ind = TRUE) # sample 156 (402007)
which(data$Organic_carbon_percent == "NA", arr.ind = TRUE) # sample 156
which(data$Ammonium_nitrogen == "NA", arr.ind = TRUE) # 44 samples, omit in models
which(data$Nitrate_nitrogen == "NA", arr.ind = TRUE) # 90 samples, omit in models
which(data$Phosphorus_Colwell == "NA", arr.ind = TRUE) # 31 samples, omit in models
which(data$Clay_percent == "NA", arr.ind = TRUE) # samples 156 and 172 (402023)
which(data$Gravel_percent == "NA", arr.ind = TRUE) # sample 156
which(data$Sand_percent == "NA", arr.ind = TRUE) # sample 156
which(data$Silt_percent == "NA", arr.ind = TRUE) # 16 samples, omit in models

# Remove rows with missing data
rownames(data)[156]
rownames(data)[172]

subset.mc1.rel = subset_samples(mc1.rel, rownames(sample_data(mc1.rel)) !="402007")
subset.mc1.rel = subset_samples(subset.mc1.rel, rownames(sample_data(subset.mc1.rel)) !="402023")
subset.mc1.rel

subset.mc2.binary = subset_samples(mc2.binary, rownames(sample_data(mc2.binary)) !="402007")
subset.mc2.binary = subset_samples(subset.mc2.binary, rownames(sample_data(subset.mc2.binary)) !="402023")
subset.mc2.binary


### 6. Species richness

# Select data (trophic modes)
subset.mc1.rel

list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")]
ecm.subset.rel=prune_taxa(list, subset.mc1.rel)

glo.subset.rel=subset_taxa(subset.mc1.rel, Phylum=="Glomeromycota")

list=row.names(guild_table)[which(guild_table$primary_lifestyle %like% "saprotroph")]   ### include all saprotroph sub-cathegories
sap.subset.rel=prune_taxa(list, subset.mc1.rel)

list=row.names(guild_table)[which(guild_table$primary_lifestyle=="plant_pathogen")]
pat.subset.rel=prune_taxa(list, subset.mc1.rel)


# Alpha-diversity estimates
diversity <-estimate_richness(ecm.subset.rel, measures=c("Observed", "Shannon"))

# Add environmental variables
env <- c("Depth","Type","State","Dom_grasses_percent","Dom_trees_percent",
         "pH_solid_H2O","Organic_carbon_percent","Clay_percent","Gravel_percent","Sand_percent")
env_table = subset(sample_data(subset.mc1.rel), select = env)

data <- cbind(diversity, env_table)

# Convert variables to numeric or factors
data = data %>% 
  mutate_at('Depth', as.factor) %>%
  mutate_at('Type', as.factor) %>%
  mutate_at('State', as.factor) %>%
  mutate_at('Dom_grasses_percent', as.numeric) %>% 
  mutate_at('Dom_trees_percent', as.numeric) %>%
  mutate_at('pH_solid_H2O', as.numeric) %>% 
  mutate_at('Organic_carbon_percent', as.numeric) %>% 
  mutate_at('Clay_percent', as.numeric) %>%
  mutate_at('Gravel_percent', as.factor) %>%
  mutate_at('Sand_percent', as.numeric) 
data %>% str()

data$Gravel_percent  = factor(data$Gravel_percent,  ## reorder
                          levels=c("0-5", "5-10", "10-15", "15-20", "25-30","35-40","45-50","55-60"))  


# Find best predictor of species richness (generalized linear model with negative binomial distribution)
all <- glm.nb(Observed ~ Depth + 
                pH_solid_H2O +
                Organic_carbon_percent +
                Type +
                State +
                Clay_percent +
                Gravel_percent +
                Sand_percent +
                Dom_grasses_percent +
                Dom_trees_percent,
                data = data)

summary(all)
car::vif(all)  ## check for multi-colinearity of factors
dropterm(all, test = "Chisq")  ## significant predictors = pH !!    
                
# Plot
p = ggplot(data, aes(x=Dom_trees_percent, y=Observed)) + geom_point() +
  geom_smooth(method=lm) +
  ylab("Number of OTUs per sample") +
  xlab("%  Dominant trees")
p


p = ggplot(data, aes(x=Type, y=Observed)) + geom_boxplot() +
  ylab("Number of OTUs per sample") +
  xlab("Vegetation type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

png(file="../Output/ITS/R_plots/pat.richness_per_vegtype.png")  ## export graphics (changing names accordingly)
p
dev.off()


### 5. Community composition

# Select data (trophic modes)
subset.mc2.binary

list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")]
ecm.subset.binary=prune_taxa(list, subset.mc2.binary)
ecm.subset.binary = prune_samples(sample_sums(ecm.subset.binary)>0, ecm.subset.binary)  ### remove samples with zero count
ecm.subset.binary = subset_samples(sap.subset.binary, sample_names(sap.subset.binary) != '402150')   ### remove outliers

glo.subset.binary=subset_taxa(subset.mc2.binary, Phylum=="Glomeromycota")
glo.subset.binary = prune_samples(sample_sums(glo.subset.binary)>0, glo.subset.binary)

list=row.names(guild_table)[which(guild_table$primary_lifestyle %like% "saprotroph")]   ### include all saprotroph sub-categories
sap.subset.binary=prune_taxa(list, subset.mc2.binary)
sap.subset.binary = prune_samples(sample_sums(sap.subset.binary)>0, sap.subset.binary)
sap.subset.binary = subset_samples(sap.subset.binary, sample_names(sap.subset.binary) != '402042')   ### remove outliers
sap.subset.binary = subset_samples(sap.subset.binary, sample_names(sap.subset.binary) != '402043')
sap.subset.binary = subset_samples(sap.subset.binary, sample_names(sap.subset.binary) != '402077')
sap.subset.binary = subset_samples(sap.subset.binary, sample_names(sap.subset.binary) != '402103')

list=row.names(guild_table)[which(guild_table$primary_lifestyle=="plant_pathogen")]
pat.subset.binary=prune_taxa(list, subset.mc2.binary)
pat.subset.binary = prune_samples(sample_sums(pat.subset.binary)>0, pat.subset.binary)


# Ordination (Raup-Crick distance of presence/absence data)
data_selected = pat.subset.binary

ordination = ordinate(data_selected, method = "NMDS", distance = "raup", trymax=500, k=2, noshare=(engine="isoMDS"))
stressplot(ordination)
str(ordination)

# Test best environmental predictors
env <- c("Depth","Type","State","Dom_grasses_percent","Dom_trees_percent",
         "pH_solid_H2O","Organic_carbon_percent","Clay_percent","Gravel_percent","Sand_percent")

data = cbind(subset(sample_data(data_selected), select = env))

# Convert variables to numeric or factors
data = data %>% 
  mutate_at('Depth', as.factor) %>%
  mutate_at('Type', as.factor) %>%
  mutate_at('State', as.factor) %>%
  mutate_at('Dom_grasses_percent', as.numeric) %>% 
  mutate_at('Dom_trees_percent', as.numeric) %>%
  mutate_at('pH_solid_H2O', as.numeric) %>% 
  mutate_at('Organic_carbon_percent', as.numeric) %>% 
  mutate_at('Clay_percent', as.numeric) %>%
  mutate_at('Gravel_percent', as.factor) %>%
  mutate_at('Sand_percent', as.numeric) 
data %>% str()

data$Gravel_percent  = factor(data$Gravel_percent,  ## reorder
                              levels=c("0-5", "5-10", "10-15", "15-20", "25-30","35-40","45-50","55-60"))  
# Best predictors
predictors <- envfit(ordination, data, permutations = 999, na.rm = TRUE)
predictors  ## all fungi: all variables are significant apart from Depth and Gravel %
            ## ecm: same (outliers removed)
            ## glo: none
            ## sap: same (outliers removed)
            ## pat: all except pH, Clay, Depth and Gravel, but strage ordination 

# NMDS plots
NMDS = data.frame(scores(ordination)$sites, data)

ggplot(NMDS, aes(NMDS1, NMDS2, color=Type)) +   ### get names to remove outliers
  labs(color = "Vegetation type") +
  theme(aspect.ratio=1) + geom_point(size=3) +
  geom_text(label=rownames(NMDS), angle =145)

p1=ggplot(NMDS, aes(NMDS1, NMDS2, color=Type)) +
  labs(color = "Vegetation type") +
  theme(aspect.ratio=1) + geom_point(size=3)
p1
p2=ggplot(NMDS, aes(NMDS1, NMDS2, color=State)) +
  labs(color = "Vegetation state") +
  theme(aspect.ratio=1) + geom_point(size=3)
p2
p3=ggplot(NMDS, aes(NMDS1, NMDS2, color=Dom_grasses_percent)) +
  labs(color = "Dominant grasses (%)") +
  theme(aspect.ratio=1) + geom_point(size=3) +
  scale_colour_gradient(low = "yellow", high = "dark blue")
p3
p4=ggplot(NMDS, aes(NMDS1, NMDS2, color=Dom_trees_percent)) +
  labs(color = "Dominant trees (%)") +
  theme(aspect.ratio=1) + geom_point(size=3) +
  scale_colour_gradient(low = "yellow", high = "dark blue")
p4
p5=ggplot(NMDS, aes(NMDS1, NMDS2, color=pH_solid_H2O)) +
  labs(color = "pH (solid H20)") +
  theme(aspect.ratio=1) + geom_point(size=3) +
  scale_colour_gradient(low = "yellow", high = "dark blue")
p5
p6=ggplot(NMDS, aes(NMDS1, NMDS2, color=Organic_carbon_percent)) +
  labs(color = "Organic C (%)") +
  theme(aspect.ratio=1) + geom_point(size=3) +
  scale_colour_gradient(low = "yellow", high = "dark blue")
p6
p7=ggplot(NMDS, aes(NMDS1, NMDS2, color=Clay_percent)) + 
  labs(color = "Clay (%)") +
  theme(aspect.ratio=1) + geom_point(size=3) +
  scale_colour_gradient(low = "yellow", high = "dark blue")
p7
p8=ggplot(NMDS, aes(NMDS1, NMDS2, color=Sand_percent)) + 
  labs(color = "Sand (%)") +
  theme(aspect.ratio=1) + geom_point(size=3) +
  scale_colour_gradient(low = "yellow", high = "dark blue")
p8

png(file="../Output/ITS/R_plots/NMDS_plots_ecm.png", width = 750, height = 1000,)  ## export graphic
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8, ncol=2, nrow=4, align = "v")
dev.off()
