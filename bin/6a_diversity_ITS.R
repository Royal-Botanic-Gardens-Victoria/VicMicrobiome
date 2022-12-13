### Diversity analyses of fungal communities in the State of Victoria

# Variables:
# Depth
# pH_solid_H2O
# Organic_carbon_percent
# Ammonium_nitrogen_mg_kg
# Nitrate_nitrogen_mg_kg
# Phosphorus_Colwell_mg_kg

# Type (vegetation)
# State (vegetation)
# Clay_percent_2_microm
# Gravel_percent_2.0_mm
# Sand_percent
# Silt_percent_2-20_microm
# Vegetation_dom._grasses_percent
# Vegetation_dom._trees_percent


### 1. Set working directory to source file


### 2. Load the data
source("5a_filter_otu_table_ITS.R")

mc1.binary
mc2.binary
mc1.rel
mc1.prc

select = subset_samples(mc1.rel, (Depth %in% "upper"))   ### selection criteria


### 3. Species diversity

# Number of OTUs per taxonomic group/treatment
plot_bar(mc1.prc, x="Type", fill="Phylum") +
  geom_bar(stat="identity")

plot_bar(mc1.binary, x="Type", fill="Phylum") +
  geom_bar(stat="identity")


# Top 20 genera
mc1.genus=tax_glom(mc1.binary, taxrank = "Genus")
top20 <- names(sort(taxa_sums(mc1.genus), TRUE)[1:20])
mc1.top20 <- prune_taxa(top20, mc1.genus)

taxa_names(mc1.top20) <- tax_table(mc1.top20)[,"Genus"]  ## change OTU names to genus
sort(taxa_sums(mc1.top20))

plot_bar(mc1.top20, fill="Genus") +
  geom_bar(stat="identity")

plot_bar(mc1.top20, x="Treatment", fill="Genus") +
  geom_bar(stat="identity")


# Alpha diversity
diversity <-estimate_richness(mc1.rel, measures=c("Observed", "Shannon"))

data <- cbind(sample_data(mc1.rel), diversity)   ### Add metadata to alpha diversity table
data$sample <- rownames(data)
data

# Plot species richness
ggplot(data, aes(x=Type, y=Observed)) + geom_boxplot() +
  ylab("Number of OTUs") +
  xlab("") +
  theme_bw()

# Anova test
anova = aov(Observed ~ Type, data = data)
summary(anova)
shapiro.test(data$Observed)   ### Normality test 



### 4. Community composition

ordination = ordinate(mc2.binary, method = "NMDS", distance = "raup", trymax=1000, k=2, noshare=(engine="isoMDS"))
stressplot(ordination)
  
# NMDS plot
plot_ordination(mc2.binary, ordination, color="Type", shape="Depth", title="") + 
    theme(aspect.ratio=1) + geom_point(size=3) 

# Adonis test
raup = distance(mc2.binary, method = "raup")   ### create distance matrix
data = cbind(sample_data(mc2.binary))

adonis2(raup ~ Type, data = data)
adonis2(raup ~ Depth, data = data)
adonis2(raup ~ Type * Depth, data = data)