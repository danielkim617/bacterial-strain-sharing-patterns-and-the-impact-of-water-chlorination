source(file.path(here::here(),"0-config.R"))
#Import data
R21_sample_list = readRDS(file.path(here::here(), "data/R21_sample_list_no_gps.rds"))
Kraken.phy.genus.filter = readRDS(file = file.path(here::here(), "data/Kraken.phy.genus.filter.rds"))

##Microbial community profiles
Kraken.phy.family = tax_glom(Kraken.phy.genus.filter, "family")
Kraken.family = psmelt(Kraken.phy.family)
# Calculate the percentage
Kraken.family$Percentage <- with(Kraken.family, Abundance / tapply(Abundance, Sample, sum)[Sample])
# Identify high abundance taxa (>=1%)
high_abundance_taxa <- with(Kraken.family, family[ave(Percentage, family, FUN = max) >= 0.01])
# Replace low abundance taxa names with 'Others'
Kraken.family$family <- with(Kraken.family, replace(family, !(family %in% high_abundance_taxa), "Others"))
average_profiles <- Kraken.family %>%
  group_by(area, tr, Age_group, family) %>%
  summarize(Mean_Abundance = mean(Percentage, na.rm = TRUE))

my_colors <- brewer.pal(12, "Paired")
my_colors <- colorRampPalette(my_colors)(39)

average_profiles$family = factor(average_profiles$family, levels = c(
  "Acidaminococcaceae", "Actinomycetaceae", "Akkermansiaceae", "Atopobiaceae", "Bacteroidaceae",
  "Bifidobacteriaceae", "Campylobacteraceae", "Clostridiaceae", "Comamonadaceae", "Coprobacillaceae",
  "Coriobacteriaceae", "Desulfovibrionaceae", "Eggerthellaceae", "Enterobacteriaceae", "Enterococcaceae",
  "Erysipelotrichaceae", "Eubacteriaceae", "Hominidae", "Lachnospiraceae", "Lactobacillaceae",
  "Methanobacteriaceae", "Micrococcaceae", "Odoribacteraceae", "Oscillospiraceae",
  "Paenibacillaceae", "Pasteurellaceae", "Peptostreptococcaceae", "Prevotellaceae", "Rikenellaceae",
  "Selenomonadaceae", "Streptococcaceae", "Succinivibrionaceae", "Suoliviridae", "Sutterellaceae",
  "Tannerellaceae", "Treponemataceae", "Veillonellaceae", "Vibrionaceae","Others"))

average_profiles$area = factor(average_profiles$area, levels = c("rural", "urban"))
average_profiles$tr = factor(average_profiles$tr, levels = c("Control","Water"))

# Plot family level microbial community plot
ggplot(average_profiles[!is.na(average_profiles$Age_group),], aes(x = tr, y = Mean_Abundance, fill = family)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  labs(x = "", y = "Mean Abundance", title = "Average Microbial Community Profiles", fill="Family") +
  facet_wrap(.~ area + Age_group , scales = "free_y", nrow = 1)


##Diffrential abundance analysis (DAA)
retain_taxa_over = function(x, t = 0.20, rel_abundance_threshold = 0.0001) {
  # Convert the OTU table to relative abundances (normalize within each sample)
  # sweep() broadcasts colSums (length = n samples) across columns; plain `/`
  # recycles column-major and misaligns whenever n taxa != n samples.
  relative_abundance <- sweep(otu_table(x), 2, colSums(otu_table(x)), "/")
  
  # Mark OTUs as present in a sample if relative abundance > rel_abundance_threshold
  presence_absence <- relative_abundance > rel_abundance_threshold
  
  # Calculate prevalence as the proportion of samples (columns) where each taxon is present
  taxa_prevalence <- rowSums(presence_absence) / ncol(otu_table(x))
  
  # Retain taxa with prevalence >= threshold t
  taxa_to_keep <- taxa_prevalence >= t
  
  # Prune taxa to keep only the relevant ones
  y <- prune_taxa(taxa_to_keep, x)
  
  return(y)
}

# Genus
#Urban
Kraken.phy.genus.filter.urban = subset_samples(Kraken.phy.genus.filter, area == "urban") %>% retain_taxa_over(., t=0.05, rel_abundance_threshold = 0.001)
#Rural
Kraken.phy.genus.filter.rural = subset_samples(Kraken.phy.genus.filter, area == "rural") %>% retain_taxa_over(., t=0.05, rel_abundance_threshold = 0.001)

# to determine the threshold
compute_taxon_stats <- function(physeq, taxon_name, threshold = 0.0001) {
  # Ensure tax_table is a data frame
  tax_table_df <- as.data.frame(tax_table(physeq))
  
  # Find OTU IDs that match the given taxon name
  taxon_otu_ids <- rownames(tax_table_df)[tax_table_df$genus == taxon_name]
  
  if (length(taxon_otu_ids) == 0) {
    stop(paste("Taxon", taxon_name, "not found in the phyloseq object."))
  }
  
  # Convert OTU table to relative abundance
  otu_rel_abund <- transform_sample_counts(physeq, function(x) x / sum(x))
  
  # Extract abundance data for the selected taxon
  taxon_abundance <- as.data.frame(otu_table(otu_rel_abund))[taxon_otu_ids, , drop = FALSE]
  taxon_abundance <- colSums(taxon_abundance)  # Aggregate if multiple OTUs belong to the taxon
  
  # Calculate statistics
  prevalence <- sum(taxon_abundance > threshold) / length(taxon_abundance)  # Fraction of samples above threshold
  max_abundance <- max(taxon_abundance, na.rm = TRUE)
  min_abundance <- ifelse(any(taxon_abundance > threshold), min(taxon_abundance[taxon_abundance > threshold], na.rm = TRUE), 0)
  #mean_abundance <- mean(taxon_abundance, na.rm = TRUE)
  mean_abundance <- mean(taxon_abundance[taxon_abundance > threshold], na.rm = TRUE)
  
  # Return as a named list
  return(list(
    Taxon = taxon_name,
    Prevalence = prevalence,
    Max_Abundance = max_abundance,
    Min_Abundance = min_abundance,
    Mean_Abundance = mean_abundance,
    Threshold = threshold
  ))
}

#Running DAA 
set.seed(617)
# urban samples
urban.All.g.diff<- corncob::differentialTest(formula = ~ tr + Age_group, formula_null = ~ Age_group, phi.formula = ~ tr + Age_group,phi.formula_null = ~ tr + Age_group, data = Kraken.phy.genus.filter.urban, test = "Wald", fdr_cutoff = 0.05, boot=FALSE)
# rural samples
rural.All.g.diff = corncob::differentialTest(formula = ~ tr + Age_group, formula_null = ~ Age_group, phi.formula = ~ tr + Age_group,phi.formula_null = ~ tr + Age_group, data = Kraken.phy.genus.filter.rural, test = "Wald", fdr_cutoff = 0.05, boot=FALSE)
# function for extracting significant results
extract_significant_results <- function(differential_results, phyloseq_obj) {
  # Initialize an empty list to store results
  results <- list()
  
  # Loop through each significant model
  for (i in seq_along(differential_results$significant_models)) {
    model <- differential_results$significant_models[[i]]
    
    # Extract coefficients table
    coef_table <- model$coefficients
    
    # Extract relevant information for abundance (mu) and dispersion (phi)
    abundance <- coef_table[grep("^mu\\.", rownames(coef_table)), ]
    dispersion <- coef_table[grep("^phi\\.", rownames(coef_table)), ]
    
    # Extract TaxaID (assuming the model name matches the taxa ID)
    taxa_id <- differential_results$significant_taxa[i]
    
    # Extract taxonomy
    taxonomy <- paste(tax_table(phyloseq_obj)[taxa_id, "genus"])
    
    # Combine into a single data frame for this taxon
    results[[i]] <- data.frame(
      TaxaID = taxa_id,
      Taxonomy = taxonomy,
      Abundance_Estimate = abundance["mu.trWater", "Estimate"],
      Abundance_Std_Error = abundance["mu.trWater", "Std. Error"],
      Abundance_t_value = abundance["mu.trWater", "t value"],
      Abundance_p_value = abundance["mu.trWater", "Pr(>|t|)"],
      Dispersion_Estimate = dispersion["phi.trWater", "Estimate"],
      Dispersion_Std_Error = dispersion["phi.trWater", "Std. Error"],
      Dispersion_t_value = dispersion["phi.trWater", "t value"],
      Dispersion_p_value = dispersion["phi.trWater", "Pr(>|t|)"]
    )
  }
  
  # Combine all results into a single table
  results_table <- do.call(rbind, results)
  
  # Return the results table
  return(results_table)
}

# All
# Urban
sig.diff.urban.all = extract_significant_results(urban.All.g.diff, Kraken.phy.genus.filter.urban) %>% mutate(area = "urban")
#Rural
sig.diff.rural.all = extract_significant_results(rural.All.g.diff, Kraken.phy.genus.filter.rural) %>% mutate(area = "rural")

#Combine into one dataframe
sig.diff.all = data.frame(rbind(sig.diff.urban.all, sig.diff.rural.all))
# Add CI
sig.diff.all$CI_upper = sig.diff.all$Abundance_Estimate+1.96*sig.diff.all$Abundance_Std_Error
sig.diff.all$CI_lower = sig.diff.all$Abundance_Estimate-1.96*sig.diff.all$Abundance_Std_Error

# Save into a file
write.csv(sig.diff.all, file.path(here::here(), "data/stats/sig.diff.all1.csv"), row.names = F, quote = F)

# Check abundances of potential pathogens each
## # #  Most abundant genera Top 10## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # 
# Transform to relative abundance
Kraken.phy.genus.filter_rel <- transform_sample_counts(Kraken.phy.genus.filter, function(x) x / sum(x))

# Extract the relative abundance matrix and taxonomy table
otu_table_genus_rel <- otu_table(Kraken.phy.genus.filter_rel) # Relative abundance data
taxonomy_table <- tax_table(Kraken.phy.genus.filter_rel)      # Taxonomy data

# Calculate the mean and standard deviation of relative abundance for each genus
average_rel_abundance <- rowMeans(otu_table_genus_rel)
sd_rel_abundance <- apply(otu_table_genus_rel, 1, sd)

# Create a data frame with genera, their average relative abundances, and SD
genus_stats_df <- data.frame(
  Genus = taxonomy_table[, "genus"],          # Extract genus names
  AvgRelAbundance = average_rel_abundance,
  SDRelAbundance = sd_rel_abundance
)

# Sort by average relative abundance
genus_stats_df <- genus_stats_df[order(-genus_stats_df$AvgRelAbundance), ]

# Select the top genera (e.g., top 10)
top_genera_stats <- head(genus_stats_df, n = 5)
top_genera_stats$genus <- factor(top_genera_stats$genus, levels = top_genera_stats$genus[order(top_genera_stats$AvgRelAbundance)])

# Print the result
print(top_genera_stats)

# Create bar plot with error bars
ggplot(top_genera_stats, aes(x = genus, y = AvgRelAbundance*100)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = AvgRelAbundance*100 - SDRelAbundance*100, 
                    ymax = AvgRelAbundance*100 + SDRelAbundance*100), 
                width = 0.2, color = "black") +
  #scale_y_continuous(limits = c(0,50)) +
  coord_flip() +  # Flip coordinates for readability
  labs(title = "Relative Abundance of Genera with Standard Deviation",
       x = "Genus", y = "Average relative abundance (%)") +
  theme_minimal()

# Plot box plot of these genera
# Select taxa of interest (replace with your taxa names or IDs)
selected_taxa <- c("Prevotella", "Faecalibacterium", "Bifidobacterium","Blautia","Bacteroides")  # Example genera

# Extract OTU table and convert to data frame
otu_df <- as.data.frame(as.matrix(otu_table(Kraken.phy.genus.filter_rel)))
otu_df$OTU_ID <- rownames(otu_df)  # Add OTU IDs as a new column

# Extract taxonomy table and convert to data frame
tax_df <- as.data.frame(as.matrix(tax_table(Kraken.phy.genus.filter_rel)))
tax_df$OTU_ID <- rownames(tax_df)  # Add OTU IDs

# Merge OTU table with taxonomy to get genus names
otu_taxa <- otu_df %>%
  left_join(tax_df, by = "OTU_ID")  # Merge taxonomy info

# Filter OTUs belonging to selected genera
otu_filtered <- otu_taxa %>% filter(genus %in% selected_taxa)

# Convert to long format for ggplot
otu_long <- otu_filtered %>%
  pivot_longer(cols = -c(OTU_ID, domain, phylum,class,order,family,genus, species), names_to = "Sample", values_to = "Abundance")

otu_long$genus = factor(otu_long$genus, levels = c("Prevotella", "Faecalibacterium", "Bifidobacterium","Blautia","Bacteroides"))

# Bar plot for the number of genomes?
n.genomes = data.frame(genus = c("Prevotella", "Faecalibacterium", "Bifidobacterium", "Blautia", "Bacteroides"),
                       num = c(60, 21, 225, 17, 80))
n.genomes$genus = factor(n.genomes$genus, levels = c("Prevotella", "Faecalibacterium", "Bifidobacterium","Blautia","Bacteroides"))

top_genera_bar = ggplot(n.genomes, aes(x = genus, y = num)) +
  geom_bar(stat = "identity", fill = "blue") +  # Box plot
  theme_minimal() +
  labs(title = "Top 5 abundant genera",
       x = "Genus",
       y = "Number of complete genomes in NCBI database") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create box plot with mean overlay
top_genera = ggplot(otu_long, aes(x = genus, y = Abundance, col = genus)) +
  geom_boxplot(outlier.shape = NA) +  # Box plot
  geom_point(position = position_jitter(width = 0.2), alpha = 0.3, stroke = 0) +  # Jittered points
  #stat_summary(fun = mean, geom = "point", shape = 21, size = 2, color = "black", fill = "red") +  # Mean points
  theme_minimal() +
  labs(title = "Top 5 abundant genera",
       x = "Genus",
       y = "Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

top_genera_all = top_genera_bar / top_genera + 
  plot_layout(heights = c(1, 2))

# Save
ggsave(file.path(here::here(),"data/figures/top_genera.pdf"), dpi = 300, scale = 0.3, width = 400, height = 600, units = "mm", plot = top_genera_all)

