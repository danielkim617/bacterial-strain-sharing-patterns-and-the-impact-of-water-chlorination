source(file.path(here::here(),"0-config.R"))
#Load R object
Kraken.phy.genus.filter <- readRDS(file = file.path(here::here(),"data/Kraken.phy.genus.filter.rds"))

#Loading straingst results
ecoli.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.ecoli.tsv"), sep = "\t")
camp.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Campylobacter.tsv"), sep = "\t")
camp_all.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Campylobacter_all.tsv"), sep = "\t")
camp_MAG.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Campylobacter_MAG.tsv"), sep = "\t")
enterob.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Enterobacter.tsv"), sep = "\t")
enterco.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Enterococcus.tsv") , sep = "\t")
kleb.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Klebsiella.tsv") , sep = "\t")
staph.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Staphylococcus.tsv") , sep = "\t")
#From Colin's DB
bacteroides.strains = read.table(file.path(here::here(), "data/colin_DB_straingst/bacteroides_straingst.tsv"), sep = "\t")
bifidobacterium.strains = read.table(file.path(here::here(), "data/colin_DB_straingst/bifidobacterium_straingst.tsv"), sep = "\t")

#Adding column names
for (i in c("ecoli.strains", "camp.strains", "camp_all.strains", "camp_MAG.strains", "enterob.strains", "enterco.strains", "kleb.strains", "staph.strains", "bacteroides.strains","bifidobacterium.strains")) {
  # Get the data frame by its name
  df <- get(i)
  # Set the new column names
  colnames(df) <- c('sample', 'i', 'strain',  'gkmers', 'ikmers', 'skmers', 'cov',  'kcov', 'gcov', 'acct', 'even', 'spec', 'rapct', 'old_rapct', 'wscore', 'score')
  # Assign the modified data frame back to its original name
  assign(i, df)
}

# Change abundances into RA
Kraken.phy.genus.RA <- transform_sample_counts(Kraken.phy.genus.filter, function(x) x / sum(x) )
tax.list = c("Escherichia", "Enterococcus", "Klebsiella", "Campylobacter", "Salmonella", "Staphylococcus", "Enterobacter", "Acinetobacter", "Pseudomonas", "Bacteroides", "Bifidobacterium")
Kraken.phy.target.RA <- subset_taxa(Kraken.phy.genus.RA, genus %in% tax.list)
#Transform data format
otu_df <- as.data.frame(otu_table(Kraken.phy.target.RA))
tax_df <- as.data.frame(tax_table(Kraken.phy.target.RA))
genus_names <- tax_df[,6]
otu_df <- cbind(otu_df, "Genus" = genus_names)
# Reshape the dataframe to a long format
Kraken.phy.genus.RA.tab <- otu_df %>%
  pivot_longer(cols = -Genus, names_to = "Sample", values_to = "Abundance")

# Prepare StrainGE data by summing up abundance by genus
straingst.ecoli = aggregate(ecoli.strains[,"rapct"], by=list(ecoli.strains$sample), FUN = sum) %>% setNames(c("Sample", "strainge_RA"))
straingst.camp = aggregate(camp.strains[,"rapct"], by=list(camp.strains$sample), FUN = sum) %>% setNames(c("Sample", "strainge_RA"))
straingst.camp.mag = aggregate(camp_MAG.strains[,"rapct"], by=list(camp_MAG.strains$sample), FUN = sum) %>% setNames(c("Sample", "strainge_RA"))
straingst.enterob = aggregate(enterob.strains[,"rapct"], by=list(enterob.strains$sample), FUN = sum) %>% setNames(c("Sample", "strainge_RA"))
straingst.enterco = aggregate(enterco.strains[,"rapct"], by=list(enterco.strains$sample), FUN = sum) %>% setNames(c("Sample", "strainge_RA"))
straingst.kleb = aggregate(kleb.strains[,"rapct"], by=list(kleb.strains$sample), FUN = sum) %>% setNames(c("Sample", "strainge_RA"))
straingst.staph = aggregate(staph.strains[,"rapct"], by=list(staph.strains$sample), FUN = sum) %>% setNames(c("Sample", "strainge_RA"))
straingst.bact = aggregate(bacteroides.strains[,"rapct"], by=list(bacteroides.strains$sample), FUN = sum) %>% setNames(c("Sample", "strainge_RA"))
straingst.bifi = aggregate(bifidobacterium.strains[,"rapct"], by=list(bifidobacterium.strains$sample), FUN = sum) %>% setNames(c("Sample", "strainge_RA"))

#Merging Kraken and StrainGE abundance
Kraken.phy.genus.RA.tab.ecoli = merge(Kraken.phy.genus.RA.tab[Kraken.phy.genus.RA.tab$Genus == "Escherichia",], straingst.ecoli, by = "Sample", all.x = T)
Kraken.phy.genus.RA.tab.camp = merge(Kraken.phy.genus.RA.tab[Kraken.phy.genus.RA.tab$Genus == "Campylobacter",], straingst.camp, by = "Sample", all.x = T)
Kraken.phy.genus.RA.tab.camp.mag = merge(Kraken.phy.genus.RA.tab[Kraken.phy.genus.RA.tab$Genus == "Campylobacter",], straingst.camp.mag, by = "Sample", all.x = T)
Kraken.phy.genus.RA.tab.camp.mag$Genus = "Campylobacter_MAG"
Kraken.phy.genus.RA.tab.enterob = merge(Kraken.phy.genus.RA.tab[Kraken.phy.genus.RA.tab$Genus == "Enterobacter",], straingst.enterob, by = "Sample", all.x = T)
Kraken.phy.genus.RA.tab.enterco = merge(Kraken.phy.genus.RA.tab[Kraken.phy.genus.RA.tab$Genus == "Enterococcus",], straingst.enterco, by = "Sample", all.x = T)
Kraken.phy.genus.RA.tab.kleb = merge(Kraken.phy.genus.RA.tab[Kraken.phy.genus.RA.tab$Genus == "Klebsiella",], straingst.kleb, by = "Sample", all.x = T)
Kraken.phy.genus.RA.tab.staph = merge(Kraken.phy.genus.RA.tab[Kraken.phy.genus.RA.tab$Genus == "Staphylococcus",], straingst.staph, by = "Sample", all.x = T)
Kraken.phy.genus.RA.tab.bact = merge(Kraken.phy.genus.RA.tab[Kraken.phy.genus.RA.tab$Genus == "Bacteroides",], straingst.bact, by = "Sample", all.x = T)
Kraken.phy.genus.RA.tab.bifi = merge(Kraken.phy.genus.RA.tab[Kraken.phy.genus.RA.tab$Genus == "Bifidobacterium",], straingst.bifi, by = "Sample", all.x = T)
#Those not detected by StrainGE
Kraken.phy.genus.RA.tab.pseudo = cbind(Kraken.phy.genus.RA.tab[Kraken.phy.genus.RA.tab$Genus == "Pseudomonas",], strainge_RA = NA)
Kraken.phy.genus.RA.tab.acine = cbind(Kraken.phy.genus.RA.tab[Kraken.phy.genus.RA.tab$Genus == "Acinetobacter",], strainge_RA = NA)
Kraken.phy.genus.RA.tab.sal = cbind(Kraken.phy.genus.RA.tab[Kraken.phy.genus.RA.tab$Genus == "Salmonella",], strainge_RA = NA)

#Combine results from all genera
Kraken.phy.genus.RA.tab = rbind(Kraken.phy.genus.RA.tab.ecoli, Kraken.phy.genus.RA.tab.camp,Kraken.phy.genus.RA.tab.camp.mag, Kraken.phy.genus.RA.tab.enterob, Kraken.phy.genus.RA.tab.enterco, Kraken.phy.genus.RA.tab.kleb, Kraken.phy.genus.RA.tab.staph, Kraken.phy.genus.RA.tab.pseudo, Kraken.phy.genus.RA.tab.acine, Kraken.phy.genus.RA.tab.sal, Kraken.phy.genus.RA.tab.bact, Kraken.phy.genus.RA.tab.bifi)
Kraken.phy.genus.RA.tab[is.na(Kraken.phy.genus.RA.tab$strainge_RA), "strainge_RA"] = 0
Kraken.phy.genus.RA.tab$Abundance = Kraken.phy.genus.RA.tab$Abundance*100 # Convert into %
#Save into rds
saveRDS(Kraken.phy.genus.RA.tab, file = file.path(here::here(), "data/Kraken.phy.genus.RA.tab.rds"))


## Enterococcus avium abundance comparison between StrainGE and Kranken2+Bracken
Kraken.phy = readRDS(file.path(here::here(), "data/Kraken.phy.rds"))

Kraken.phy.RA <- transform_sample_counts(Kraken.phy, function(x) x / sum(x) )

tax.list = c("Enterococcus avium")

Kraken.phy.avium.RA <- subset_taxa(Kraken.phy.RA, species %in% tax.list)

#Transform data format
otu_df <- as.data.frame(otu_table(Kraken.phy.avium.RA))
tax_df <- as.data.frame(tax_table(Kraken.phy.avium.RA))
sp_names <- tax_df[,7]
otu_df <- cbind(otu_df, "sp" = sp_names)

# Reshape the dataframe to a long format
Kraken.phy.avium.RA.tab <- otu_df %>%
  pivot_longer(cols = -sp, names_to = "Sample", values_to = "Abundance")
#
enterco.strains.avium = enterco.strains %>% filter(grepl("avium", strain))

straingst.enterco.avium = aggregate(enterco.strains.avium %>% select(rapct), by=list(enterco.strains.avium$sample), FUN = sum) %>% setNames(c("Sample", "strainge_RA"))

Kraken.phy.genus.RA.tab.enterco.avium = merge(Kraken.phy.avium.RA.tab[Kraken.phy.avium.RA.tab$sp == "Enterococcus avium",], straingst.enterco.avium, by = "Sample", all.x = T)
Kraken.phy.genus.RA.tab.enterco.avium$strainge_RA

Kraken.phy.genus.RA.tab.enterco.avium[is.na(Kraken.phy.genus.RA.tab.enterco.avium)] <- 0
