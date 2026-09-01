source(file.path(here::here(), "0-config.R"))
#Load R object
R21_sample_list = readRDS(file.path(here::here(), "data/R21_sample_list_no_gps.rds"))
R21_sample_list =  R21_sample_list %>% filter(sample_id != "stools_12093101_SM-NA1P5")

#Microbial community
#Load Bracken summary
R21_Bracken = file.path(here::here(), "data/bracken_summary.tsv") %>% read.table(., header = T, sep = "\t", fill = T)

#Make tables
R21_Bracken.count = R21_Bracken %>% dplyr::select("domain", "phylum", "class","order","family","genus", "species", "tax_id",ends_with("_num")) # nolint: line_length_linter.
R21_Bracken.count = R21_Bracken.count %>% rename_with(~ gsub("_num$", "", .x), ends_with("_num"))
R21_Bracken.count.sp = R21_Bracken.count %>% select(-c("domain", "phylum", "class","order","family","genus", "species", "tax_id"))

rownames(R21_Bracken.count.sp) = R21_Bracken.count$tax_id
colnames(R21_Bracken.count.sp) = gsub("SM.", "SM-", colnames(R21_Bracken.count.sp)) #Correcting sample names
R21_Bracken.count.sp = R21_Bracken.count.sp[,!colnames(R21_Bracken.count.sp) %in% "stools_12093101_SM-NA1P5"] # Added later

#Taxa rank (Genus/Species)
tax.info = R21_Bracken.count %>% select(c("domain", "phylum", "class","order","family","genus", "species"))
rownames(tax.info) = R21_Bracken.count$tax_id

#Make Phyloseq object
R21_sample_list.row = R21_sample_list
rownames(R21_sample_list.row) = R21_sample_list.row$sample_id
Kraken.phy <- phyloseq(otu_table(R21_Bracken.count.sp, taxa_are_rows = T), sample_data(R21_sample_list.row), tax_table(as.matrix(tax.info)))

saveRDS(Kraken.phy, file = file.path(here::here(), "data/Kraken.phy.rds"))

#Retain species annotation with RA >= 0.01% Aitor Blanco-Míguez et al., 2023
#Agglomerate by Genus
Kraken.phy.genus = tax_glom(Kraken.phy, "genus")
rel_abundances <- transform_sample_counts(Kraken.phy.genus, function(x) x / sum(x) * 100)
taxid_to_keep <- apply(rel_abundances@otu_table, 1, function(x) any(x >= 0.01))
Kraken.phy.genus.filter <- prune_taxa(taxid_to_keep, Kraken.phy.genus) # filtering out genus with RA < 0.01%

# Extract taxonomy table and convert to dataframe
tax_df <- as.data.frame(as(tax_table(Kraken.phy.genus.filter), "matrix"))
# Retain only rows where the "Domain" column is "Bacteria" or "Archaea"
filtered_taxa <- rownames(tax_df)[tax_df$domain %in% c("Bacteria", "Archaea")]
# Subset the phyloseq object to keep only those taxa
Kraken.phy.genus.filter <- prune_taxa(filtered_taxa, Kraken.phy.genus.filter)

#Save into a file
saveRDS(Kraken.phy.genus.filter, file = file.path(here::here(), "data/Kraken.phy.genus.filter.rds"))

############
#Calculate alpha diversities
Kraken.phy.genus.alpha = estimate_richness(Kraken.phy.genus.filter)
Kraken.phy.genus.alpha$sample_id = rownames(Kraken.phy.genus.alpha) %>% gsub("SM.", "SM-", .)

sample_df <- as(sample_data(Kraken.phy.genus.filter), "data.frame")

Kraken.phy.genus.alpha = merge(Kraken.phy.genus.alpha, sample_df, by= "sample_id")

#save into a file 
saveRDS(Kraken.phy.genus.alpha, file = file.path(here::here(), "data/Kraken.phy.genus.alpha.rds"))
write.csv(Kraken.phy.genus.alpha, file.path(here::here(), "data/Kraken.phy.genus.alpha.csv"), row.names = F, quote = F)

###################################################################################################
#Test the impact of treated water on Shannon diversity between treatment and control groups (by Age_group)
shannon.test.area.age.tr = data.frame(matrix(ncol=7))
colnames(shannon.test.area.age.tr) = c("area", "Age_group", "mean_w", "mean_c","median_w", "median_c", "p.value")
n=1
for (i in c("rural", "urban")) {
  for (j in c("0 - 19 months", "19 - 30 months", "> 30 months", "Adults > 15 years")) {
    k1 = Kraken.phy.genus.alpha[Kraken.phy.genus.alpha$area == i & Kraken.phy.genus.alpha$Age_group == j & Kraken.phy.genus.alpha$tr == "Water","Shannon"] %>% na.omit()
    k2 = Kraken.phy.genus.alpha[Kraken.phy.genus.alpha$area == i & Kraken.phy.genus.alpha$Age_group == j & Kraken.phy.genus.alpha$tr == "Control","Shannon"] %>% na.omit()
    if(i == "rural" & j %in% c("0 - 19 months","Adults > 15 years")){ 
      next
    }
    t = wilcox.test(k1, k2)
    shannon.test.area.age.tr[n,1] = i
    shannon.test.area.age.tr[n,2] = j
    shannon.test.area.age.tr[n,3] = mean(k1)
    shannon.test.area.age.tr[n,4] = mean(k2)
    shannon.test.area.age.tr[n,5] = median(k1)
    shannon.test.area.age.tr[n,6] = median(k2)
    shannon.test.area.age.tr[n,7] = t$p.value
    n=n+1
  }
}

#Add columns for the adjusted p values
shannon.test.area.age.tr[shannon.test.area.age.tr$area=="rural", "adj.p"] = shannon.test.area.age.tr %>% filter(area=="rural") %>% pull(p.value) %>% p.adjust(.,method = "BH")
shannon.test.area.age.tr[shannon.test.area.age.tr$area=="urban", "adj.p"] = shannon.test.area.age.tr %>% filter(area=="urban") %>% pull(p.value) %>% p.adjust(.,method = "BH")

#Write into a file
write.csv(shannon.test.area.age.tr, file.path(here::here(), "data/stats/shannon.test.area.age.tr.csv"), quote = F, row.names = F)

# Compare between age groups within each study
Kraken.phy.genus.alpha.urban = Kraken.phy.genus.alpha %>% filter(area == "urban")
Kraken.phy.genus.alpha.rural = Kraken.phy.genus.alpha %>% filter(area == "rural")

test1 = pairwise.wilcox.test(Kraken.phy.genus.alpha.urban$Shannon, Kraken.phy.genus.alpha.urban$Age_group, p.adjust.method = "BH")
test2 = pairwise.wilcox.test(Kraken.phy.genus.alpha.rural$Shannon, Kraken.phy.genus.alpha.rural$Age_group, p.adjust.method = "BH")

pvalues1 <- as.data.frame(as.table(test1$p.value)) %>% filter(!is.na(Freq))
pvalues2 <- as.data.frame(as.table(test2$p.value)) %>% filter(!is.na(Freq))

# Add a column to identify the tests
pvalues1$Test <- "urban"
pvalues2$Test <- "rural"

# Combine the results into a single dataframe
compare.shannon.age <- rbind(pvalues1, pvalues2)
colnames(compare.shannon.age) = c("Age1", "Age2", "adj.p", "area")

write.csv(compare.shannon.age, file.path(here::here(), "data/stats/compare.shannon.age.csv"), quote = F, row.names = F)

#ratio of Prevotella to Bacteroides
#Faecalibacterium / Bacteroides
otu_mat <- otu_table(Kraken.phy.genus.filter)
bact_abundance <- otu_mat[rownames(otu_mat) == "817", drop=FALSE] # 817: Bacteroides
prev_abundance <- otu_mat[rownames(otu_mat) == "165179", drop=FALSE] # 165179: Prevotella
#faecal_abundance <- otu_mat[rownames(otu_mat) == "853" , drop=FALSE] # 853: Faecalibacterium

Prevotella_Bacteroides <- t(prev_abundance / (bact_abundance + 1e-6)) %>% as.data.frame() # Avoiding zero division by adding an offset (1e-6).
colnames(Prevotella_Bacteroides) = "Prevotella/Bacteroides"
Prevotella_Bacteroides$sample_id = rownames(Prevotella_Bacteroides)


Kraken.gen.ratio = merge(Kraken.phy.genus.alpha, Prevotella_Bacteroides, by = "sample_id")

################ Figure 2D
Pre_Bac.area = ggplot(Kraken.gen.ratio[Kraken.gen.ratio$Age_group %in% c("19 - 30 months", "> 30 months"),], aes(x=Age_group, y=log2(`Prevotella/Bacteroides`), col=area)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.3) +
  scale_color_manual(values =  c("red", "black")) +
  #facet_grid(.~ area) +
  labs(x="Age group", y="log2(Prevotella / Bacteroides)", title="Prevotella ratio / Bacteroides", col="Area") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10), text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=16))

# Save into a file
ggsave(file.path(here::here(),"data/figures/Pre_Bac_area.pdf"), dpi = 300, scale = 0.3, width = 700, height = 400, units = "mm", plot = ggarrange(Pre_Bac.area, align = c("hv"), nrow=1))

#Comparison between urban and rural areas
genus.ratio.stat.area = data.frame(matrix(ncol=7))
colnames(genus.ratio.stat.area) = c("Genus_Genus","Age_group", "mean_u", "mean_r", "median_u", "median_r", "p.value")
n=1
for (i in c("Prevotella/Bacteroides")) {
  for(k in c("19 - 30 months", "> 30 months")){
    urb = Kraken.gen.ratio[Kraken.gen.ratio$area == "urban" & Kraken.gen.ratio$Age_group == k, i] %>% na.omit() %>% log2()
    rur = Kraken.gen.ratio[Kraken.gen.ratio$area == "rural" & Kraken.gen.ratio$Age_group == k, i] %>% na.omit() %>% log2()
    
    test = wilcox.test(urb, rur)
    genus.ratio.stat.area[n,1] = i
    genus.ratio.stat.area[n,2] = k
    genus.ratio.stat.area[n,3] = mean(urb)
    genus.ratio.stat.area[n,4] = mean(rur)
    genus.ratio.stat.area[n,5] = median(urb)
    genus.ratio.stat.area[n,6] = median(rur)
    genus.ratio.stat.area[n,7] = test$p.value
    n=n+1
  }
}

write.csv(genus.ratio.stat.area, file.path(here::here(),"data/stats/Pre_Bac_area.csv"), row.names = F, quote = F)

