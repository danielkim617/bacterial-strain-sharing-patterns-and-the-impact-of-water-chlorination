source(file.path(here::here(),"0-config.R"))
#Load R object
R21_sample_list = readRDS(file.path(here::here(), "data/R21_sample_list_no_gps.rds"))
R21_sample_list =  R21_sample_list %>% filter(sample_id != "stools_12093101_SM-NA1P5")
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

# Add if each target species was detected in each sample
R21_sample_list_detection = R21_sample_list
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(ecoli_detect = ifelse(sample_id %in% unique(ecoli.strains$sample), "Yes", "No"))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(camp_detect = ifelse(sample_id %in% unique(camp_MAG.strains$sample), "Yes", "No"))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(enterob_detect = ifelse(sample_id %in% unique(enterob.strains$sample), "Yes", "No"))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(enterco_detect = ifelse(sample_id %in% unique(enterco.strains$sample), "Yes", "No"))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(kleb_detect = ifelse(sample_id %in% unique(kleb.strains$sample), "Yes", "No"))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(staph_detect = ifelse(sample_id %in% unique(staph.strains$sample), "Yes", "No"))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(bacteroides_detect = ifelse(sample_id %in% unique(bacteroides.strains$sample), "Yes", "No"))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(bifidobacterium_detect = ifelse(sample_id %in% unique(bifidobacterium.strains$sample), "Yes", "No"))

# Extract columns ending with "_detect"
R21_sample_list_detection %>% dplyr::select(ecoli_detect, tr, area) %>% table
R21_sample_list_detection %>% dplyr::select(enterob_detect, tr, area) %>% table
R21_sample_list_detection %>% dplyr::select(enterco_detect, tr, area) %>% table
R21_sample_list_detection %>% dplyr::select(camp_detect, tr, area) %>% table
R21_sample_list_detection %>% dplyr::select(kleb_detect, tr, area) %>% table
R21_sample_list_detection %>% dplyr::select(staph_detect, tr, area) %>% table
R21_sample_list_detection %>% dplyr::select(bacteroides_detect, tr, area) %>% table
R21_sample_list_detection %>% dplyr::select(bifidobacterium_detect, tr, area) %>% table
#
prev = R21_sample_list_detection %>% dplyr::select(area, tr, camp_detect, enterob_detect, enterco_detect, kleb_detect, ecoli_detect, bacteroides_detect, bifidobacterium_detect) 
# Convert detection columns ("Yes"/"No") to binary 1/0
prev <- prev %>%
  mutate(across(3:ncol(.), ~ ifelse(. == "Yes", 1, 0)))

# Reshape to long format
prev_long <- prev %>%
  pivot_longer(cols = 3:ncol(.), names_to = "organism", values_to = "present")


# Without grouping by water status
prevalence.all <- prev_long %>%
  group_by(area, organism) %>%
  summarise(
    positive = sum(present),         # number of "Yes"
    total = n(),                     # total number of samples
    prevalence = mean(present) * 100,
    .groups = "drop"
  )

prevalence.all$organism = factor(prevalence.all$organism, levels = c("ecoli_detect", "camp_detect", "kleb_detect" , "enterco_detect", "enterob_detect", "bacteroides_detect", "bifidobacterium_detect"),
                             labels = c("Escherichia", "Campylobacter","Klebsiella", "Enterococcus", "Enterobacter", "Bacteroides", "Bifidobacterium"))

ci <- binom.confint(prevalence.all$positive, prevalence.all$total, method = "wilson")
prevalence.all$lower <- ci$lower*100
prevalence.all$upper <- ci$upper*100

# Plot all Figure 2E
prev.all.fig = ggplot(prevalence.all, aes(x = organism, y = prevalence, fill = area)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  #facet_wrap(~ area) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  scale_fill_brewer(palette = "Set2") +  # Softer, colorblind-friendly palette
  labs(
    title = "Prevalence of Organisms by Area and Treatment",
    x = NULL,
    y = "Prevalence (%)",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    legend.position = "top"
  )

#Save into a file
ggsave(file.path(here::here(),"data/figures/prev.all.fig.pdf"), dpi = 300, scale = 0.3, width = 500, height = 550, units = "mm", plot = prev.all.fig)


## Comparison of the number of strains
# Add number of strains
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(ecoli_num = ifelse(sum(ecoli.strains$sample == sample_id) > 0, sum(ecoli.strains$sample == sample_id), 0))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(camp_num = ifelse(sum(camp_MAG.strains$sample == sample_id) > 0, sum(camp_MAG.strains$sample == sample_id), 0))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(enterob_num = ifelse(sum(enterob.strains$sample == sample_id) > 0, sum(enterob.strains$sample == sample_id), 0))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(enterco_num = ifelse(sum(enterco.strains$sample == sample_id) > 0, sum(enterco.strains$sample == sample_id), 0))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(kleb_num = ifelse(sum(kleb.strains$sample == sample_id) > 0, sum(kleb.strains$sample == sample_id), 0))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(staph_num = ifelse(sum(staph.strains$sample == sample_id) > 0, sum(staph.strains$sample == sample_id), 0))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(bacteroides_num = ifelse(sum(bacteroides.strains$sample == sample_id) > 0, sum(bacteroides.strains$sample == sample_id), 0))
R21_sample_list_detection = R21_sample_list_detection %>% rowwise() %>% mutate(bifidobacterium_num = ifelse(sum(bifidobacterium.strains$sample == sample_id) > 0, sum(bifidobacterium.strains$sample == sample_id), 0))

write.csv(R21_sample_list_detection, file.path(here::here(), "data/R21_sample_list_detection.csv"), row.names = F, quote = F)

# Plot data
R21_sample_list_detection.p = R21_sample_list_detection %>% select(area, tr, Age_group, ecoli_num, camp_num, enterob_num, enterco_num, 
                                                                   kleb_num, bacteroides_num, bifidobacterium_num)

R21_sample_list_detection.p = melt(R21_sample_list_detection.p, id.vars = c("area", "Age_group", "tr"))

num_strain_comp = data.frame(matrix(ncol=6))
colnames(num_strain_comp) = c("area", "Age_group", "target", "mean_water", "mean_control", "perm_p_value")

n <- 1
for(a in c("urban", "rural")){
  for(age in unique(R21_sample_list_detection$Age_group)){
    
    if(a == "rural" & age %in% c("0 - 19 months", "Adults > 15 years")){
      next
    }
    
    for (g in c("bacteroides_num", "bifidobacterium_num", "ecoli_num", 
                "enterco_num", "enterob_num", "camp_num", "kleb_num")) {
      
      a1 <- R21_sample_list_detection %>%
        filter(area == a, Age_group == age, tr == "Water") %>%
        pull(!!sym(g))
      
      a2 <- R21_sample_list_detection %>%
        filter(area == a, Age_group == age, tr == "Control") %>%
        pull(!!sym(g))
      
      if(length(a1) == 0 | length(a2) == 0) next
      
      num_strain_comp[n, 1] <- a
      num_strain_comp[n, 2] <- age
      num_strain_comp[n, 3] <- g
      num_strain_comp[n, 4] <- mean(a1, na.rm = TRUE)
      num_strain_comp[n, 5] <- mean(a2, na.rm = TRUE)
      num_strain_comp[n, 6] <- permutation.test(a1, a2)
      n <- n + 1
    }
  }
}
# Save into a file
write.csv(num_strain_comp, file.path(here::here(), "data/stats/num_strain_comp.csv"), quote = F, row.names = F)
