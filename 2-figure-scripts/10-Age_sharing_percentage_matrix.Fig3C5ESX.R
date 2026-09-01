# Figure 3C
source(file.path(here::here(),"0-config.R"))
#Load data
load(file = file.path(here::here(), "data/all.pairs.bet.age.mat.RData"))

urban.all.pairs.bet.age.mat1$genus = factor(urban.all.pairs.bet.age.mat1$genus, levels = c("All", "commensals","non_commensals","Bact", "Bifi", "Esch","Camp_MAG", "Kleb", "Entc", "Entb"), 
                              labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))
rural.all.pairs.bet.age.mat1$genus = factor(rural.all.pairs.bet.age.mat1$genus, levels = c("All", "commensals","non_commensals","Bact", "Bifi", "Esch","Camp_MAG", "Kleb", "Entc", "Entb"), 
                                       labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))
urban.pairs.bet.age.mat$genus = factor(urban.pairs.bet.age.mat$genus, levels = c("All", "commensals","non_commensals","Bact", "Bifi", "Esch","Camp_MAG", "Kleb", "Entc", "Entb"), 
                                       labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))
urban.pairs.bet.age.mat$tr = factor(urban.pairs.bet.age.mat$tr, levels = c( "Overall", "Water", "Control"))

rural.pairs.bet.age.mat$genus = factor(rural.pairs.bet.age.mat$genus, levels = c("All", "commensals","non_commensals","Bact", "Bifi", "Esch","Camp_MAG", "Kleb", "Entc", "Entb"), 
                                       labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))
rural.pairs.bet.age.mat$tr = factor(rural.pairs.bet.age.mat$tr, levels = c( "Overall", "Water", "Control"))

urban.pairs.bet.age.mat %>% filter(genus == "Commensal", tr == "Overall")

# Pairwise Fisher's exact test of strain-sharing percentage among age-group-pair categories
urban.commensal.overall = urban.pairs.bet.age.mat %>% filter(genus == "Commensal", tr == "Overall")
urban.commensal.overall.fisher = pairwise_fisher_test(urban.commensal.overall)
urban.commensal.overall.fisher$genus = "Commensal"
urban.commensal.overall.fisher$study = "Nairobi"
urban.commensal.overall.fisher %>% filter(p_adj < 0.05)

urban.noncommensal.overall = urban.pairs.bet.age.mat %>% filter(genus == "Non-commensal", tr == "Overall")
urban.noncommensal.overall.fisher = pairwise_fisher_test(urban.noncommensal.overall)
urban.noncommensal.overall.fisher$genus = "Target"
urban.noncommensal.overall.fisher$study = "Nairobi"
urban.noncommensal.overall.fisher %>% filter(p_adj < 0.05)

rural.commensal.overall = rural.pairs.bet.age.mat %>% filter(genus == "Commensal", tr == "Overall")
rural.commensal.overall.fisher = pairwise_fisher_test(rural.commensal.overall)
rural.commensal.overall.fisher$genus = "Commensal"
rural.commensal.overall.fisher$study = "Western Kenya"
rural.commensal.overall.fisher %>% filter(p_adj < 0.05)

rural.noncommensal.overall = rural.pairs.bet.age.mat %>% filter(genus == "Non-commensal", tr == "Overall")
rural.noncommensal.overall.fisher = pairwise_fisher_test(rural.noncommensal.overall)
rural.noncommensal.overall.fisher$genus = "Target"
rural.noncommensal.overall.fisher$study = "Western Kenya"
rural.noncommensal.overall.fisher %>% filter(p_adj < 0.05)

write.table(rbind(urban.commensal.overall.fisher, urban.noncommensal.overall.fisher, rural.commensal.overall.fisher, rural.noncommensal.overall.fisher), file = file.path(here::here(), "data/stats/age_pairwise_fisher_test.tsv"), sep = "\t", row.names = F, quote = F)

# For plotting a lower triangle figure of strain-sharing pair percentage
age_sharing_fig = function(x){
  # Define the factor levels
  age_levels <- c("0 - 19 months", "19 - 30 months", "> 30 months","Adults > 15 years")
  
  # Convert `Group1` and `Group2` to ordered factors
  x <- x %>%
    mutate(Group1 = factor(Group1, levels = age_levels),
           Group2 = factor(Group2, levels = age_levels))
  
  # Reposition values to the lower triangle
  lower_half_data <- x %>%
    rowwise() %>%
    mutate(Group1_lower = ifelse(as.integer(Group1) >= as.integer(Group2), as.character(Group1), as.character(Group2)),
           Group2_lower = ifelse(as.integer(Group1) >= as.integer(Group2), as.character(Group2), as.character(Group1))) %>%
    ungroup() %>%
    select(Group1 = Group1_lower, Group2 = Group2_lower, percent)
  
  # Convert `Group1` and `Group2` back to ordered factors
  lower_half_data$Group1 <- factor(lower_half_data$Group1, levels = age_levels)
  lower_half_data$Group2 <- factor(lower_half_data$Group2, levels = age_levels)
  
  # Create the heatmap with lower triangle only
  fig = ggplot(data = lower_half_data, aes(x = Group1, y = Group2, fill = percent)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         space = "Lab", name = "Strain-sharing\npairs (%)", limits = c(0, 4)) +
    theme_minimal() + 
    labs(x = "", y = "") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    coord_fixed()
  
  return(fig)
}

# Heatmap showing strain-sharing pair percentage between HHs
# urban
urban_age_pair_mat  = age_sharing_fig(urban.all.pairs.bet.age.mat1 %>% filter(genus == "All genera")) + labs(title = "Nairobi (urban)")
urban_age_pair_mat.path  = age_sharing_fig(urban.all.pairs.bet.age.mat1 %>% filter(genus == "Non-commensal")) + labs(title = "Nairobi (urban)")
urban_age_pair_mat.comm  = age_sharing_fig(urban.all.pairs.bet.age.mat1 %>% filter(genus == "Commensal")) + labs(title = "Nairobi (urban)")

# rural
rural_age_pair_mat  = age_sharing_fig(rural.all.pairs.bet.age.mat1 %>% filter(genus == "All genera")) + labs(title = "WASH B (rural)")
rural_age_pair_mat.path  = age_sharing_fig(rural.all.pairs.bet.age.mat1 %>% filter(genus == "Non-commensal")) + labs(title = "WASH B (rural)")
rural_age_pair_mat.comm  = age_sharing_fig(rural.all.pairs.bet.age.mat1 %>% filter(genus == "Commensal")) + labs(title = "WASH B (rural)")

#Save into files
ggsave(file.path(here::here(),"data/figures/urban_age_pair_mat.pdf"), dpi = 300, scale = 0.3, width = 600, height = 500, units = "mm", plot = urban_age_pair_mat)
ggsave(file.path(here::here(),"data/figures/rural_age_pair_mat.pdf"), dpi = 300, scale = 0.2, width = 600, height = 500, units = "mm", plot = rural_age_pair_mat)
ggsave(file.path(here::here(),"data/figures/All_age_pair_mat.pdf"), dpi = 300, scale = 0.3, width = 900, height = 900, units = "mm", plot = ggarrange(urban_age_pair_mat.path, urban_age_pair_mat.comm, rural_age_pair_mat.path, rural_age_pair_mat.comm))

urban.pairs.bet.age.mat$cat = factor(urban.pairs.bet.age.mat$cat, levels = c("0 - 19 months — 0 - 19 months","19 - 30 months — 19 - 30 months", "> 30 months — > 30 months", "Adults > 15 years — Adults > 15 years",
                                                                                                         "0 - 19 months — 19 - 30 months", "0 - 19 months — > 30 months", "0 - 19 months — Adults > 15 years",
                                                                                                         "19 - 30 months — > 30 months", "19 - 30 months — Adults > 15 years", "> 30 months — Adults > 15 years"),
                                                   labels = c("Younger child — Younger child", "Child — Child", "Older child — Older child", "Mother — Mother",
                                                              "Yonger child — Child", "Younger child — Older child", "Younger child — Mother",
                                                              "Child — Older child", "Child — Mother",
                                                              "Older child — Mother"))

ci <- binom.confint(urban.pairs.bet.age.mat$Yes, urban.pairs.bet.age.mat$Yes + urban.pairs.bet.age.mat$No, method = "wilson")
urban.pairs.bet.age.mat$lower <- ci$lower*100
urban.pairs.bet.age.mat$upper <- ci$upper*100

urban_age_pair_sharing_tr = ggplot(urban.pairs.bet.age.mat %>% filter(tr != "Overall" & genus %in% c("Commensal", "Non-commensal")),aes(x = cat, y = percent, fill = tr)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Control" = "#69b3a2", "Water" = "#404080")) +  # Custom colors
  theme_minimal(base_size = 14) +
  facet_wrap(.~genus) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "top"
  ) +
  labs(
    title = "Nairobi (urban)",
    x = "",
    y = "Strain-sharing pairs (%)",
    fill = "tr"
  )

# Grouping 0 - 19 months and 19- 30 months to 0 - 30 months
urban.pairs.bet.age.mat.g <- urban.pairs.bet.age.mat %>%
  mutate(
    Group1 = case_when(
      Group1 %in% c("0 - 19 months", "19 - 30 months") ~ "0 - 30 months",
      TRUE ~ Group1
    ),
    Group2 = case_when(
      Group2 %in% c("0 - 19 months", "19 - 30 months") ~ "0 - 30 months",
      TRUE ~ Group2
    )
  ) %>% select(-cat, - percent)

urban.pairs.bet.age.mat.g1 <- urban.pairs.bet.age.mat.g %>%
  group_by(Group1, Group2, genus, tr) %>%
  summarise(
    No = sum(No, na.rm = TRUE),
    Yes = sum(Yes, na.rm = TRUE),
    .groups = "drop"  # Optional: drops grouping after summarizing
  ) %>%
  mutate(
    percent = Yes / (Yes + No) * 100
  )

urban.pairs.bet.age.mat.g1$cat = paste(urban.pairs.bet.age.mat.g1$Group1, " - ", urban.pairs.bet.age.mat.g1$Group2)
urban.pairs.bet.age.mat.g1$cat = factor(urban.pairs.bet.age.mat.g1$cat, levels = c("0 - 30 months  -  0 - 30 months", "> 30 months  -  > 30 months",
                                                                                   "Adults > 15 years  -  Adults > 15 years", "0 - 30 months  -  > 30 months",
                                                                                   "0 - 30 months  -  Adults > 15 years", "> 30 months  -  Adults > 15 years"))

ci <- binom.confint(urban.pairs.bet.age.mat.g1$Yes, urban.pairs.bet.age.mat.g1$Yes + urban.pairs.bet.age.mat.g1$No, method = "wilson")
urban.pairs.bet.age.mat.g1$lower <- ci$lower*100
urban.pairs.bet.age.mat.g1$upper <- ci$upper*100

urban_age_pair_sharing_tr_g = ggplot(urban.pairs.bet.age.mat.g1 %>% filter(tr != "Overall" & genus %in% c("Commensal", "Non-commensal")),aes(x = cat, y = percent, fill = tr)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Control" = "#69b3a2", "Water" = "#404080")) +  # Custom colors
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) + 
  theme_minimal(base_size = 14) +
  scale_y_continuous(limits = c(0,4)) +
  facet_wrap(.~genus) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "top"
  ) +
  labs(
    title = "Nairobi (urban)",
    x = "",
    y = "Strain-sharing pairs (%)",
    fill = "tr"
  )

ci <- binom.confint(rural.pairs.bet.age.mat$Yes, rural.pairs.bet.age.mat$Yes + rural.pairs.bet.age.mat$No, method = "wilson")
rural.pairs.bet.age.mat$lower <- ci$lower*100
rural.pairs.bet.age.mat$upper <- ci$upper*100


rural.pairs.bet.age.mat$cat = factor(rural.pairs.bet.age.mat$cat, levels = c("0 - 19 months — 0 - 19 months","19 - 30 months — 19 - 30 months", "> 30 months — > 30 months", "Adults > 15 years — Adults > 15 years",
                                                                             "0 - 19 months — 19 - 30 months", "0 - 19 months — > 30 months", "0 - 19 months — Adults > 15 years",
                                                                             "19 - 30 months — > 30 months", "19 - 30 months — Adults > 15 years", "> 30 months — Adults > 15 years"),
                                     labels = c("Younger child — Younger child", "Child — Child", "Older child — Older child", "Mother — Mother",
                                                "Yonger child — Child", "Younger child — Older child", "Younger child — Mother",
                                                "Child — Older child", "Child — Mother",
                                                "Older child — Mother"))

rural_age_pair_sharing_tr = ggplot(rural.pairs.bet.age.mat %>% filter(tr != "Overall" & genus %in% c("Commensal", "Non-commensal")),aes(x = cat, y = percent, fill = tr)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Control" = "#69b3a2", "Water" = "#404080")) +  # Custom colors
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) + 
  theme_minimal(base_size = 14) +
  scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
  facet_wrap(.~genus) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "top"
  ) +
  labs(
    title = "WASH B (rural)",
    x = "",
    y = "Strain-sharing pairs (%)",
    fill = "tr"
  )
# Fig 5E
#Save into files
ggsave(file.path(here::here(),"data/figures/age_pair_sharing_tr.pdf"), dpi = 300, scale = 0.3, width = 830, height = 500, units = "mm", plot = ggarrange(rural_age_pair_sharing_tr, urban_age_pair_sharing_tr, ncol=2, widths = c(1:2), align = "hv"))
ggsave(file.path(here::here(),"data/figures/age_pair_sharing_tr_g.pdf"), dpi = 300, scale = 0.3, width = 850, height = 550, units = "mm", plot = ggarrange(rural_age_pair_sharing_tr, urban_age_pair_sharing_tr_g, ncol=2, widths = c(2:3), align = "hv"))
ggsave(file.path(here::here(),"data/figures/urban_age_pair_sharing_tr.pdf"), dpi = 300, scale = 0.3, width = 800, height = 600, units = "mm", plot = urban_age_pair_sharing_tr)
ggsave(file.path(here::here(),"data/figures/rural_age_pair_sharing_tr.pdf"), dpi = 300, scale = 0.3, width = 600, height = 600, units = "mm", plot = rural_age_pair_sharing_tr)

# Comparing Child - Child and Mother - Child strain-sharing in Nairobi (urban) study
# Load data
urban.all.pairs.wit.tab = read.csv(file.path(here::here(), "data/stats/urban.all.pairs.wit.tab.csv"))
urban.all.pairs.wit.tab$tr = factor(urban.all.pairs.wit.tab$tr, levels = c("Overall", "Water", "Control"))
urban.all.pairs.wit.tab$genus = factor(urban.all.pairs.wit.tab$genus, levels = c("all", "commensals", "non_commensals"), labels = c("All genera", "Commensal", "Non-commensal"))
# Add CI
ci <- binom.confint(urban.all.pairs.wit.tab$Yes, urban.all.pairs.wit.tab$Yes + urban.all.pairs.wit.tab$No, method = "wilson")
urban.all.pairs.wit.tab$lower <- ci$lower*100
urban.all.pairs.wit.tab$upper <- ci$upper*100

urban.all.pairs.wit.tab.fig = ggplot(urban.all.pairs.wit.tab, aes(x = genus, y = percent, fill = Age_group_pair)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Child - Child" = "#fdae61", "Mother - Child" = "#2b83ba")) +  # New custom colors
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  facet_wrap(.~ tr) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "top"
  ) +
  labs(
    title = "Within households in urban area",
    x = "",
    y = "Strain-sharing pairs (%)",
    fill = "Genus"
  )

# Add asterisks
urban.all.pairs.wit.tab %>% filter(tr == "Overall" & genus == "All genera") %>% select(Age_group_pair, Yes, No) 

urban.all.pairs.wit.tab.fisher.test = data.frame(matrix(ncol=4))
colnames(urban.all.pairs.wit.tab.fisher.test) = c("tr", "genus", "p_value", "label")
n=1
for (i in c("Overall", "Water", "Control")) {
  for (j in c("All genera", "Commensal", "Non-commensal")) {
    fisher_p = urban.all.pairs.wit.tab %>% filter(tr == i & genus == j) %>% select(Yes, No) %>% as.matrix.data.frame() %>% fisher.test()
    
    urban.all.pairs.wit.tab.fisher.test[n,"tr"] = i
    urban.all.pairs.wit.tab.fisher.test[n,"genus"] = j
    urban.all.pairs.wit.tab.fisher.test[n,"p_value"] = fisher_p$p.value
    if(urban.all.pairs.wit.tab.fisher.test[n,"p_value"] < 0.001){
      urban.all.pairs.wit.tab.fisher.test[n,"label"] = "***" 
    }else if(urban.all.pairs.wit.tab.fisher.test[n,"p_value"] < 0.01){
      urban.all.pairs.wit.tab.fisher.test[n,"label"] = "**" 
    }else if(urban.all.pairs.wit.tab.fisher.test[n,"p_value"] < 0.05){
      urban.all.pairs.wit.tab.fisher.test[n,"label"] = "*" 
    }else{
      urban.all.pairs.wit.tab.fisher.test[n,"label"] = "" 
    }
    n=n+1
  }
}

# None of the test results were significant for Child - Child vs. Mother - Child strain-sharing

#Save into files
ggsave(file.path(here::here(),"data/figures/urban.all.pairs.wit.tab.fig.pdf"), dpi = 300, scale = 0.3, width = 600, height = 500, units = "mm", plot = urban.all.pairs.wit.tab.fig)


