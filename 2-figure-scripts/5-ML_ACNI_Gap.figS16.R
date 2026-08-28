# Fig S16
source(file.path(here::here(),"0-config.R"))
#Load the output from ML strain sharing calls
ML_test=read.table(file.path(here::here(), "data/StrainGE_ML/ML_Esch_Bact_Camp_Bifi_Camp_MAG_Kleb_Entb_Entc.tsv"), header = T, sep = "\t")
R21_sample_list = readRDS(file.path(here::here(), "data/Kraken.phy.genus.alpha.rds"))
#Add Age information
ML_test$sample1_Age_group = R21_sample_list[match(ML_test$sample1, R21_sample_list$sample_id),"Age_group"]
ML_test$sample2_Age_group = R21_sample_list[match(ML_test$sample2, R21_sample_list$sample_id),"Age_group"]

#Make the class of genus name to factors
ML_test$genus = factor(ML_test$genus, levels = c("Bact", "Bifi", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb"), 
                       labels = c("Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter", "Campylobacter_MAG", "Klebsiella", "Enterococcus", "Enterobacter"))
ML_test$HH_type = factor(ML_test$HH_type, levels = c("Within", "Between"))
ML_test$combined_rf_call = factor(ML_test$combined_rf_call, levels = c(1, 0), labels = c("Same", "Different"))

ML_test = ML_test %>% filter(sample1 != "stools_12093101_SM-NA1P5" & sample2 != "stools_12093101_SM-NA1P5")
#ACNI distribution baserd on ML call (Fig. S16)
ML_ACNI_dist =  ggplot(ML_test %>% filter(genus != "Campylobacter"), aes(x = singleAgreePct, fill = as.character(combined_rf_call))) +
  geom_density(alpha = 0.5) +  # Adjust transparency with alpha
  scale_fill_brewer(palette = "Set1") +  # Better color palette
  labs(
    title = "Density of Single Agreement Percentage",
    x = "Single Agreement Percentage (ACNI)",
    fill = "RF Call",
    y = "Density"
  ) +
  scale_x_continuous(limits = c(95,100)) +
  theme_minimal() +  # Cleaner theme
  facet_wrap(.~genus) +
  theme(
    text = element_text(size = 12), # Adjust text size
    legend.position = "right"  # Adjust legend position
  )

#ACNI vs. Gap
ML_ACNI_Gap = ggplot(ML_test %>% filter(genus != "Campylobacter"), aes(x = gapJaccardSim, y = singleAgreePct, color = HH_type, size = commonPct)) +
  geom_point(alpha = 0.2) +
  facet_wrap(.~genus) +
  labs(
    x = "Gap Jaccard Similarity", 
    y = "ACNI (%)", 
    color = "HH type", 
    size = "Common Percentage (%)"
  ) +
  #geom_hline(yintercept = 99.95, linetype = "dashed", color = "blue") +  # Adding the dotted line
  theme_minimal(base_size = 12) +  # Minimal theme with base font size set for clarity
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.caption = element_text(size = 8)
  ) +
  scale_color_brewer(palette = "Set1") +  # Adjust color palette as needed
  guides(
    size = guide_legend(order = 1),
    color = guide_legend(order = 2)
  )

#Save into a file
ggsave(file.path(here::here(),"data/figures/ML_ACNI_dist.pdf"), dpi = 300, scale = 0.3, width = 590, height = 520, units = "mm", plot = ML_ACNI_dist)
ggsave(file.path(here::here(),"data/figures/ML_ACNI_Gap.pdf"), dpi = 300, scale = 0.3, width = 750, height = 580, units = "mm", plot = ML_ACNI_Gap)


# Stat test
shapiro_results <- ML_test %>%
  filter(area != "" &  genus != "Campylobacter") %>%
  group_by(genus) %>%
  sample_n(size = min(5000, n()), replace = FALSE) %>%  # Subsample if > 5000
  shapiro_test(singleAgreePct)

test <- ML_test %>%
  filter(area != "" & genus != "Campylobacter") %>%
  group_by(genus) %>%
  wilcox_test(singleAgreePct ~ combined_rf_call, p.adjust.method = "BH") %>%
  ungroup()
