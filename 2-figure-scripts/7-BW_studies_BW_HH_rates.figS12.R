source(file.path(here::here(),"0-config.R"))
#Load source data
bw.studies.test = read.csv(file.path(here::here(), "data/stats/bw.studies.test.csv"))

bw.studies.cleaned <- bw.studies.test %>%
  select(group = study1, mean_value = study1_mean, genus) %>%
  bind_rows(
    bw.studies.test %>%
      select(group = study2, mean_value = study2_mean, genus)
  ) %>% unique

bw.studies.cleaned$group = factor(bw.studies.cleaned$group, levels = c("urban", "Between studies", "rural"), labels = c("Nairobi", "Between studies", "Western Kenya"))
bw.studies.cleaned$genus = factor(bw.studies.cleaned$genus, levels = c("All", "comm","non_comm","Bact", "Bifi", "Esch", "Camp_MAG", "Kleb", "Entc", "Entb"),
                              labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))

#plot figure (Bar plot)
# Between studies
bw.studies.rates.area = ggplot(bw.studies.cleaned, aes(x = group, y = mean_value, fill = group)) +
  # Bars with fresh color palette
  geom_bar(stat = "identity", color = "black") +
  # Points with matching group color
  geom_point(aes(color = group), size = 3) +
  # Line connecting points: thicker, dotted, distinct color
  geom_line(aes(group = 1), color = "darkgrey", size = 1, linetype = "longdash") +
  # Facet by genus
  facet_wrap(~ genus, scales = "free", nrow = 2) +
  # Updated fill colors for bars
  scale_fill_manual(values = c("Nairobi" = "#1f78b4",      # blue
                               "Western Kenya" = "#33a02c", # green
                               "Between studies" = "#e31a1c")) + # red
  # Updated color for points
  scale_color_manual(values = c("Nairobi" = "#1f78b4",
                                "Western Kenya" = "#33a02c",
                                "Between studies" = "#e31a1c")) +
  labs(x = "", y = "Strain-sharing rate", fill = "", color = "", 
       title = "Between-household strain-sharing") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))


# Plot only Bifidobacterium
bw.studies.rates.area.bifi = ggplot(bw.studies.cleaned %>% filter(genus == "Bifidobacterium"), aes(x = group, y = mean_value, fill = group)) +
  # Bars with fresh color palette
  geom_bar(stat = "identity", color = "black") +
  # Points with matching group color
  geom_point(aes(color = group), size = 3) +
  # Line connecting points: thicker, dotted, distinct color
  geom_line(aes(group = 1), color = "darkgrey", size = 1, linetype = "longdash") +
  # Updated fill colors for bars
  scale_fill_manual(values = c("Nairobi" = "#1f78b4",      # blue
                               "Western Kenya" = "#33a02c", # green
                               "Between studies" = "#e31a1c")) + # red
  # Updated color for points
  scale_color_manual(values = c("Nairobi" = "#1f78b4",
                                "Western Kenya" = "#33a02c",
                                "Between studies" = "#e31a1c")) +
  labs(x = "", y = "Strain-sharing rate", fill = "", color = "", 
       title = "Bifidobacterium strain-sharing\n between-household") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

#Save into files
ggsave(file.path(here::here(),"data/figures/bw.studies.rates.area.bifi.pdf"), dpi = 300, scale = 0.3, width = 500, height = 500, units = "mm", plot = bw.studies.rates.area.bifi)

#Load source data with CI
bw.studies.test.CI = read.csv(file.path(here::here(), "data/stats/bw.studies.test.CI.csv"))

bw.studies.CI.cleaned <- bw.studies.test.CI %>%
  select(group = study1, mean_value = study1_mean, lower = study1_lower, upper = study1_upper, genus) %>%
  bind_rows(
    bw.studies.test.CI %>%
      select(group = study2, mean_value = study2_mean, lower = study2_lower, upper = study2_upper, genus)
  ) %>% unique

bw.studies.CI.cleaned$group = factor(bw.studies.CI.cleaned$group, levels = c("urban", "Between studies", "rural"), labels = c("Nairobi", "Between studies", "Western Kenya"))
bw.studies.CI.cleaned$genus = factor(bw.studies.CI.cleaned$genus, levels = c("All", "comm","non_comm","Bact", "Bifi", "Esch", "Camp_MAG", "Kleb", "Entc", "Entb"),
                              labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))


# Plot only Bifidobacterium with CI (Fig S12)
bw.studies.rates.area.bifi.CI = ggplot(bw.studies.CI.cleaned %>% filter(genus == "Bifidobacterium"), aes(x = group, y = mean_value, fill = group)) +
  # Bars with fresh color palette
  geom_bar(stat = "identity", color = "black") +
  # Points with matching group color
  geom_point(aes(color = group), size = 3) +
  # Line connecting points: thicker, dotted, distinct color
  geom_line(aes(group = 1), color = "darkgrey", size = 1, linetype = "longdash") +
  # Facet by genus
  #facet_wrap(~ genus, scales = "free", nrow = 2) +
  # Updated fill colors for bars
  scale_fill_manual(values = c("Nairobi" = "#1f78b4",      # blue
                               "Western Kenya" = "#33a02c", # green
                               "Between studies" = "#e31a1c")) + # red
  # Updated color for points
  scale_color_manual(values = c("Nairobi" = "#1f78b4",
                                "Western Kenya" = "#33a02c",
                                "Between studies" = "#e31a1c")) +
  labs(x = "", y = "Strain-sharing rate", fill = "", color = "", 
       title = "Bifidobacterium strain-sharing\n between-household") +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

#Save into files
ggsave(file.path(here::here(),"data/figures/bw.studies.rates.area.bifi.CI.pdf"), dpi = 300, scale = 0.3, width = 500, height = 500, units = "mm", plot = bw.studies.rates.area.bifi.CI)
