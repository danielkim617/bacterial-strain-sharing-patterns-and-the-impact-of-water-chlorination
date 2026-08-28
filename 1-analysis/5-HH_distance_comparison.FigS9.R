source(file.path(here::here(),"0-config.R"))
#Load data
R21_sample_list = readRDS(file.path(here::here(), "data/R21_sample_list_no_gps.rds"))
R21_sample_list = R21_sample_list %>% filter(sample_id != "stools_12093101_SM-NA1P5")

HH_list = R21_sample_list %>% select(hhid, area, tr, location, gps_hhlatitude, gps_hhlongitude) %>% unique()


# Get household combinations
# rural
rural_hh_comb = HH_list %>% filter(area == "rural") %>% pull(hhid) %>% combn(., 2) %>% t() %>% as.data.frame() %>%  mutate(area = "rural")
colnames(rural_hh_comb) = c("HH1", "HH2", "area")

urban_hh_comb = HH_list %>% filter(area == "urban") %>% pull(hhid) %>% combn(., 2) %>% t() %>% as.data.frame() %>%  mutate(area = "urban")
colnames(urban_hh_comb) = c("HH1", "HH2", "area")

hh_comb = rbind(rural_hh_comb, urban_hh_comb)

for (i in 1:nrow(hh_comb)) {
  hh1_info = HH_list %>% filter(hhid == hh_comb[i,1])
  hh2_info = HH_list %>% filter(hhid == hh_comb[i,2])
  
  hh_comb[i, "hh1_location"] = hh1_info$location
  hh_comb[i, "hh2_location"] = hh2_info$location
  hh_comb[i, "hh1_tr"] = hh1_info$tr
  hh_comb[i, "hh2_tr"] = hh2_info$tr
  hh_comb[i, "distance"] = distm(hh1_info[, c("gps_hhlongitude", "gps_hhlatitude")], hh2_info[, c("gps_hhlongitude", "gps_hhlatitude")])
  
}

hh_comb = hh_comb %>% rowwise() %>% mutate(tr = ifelse(hh1_tr == hh2_tr, hh1_tr, "Different"))
hh_comb = hh_comb %>% rowwise() %>% mutate(location = ifelse(hh1_location == hh2_location, "Same", "Different"))

# Adjust location value for urban households
hh_comb = hh_comb %>% rowwise() %>% mutate(location = ifelse(hh1_tr == hh2_tr & area == "urban", "Same", location))

hh_comb$location = factor(hh_comb$location, levels = c("Same", "Different"))

# Plot data
# (Within versus Between clusters) and (within vs. between villages)
hh_dist_comp_vil = ggplot(hh_comb, aes(x = location, y = distance)) +
  geom_violin() +
  geom_boxplot(col = "blue", fill = NA, width = 0.1) +
  facet_wrap(. ~ area, scales = "free") +
  geom_hline(yintercept = 500, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Household distance location",
       x = "Location",
       y = "Household distance (m)",
  ) +
  theme_minimal(base_size = 15) + # Apply minimal theme with larger base font size
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
        axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis text for readability
  ) # Hide legend if not needed

# Save into a file
ggsave(file.path(here::here(),"data/figures/hh_dist_comp_vil.pdf"), dpi = 300, scale = 0.3, width = 500, height = 500, units = "mm", plot = hh_dist_comp_vil)


# pop up figure for rural study
# Zoomed-in version for rural
p_zoom <- 
  ggplot(hh_comb %>% filter(area == "rural"), aes(x = location, y = distance)) +
  geom_violin() +
  geom_boxplot(col = "blue", fill = NA, width = 0.1) +
  #facet_wrap(. ~ area, scales = "free") +
  geom_hline(yintercept = 500, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Household distance location",
       x = "Location",
       y = "Household distance (m)",
  ) +
  theme_minimal(base_size = 15) + # Apply minimal theme with larger base font size
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
        axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis text for readability
  )  +
  coord_cartesian(ylim = c(0, 6000)) +   # limit y-axis to 0–5000 m
  labs(title = NULL) +                   # remove title for inset
  theme(
    text = element_text(size = 10),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7),
    strip.text = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
  )

# Save into a file
ggsave(file.path(here::here(),"data/figures/p_zoom.pdf"), dpi = 300, scale = 0.3, width = 300, height = 400, units = "mm", plot = p_zoom)
