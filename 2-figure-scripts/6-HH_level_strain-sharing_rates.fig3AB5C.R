source(file.path(here::here(),"0-config.R"))
# Within vs. Between strain sharing rates
#Load source data
wit.bet.test = read.csv(file.path(here::here(), "data/stats/wit.bet.test.csv"))
wit.bet.test.m = wit.bet.test[,c("area", "tr", "within", "between", "genus")] %>% melt()
wit.bet.test.m$tr = factor(wit.bet.test.m$tr, levels = c("Overall", "Water", "Control"))
wit.bet.test.m$area = factor(wit.bet.test.m$area, levels = c("rural","urban"), labels = c("Rural area", "Urban area"))
wit.bet.test.m$variable = factor(wit.bet.test.m$variable, levels = c("within", "between"), labels = c("Within households", "Between households"))
# wit.bet.test.m$genus = factor(wit.bet.test.m$genus, levels = c("All", "comm","non_comm","Bact", "Bifi", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb"), 
#                               labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))
wit.bet.test.m$genus = factor(wit.bet.test.m$genus, levels = c("All", "non_comm","comm", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb", "Bact", "Bifi"), 
                                labels = c("All genera", "Non-commensal", "Commensal", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter", "Bacteroides", "Bifidobacterium"))

wit.bet.test.m[, "split"] = "each"
wit.bet.test.m[wit.bet.test.m$genus %in% c("All genera", "Commensal", "Non-commensal"), "split"] = "combined"

#plot figure (Bar plot)
# Within vs. Between
hh.rates = ggplot(wit.bet.test.m %>% filter(tr == "Overall" & genus == "All genera"), aes(x=variable, y=value, fill=area)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  scale_fill_manual(values=c("Urban area"="darkblue", "Rural area"="orange")) +  # Custom colors for Urban and Rural
  labs(x="", y="Strain-sharing rate", fill="", title = "within vs. between households") +
  facet_wrap(~ area , scales = "free", nrow = 1) +
  scale_y_continuous(limits = c(0,0.4)) +
  #facet_grid(~ variable + split, scales = "free_x", space = "free_x") + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=12),
        plot.title = element_text(hjust=0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)) +
  guides(fill=guide_legend(title.position="top", title.hjust=0.5)) 

ggsave(file.path(here::here(),"data/figures/hh.rates.pdf"), dpi = 300, scale = 0.3, width = 450, height = 500, units = "mm", plot = hh.rates)

# Within vs. Between strain sharing rates with CI
wit.bet.test.CI = read.csv(file.path(here::here(), "data/stats/wit.bet.test.CI.csv"))

wit.bet.test.CI.m = wit.bet.test.CI %>%
  rename(within_value = within, between_value = between) %>%
  select(area, tr, genus, within_value, within_lower, within_upper, between_value, between_lower, between_upper) %>%
  pivot_longer(
    cols = -c(area, tr, genus),
    names_to = c("variable", ".value"),
    names_pattern = "(within|between)_(value|lower|upper)"
  )

wit.bet.test.CI.m$tr = factor(wit.bet.test.CI.m$tr, levels = c("Overall", "Water", "Control"))
wit.bet.test.CI.m$area = factor(wit.bet.test.CI.m$area, levels = c("rural","urban"), labels = c("Rural area", "Urban area"))
wit.bet.test.CI.m$variable = factor(wit.bet.test.CI.m$variable, levels = c("within", "between"), labels = c("Within households", "Between households"))
# wit.bet.test.CI.m$genus = factor(wit.bet.test.CI.m$genus, levels = c("All", "comm","non_comm","Bact", "Bifi", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb"), 
#                               labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))
wit.bet.test.CI.m$genus = factor(wit.bet.test.CI.m$genus, levels = c("All", "non_comm","comm", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb", "Bact", "Bifi"), 
                                labels = c("All genera", "Non-commensal", "Commensal", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter", "Bacteroides", "Bifidobacterium"))

wit.bet.test.CI.m[, "split"] = "each"
wit.bet.test.CI.m[wit.bet.test.CI.m$genus %in% c("All genera", "Commensal", "Non-commensal"), "split"] = "combined"

#plot figure (Bar plot)
# Within vs. Between (Figure 3A)
hh.rates.CI = ggplot(wit.bet.test.CI.m %>% filter(tr == "Overall" & genus == "All genera"), aes(x=variable, y=value, fill=area)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  scale_fill_manual(values=c("Urban area"="darkblue", "Rural area"="orange")) +  # Custom colors for Urban and Rural
  labs(x="", y="Strain-sharing rate", fill="", title = "within vs. between households") +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.7),
    width = 0.2
  )+
  facet_wrap(~ area , scales = "free", nrow = 1) +
  scale_y_continuous(limits = c(0,0.45)) +
  #facet_grid(~ variable + split, scales = "free_x", space = "free_x") + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=12),
        plot.title = element_text(hjust=0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)) +
  guides(fill=guide_legend(title.position="top", title.hjust=0.5)) 

ggsave(file.path(here::here(),"data/figures/hh.rates.CI.pdf"), dpi = 300, scale = 0.3, width = 450, height = 500, units = "mm", plot = hh.rates.CI)

# Urban vs. Rural
# V1
hh.rates.area = ggplot(wit.bet.test.m %>% filter(tr == "Overall" & genus != "Campylobacter_dep"), aes(x=genus, y=value, fill=area)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  scale_fill_manual(values=c("Urban area"="darkblue", "Rural area"="orange")) +  # Custom colors for Urban and Rural
  labs(x="Genera", y="Strain-sharing rate", fill="", title = "Urban vs. Rural areas") +
  facet_wrap(~ variable , scales = "free", nrow = 1) +
  #facet_grid(~ variable + split, scales = "free_x", space = "free_x") + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=12),
        plot.title = element_text(hjust=0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)) +
  guides(fill=guide_legend(title.position="top", title.hjust=0.5)) 

# Urban vs. Rural with CI
hh.rates.area.CI = ggplot(wit.bet.test.CI.m %>% filter(tr == "Overall" & genus != "Campylobacter_dep"), aes(x=genus, y=value, fill=area)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  scale_fill_manual(values=c("Urban area"="darkblue", "Rural area"="orange")) +  # Custom colors for Urban and Rural
  labs(x="Genera", y="Strain-sharing rate", fill="", title = "Urban vs. Rural areas") +
  facet_wrap(~ variable , scales = "free", nrow = 1) +
  #facet_grid(~ variable + split, scales = "free_x", space = "free_x") + 
   geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.7),
    width = 0.2
  )+
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=12),
        plot.title = element_text(hjust=0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)) +
  guides(fill=guide_legend(title.position="top", title.hjust=0.5)) 

# For asterisks
rural.urban.test = read.csv(file.path(here::here(), "data/stats/rural.urban.test.csv")) %>% filter(tr == "Overall")
rural.urban.test$variable = factor(rural.urban.test$HH_type, levels = c("Within", "Between"), labels = c("Within households", "Between households"))
#rural.urban.test$genus = factor(rural.urban.test$genus, levels = c("All", "comm","non_comm","Bact", "Bifi", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb"), 
#                                labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))
rural.urban.test$genus = factor(rural.urban.test$genus, levels = c("All", "non_comm","comm", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb", "Bact", "Bifi"), 
                                                                labels = c("All genera", "Non-commensal", "Commensal", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter", "Bacteroides", "Bifidobacterium"))
                                
rural.urban.test = rural.urban.test %>% mutate(label = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*", "")))) 

# version with CI (Figure 3B)
hh.rates.area.CI = hh.rates.area.CI +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_text(data = rural.urban.test %>% filter(tr == "Overall" & genus != "Campylobacter_dep"),
            aes(x = genus, y = pmax(rural, urban), label = label),
            vjust = -0.1, col = "red",
            size = 6,
            inherit.aes = FALSE) + 
  facet_wrap(~variable, scales = "free")

ggsave(file.path(here::here(),"data/figures/hh.rates.area.pdf"), dpi = 300, scale = 0.3, width = 700, height = 400, units = "mm", plot = hh.rates.area)
ggsave(file.path(here::here(),"data/figures/hh.rates.area.CI.pdf"), dpi = 300, scale = 0.3, width = 700, height = 400, units = "mm", plot = hh.rates.area.CI)

# Rural
hh.rates.rural = ggplot(wit.bet.test.m %>% filter(area == "Rural area" & tr != "Overall" & genus != "Campylobacter_dep"), aes(x=genus, y=value, fill=tr)) +
  geom_bar(stat="identity", position="dodge", color="black") +  # Adding outline to bars for clarity
  scale_fill_brewer(palette="Paired") +  # Using a more appealing color palette
  labs(x="", y="Strain-sharing rate", fill="Treatment", title = "Rural area") +  # Refining axis labels
  facet_wrap(~variable, scales = "free") +
  theme_minimal(base_size = 14) +  # Using a minimal theme as a base, adjust text size as needed
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),  # Adjusting x-axis text for better readability
        axis.text.y = element_text(size=12),  # Adjusting y-axis text size
        plot.title = element_text(hjust=0.5),  # Centering the plot title
        legend.title = element_text(size=12),  # Adjusting legend title size
        legend.text = element_text(size=10)) +  # Adjusting legend text size
  guides(fill=guide_legend(title.position="top", title.hjust=0.5))  # Adjusting legend position and alignment

# For asterisks
wat.cont.test_rural = read.csv(file.path(here::here(), "data/stats/wat.cont.test.csv")) %>% filter(area == "rural")
wat.cont.test_rural$variable = factor(wat.cont.test_rural$HH_type, levels = c("Within", "Between"), labels = c("Within households", "Between households"))
wat.cont.test_rural$genus = factor(wat.cont.test_rural$genus, levels = c("All", "comm","non_comm","Bact", "Bifi", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb"), 
                                labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))
wat.cont.test_rural = wat.cont.test_rural %>% mutate(label = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*", "")))) 

hh.rates.rural = hh.rates.rural +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_text(data = wat.cont.test_rural %>% filter(genus != "Campylobacter_dep"),
            aes(x = genus, y = pmax(water, control), label = label),
            vjust = -0.1, col = "red",
            size = 6,
            inherit.aes = FALSE) +
  facet_wrap(~variable, scales = "free")

#Urban
hh.rates.urban = ggplot(wit.bet.test.m %>% filter(area == "Urban area" & tr != "Overall" & genus != "Campylobacter_dep"), aes(x=genus, y=value, fill=tr)) +
  geom_bar(stat="identity", position="dodge", color="black") +  # Adding outline to bars for clarity
  scale_fill_brewer(palette="Paired") +  # Using a more appealing color palette
  labs(x="", y="Strain-sharing rate", fill="Treatment", title = "Urban area") +  # Refining axis labels
  facet_wrap(~variable, scales = "free") +
  theme_minimal(base_size = 14) +  # Using a minimal theme as a base, adjust text size as needed
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),  # Adjusting x-axis text for better readability
        axis.text.y = element_text(size=12),  # Adjusting y-axis text size
        plot.title = element_text(hjust=0.5),  # Centering the plot title
        legend.title = element_text(size=12),  # Adjusting legend title size
        legend.text = element_text(size=10)) +  # Adjusting legend text size
  guides(fill=guide_legend(title.position="top", title.hjust=0.5))  # Adjusting legend position and alignment

# For asterisks
wat.cont.test_urban = read.csv(file.path(here::here(), "data/stats/wat.cont.test.csv")) %>% filter(area == "urban")
wat.cont.test_urban$variable = factor(wat.cont.test_urban$HH_type, levels = c("Within", "Between"), labels = c("Within households", "Between households"))
wat.cont.test_urban$genus = factor(wat.cont.test_urban$genus, levels = c("All", "comm","non_comm","Bact", "Bifi", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb"), 
                                   labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))
wat.cont.test_urban = wat.cont.test_urban %>% mutate(label = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*", "")))) 
hh.rates.urban = hh.rates.urban +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_text(data = wat.cont.test_urban %>% filter(genus != "Campylobacter_dep"),
            aes(x = genus, y = pmax(water, control), label = label),
            vjust = -0.1, col = "red",
            size = 6,
            inherit.aes = FALSE) +
  facet_wrap(~variable, scales = "free")

# Save into files
ggsave(file.path(here::here(),"data/figures/hh.rates.tr.pdf"), dpi = 300, scale = 0.3, width = 700, height = 750, units = "mm", plot = ggarrange(hh.rates.urban, hh.rates.rural, nrow = 2))

#### Dumbbell style plots ####################
# Prepare data in wide format
wit.bet.test.wide.rural <- wit.bet.test.m %>%
  filter(area == "Rural area", tr != "Overall", genus != "Campylobacter_dep") %>%
  select(genus, variable, tr, value) %>%
  pivot_wider(names_from = tr, values_from = value)

# Create dumbbell plot
hh.rates.rural.dub <- ggplot(wit.bet.test.wide.rural, aes(x = genus)) +
  geom_segment(aes(y = `Control`, yend = `Water`, xend = genus),
               color = "gray70", linewidth = 1) +
  geom_point(aes(y = `Control`), color = "orange", size = 3) +
  geom_point(aes(y = `Water`), color = "darkblue", size = 3) +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  labs(
    x = "Genera",
    y = "Strain-sharing rate",
    title = "Western Kenya",
    subtitle = "Comparison between control and chlorinated water groups"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 12),
    plot.title   = element_text(hjust = 0.5),
    legend.position = "none"
  )


# formatting into a wider format
wit.bet.test.wide.urban <- wit.bet.test.m %>%
  filter(area == "Urban area", tr != "Overall", genus != "Campylobacter_dep") %>%
  select(genus, variable, tr, value) %>%
  pivot_wider(names_from = tr, values_from = value)

# plot
hh.rates.urban.dub <- ggplot(wit.bet.test.wide.urban, aes(x = genus)) +
  geom_segment(aes(y = `Control`, yend = `Water`, xend = genus),
               color = "gray70", linewidth = 1) +
  geom_point(aes(y = `Control`), color = "orange", size = 3) +
  geom_point(aes(y = `Water`), color = "darkblue", size = 3) +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  labs(
    x = "Genera",
    y = "Strain-sharing rate",
    title = "Nairobi",
    subtitle = "Comparison between control and chlorinated water groups"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 12),
    plot.title   = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Save into files
ggsave(file.path(here::here(),"data/figures/hh.rates.tr.dub.pdf"), dpi = 300, scale = 0.3, width = 600, height = 800, units = "mm", plot = ggarrange(hh.rates.rural.dub, hh.rates.urban.dub, nrow = 2))

#### Dumbbell style plots #################### CI version
# Prepare data in wide format
wit.bet.test.wide.rural.CI <- wit.bet.test.CI.m %>%
  filter(area == "Rural area", tr != "Overall", genus != "Campylobacter_dep") %>%
  mutate(genus = droplevels(genus)) %>%
  select(genus, variable, tr, value, lower, upper) %>%
  pivot_wider(names_from = tr, values_from = c(value, lower, upper)) %>%
  mutate(genus_num = as.numeric(genus))

rural_x_nudge <- 0.2

# Create dumbbell plot (Figure 5C)
hh.rates.rural.CI.dub <- ggplot(wit.bet.test.wide.rural.CI) +
  # geom_segment(aes(x = genus_num - rural_x_nudge, xend = genus_num + rural_x_nudge,
  #                   y = value_Control, yend = value_Water),
  #              color = "gray70", linewidth = 1) +
  geom_errorbar(aes(x = genus_num - rural_x_nudge, ymin = lower_Control, ymax = upper_Control),
                color = "orange", width = 0.1) +
  geom_errorbar(aes(x = genus_num + rural_x_nudge, ymin = lower_Water, ymax = upper_Water),
                color = "darkblue", width = 0.1) +
  geom_point(aes(x = genus_num - rural_x_nudge, y = value_Control), color = "orange", size = 3) +
  geom_point(aes(x = genus_num + rural_x_nudge, y = value_Water), color = "darkblue", size = 3) +
  scale_x_continuous(breaks = seq_along(levels(wit.bet.test.wide.rural.CI$genus)),
                      labels = levels(wit.bet.test.wide.rural.CI$genus)) +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  labs(
    x = "Genera",
    y = "Strain-sharing rate",
    title = "Western Kenya",
    subtitle = "Comparison between control and chlorinated water groups"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 12),
    plot.title   = element_text(hjust = 0.5),
    legend.position = "none"
  )


# formatting into a wider format
wit.bet.test.wide.urban.CI <- wit.bet.test.CI.m %>%
  filter(area == "Urban area", tr != "Overall", genus != "Campylobacter_dep") %>%
  mutate(genus = droplevels(genus)) %>%
  select(genus, variable, tr, value, lower, upper) %>%
  pivot_wider(names_from = tr, values_from = c(value, lower, upper)) %>%
  mutate(genus_num = as.numeric(genus))

urban_x_nudge <- 0.2

# plot (Figure 5C)
hh.rates.urban.CI.dub <- ggplot(wit.bet.test.wide.urban.CI) +
  # geom_segment(aes(x = genus_num - urban_x_nudge, xend = genus_num + urban_x_nudge,
  #                   y = value_Control, yend = value_Water),
  #              color = "gray70", linewidth = 1) +
  geom_errorbar(aes(x = genus_num - urban_x_nudge, ymin = lower_Control, ymax = upper_Control),
                color = "orange", width = 0.1) +
  geom_errorbar(aes(x = genus_num + urban_x_nudge, ymin = lower_Water, ymax = upper_Water),
                color = "darkblue", width = 0.1) +
  geom_point(aes(x = genus_num - urban_x_nudge, y = value_Control), color = "orange", size = 3) +
  geom_point(aes(x = genus_num + urban_x_nudge, y = value_Water), color = "darkblue", size = 3) +
  scale_x_continuous(breaks = seq_along(levels(wit.bet.test.wide.urban.CI$genus)),
                      labels = levels(wit.bet.test.wide.urban.CI$genus)) +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  labs(
    x = "Genera",
    y = "Strain-sharing rate",
    title = "Nairobi",
    subtitle = "Comparison between control and chlorinated water groups"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 12),
    plot.title   = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Save into files
ggsave(file.path(here::here(),"data/figures/hh.rates.tr.CI.dub.pdf"), dpi = 300, scale = 0.3, width = 600, height = 800, units = "mm", plot = ggarrange(hh.rates.rural.CI.dub, hh.rates.urban.CI.dub, nrow = 2))
