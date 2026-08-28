# Figure S6
source(file.path(here::here(),"0-config.R"))
#Load source data
wit.bet.test.filter = read.csv(file.path(here::here(), "data/stats/wit.bet.test.filter.csv"))

wit.bet.test.filter.m = wit.bet.test.filter[,c("area", "tr", "within", "between", "genus")] %>% melt()
wit.bet.test.filter.m$tr = factor(wit.bet.test.filter.m$tr, levels = c("Overall", "Water", "Control"))
wit.bet.test.filter.m$area = factor(wit.bet.test.filter.m$area, levels = c("rural","urban"), labels = c("Rural area", "Urban area"))
wit.bet.test.filter.m$variable = factor(wit.bet.test.filter.m$variable, levels = c("within", "between"), labels = c("Within households", "Between households"))
# wit.bet.test.m$genus = factor(wit.bet.test.m$genus, levels = c("All", "comm","non_comm","Bact", "Bifi", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb"), 
#                               labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))
wit.bet.test.filter.m$genus = factor(wit.bet.test.filter.m$genus, levels = c("All", "non_comm","comm", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb", "Bact", "Bifi"), 
                              labels = c("All genera", "Non-commensal", "Commensal", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter", "Bacteroides", "Bifidobacterium"))

wit.bet.test.filter.m[, "split"] = "each"
wit.bet.test.filter.m[wit.bet.test.filter.m$genus %in% c("All genera", "Commensal", "Non-commensal"), "split"] = "combined"

# Urban vs. Rural
# V1
hh.rates.area.filtered = ggplot(wit.bet.test.filter.m %>% filter(tr == "Overall" & genus != "Campylobacter_dep"), aes(x=genus, y=value, fill=area)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  scale_fill_manual(values=c("Urban area"="darkblue", "Rural area"="orange")) +  # Custom colors for Urban and Rural
  labs(x="Genera", y="Strain-sharing rate", fill="", title = "Urban vs. Rural areas\n(chilren & older children)") +
  facet_wrap(~ variable , scales = "free", nrow = 1) +
  #facet_grid(~ variable + split, scales = "free_x", space = "free_x") + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=12),
        plot.title = element_text(hjust=0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)) +
  guides(fill=guide_legend(title.position="top", title.hjust=0.5)) 


# For asterisks
rural.urban.test.filter = read.csv(file.path(here::here(), "data/stats/rural.urban.test.filter.csv")) %>% filter(tr == "Overall")
rural.urban.test.filter$variable = factor(rural.urban.test.filter$HH_type, levels = c("Within", "Between"), labels = c("Within households", "Between households"))
#rural.urban.test$genus = factor(rural.urban.test$genus, levels = c("All", "comm","non_comm","Bact", "Bifi", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb"), 
#                                labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))
rural.urban.test.filter$genus = factor(rural.urban.test.filter$genus, levels = c("All", "non_comm","comm", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb", "Bact", "Bifi"), 
                                labels = c("All genera", "Non-commensal", "Commensal", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter", "Bacteroides", "Bifidobacterium"))

rural.urban.test.filter = rural.urban.test.filter %>% mutate(label = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*", "")))) 

hh.rates.area.filtered = hh.rates.area.filtered +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_text(data = rural.urban.test.filter %>% filter(tr == "Overall" & genus != "Campylobacter_dep"),
            aes(x = genus, y = pmax(rural, urban), label = label),
            vjust = -0.1, col = "red",
            size = 6,
            inherit.aes = FALSE) + 
  facet_wrap(~variable, scales = "free")

# Save into files
ggsave(file.path(here::here(),"data/figures/hh.rates.area.filtered.pdf"), dpi = 300, scale = 0.3, width = 700, height = 400, units = "mm", plot = hh.rates.area.filtered)

############### CI Version ##########################################################################
#Load source data
wit.bet.test.CI.filter = read.csv(file.path(here::here(), "data/stats/wit.bet.test.CI.filter.csv"))

wit.bet.test.CI.filter.m = wit.bet.test.CI.filter %>%
  rename(within_value = within, between_value = between) %>%
  select(area, tr, genus, within_value, within_lower, within_upper, between_value, between_lower, between_upper) %>%
  pivot_longer(
    cols = -c(area, tr, genus),
    names_to = c("variable", ".value"),
    names_pattern = "(within|between)_(value|lower|upper)"
  )

wit.bet.test.CI.filter.m$tr = factor(wit.bet.test.CI.filter.m$tr, levels = c("Overall", "Water", "Control"))
wit.bet.test.CI.filter.m$area = factor(wit.bet.test.CI.filter.m$area, levels = c("rural","urban"), labels = c("Rural area", "Urban area"))
wit.bet.test.CI.filter.m$variable = factor(wit.bet.test.CI.filter.m$variable, levels = c("within", "between"), labels = c("Within households", "Between households"))
wit.bet.test.CI.filter.m$genus = factor(wit.bet.test.CI.filter.m$genus, levels = c("All", "non_comm","comm", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb", "Bact", "Bifi"), 
                                labels = c("All genera", "Non-commensal", "Commensal", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter", "Bacteroides", "Bifidobacterium"))

wit.bet.test.CI.filter.m[, "split"] = "each"
wit.bet.test.CI.filter.m[wit.bet.test.CI.filter.m$genus %in% c("All genera", "Commensal", "Non-commensal"), "split"] = "combined"

# Urban vs. Rural
# V1
hh.rates.area.filtered.CI = ggplot(wit.bet.test.CI.filter.m %>% filter(tr == "Overall" & genus != "Campylobacter_dep"), aes(x=genus, y=value, fill=area)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  scale_fill_manual(values=c("Urban area"="darkblue", "Rural area"="orange")) +  # Custom colors for Urban and Rural
  labs(x="Genera", y="Strain-sharing rate", fill="", title = "Urban vs. Rural areas\n(chilren & older children)") +
  facet_wrap(~ variable , scales = "free", nrow = 1) +
  #facet_grid(~ variable + split, scales = "free_x", space = "free_x") + 
  theme_minimal(base_size = 14) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.7),
    width = 0.2
  )+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=12),
        plot.title = element_text(hjust=0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)) +
  guides(fill=guide_legend(title.position="top", title.hjust=0.5)) 


# For asterisks
rural.urban.test.filter = read.csv(file.path(here::here(), "data/stats/rural.urban.test.filter.csv")) %>% filter(tr == "Overall")
rural.urban.test.filter$variable = factor(rural.urban.test.filter$HH_type, levels = c("Within", "Between"), labels = c("Within households", "Between households"))
#rural.urban.test$genus = factor(rural.urban.test$genus, levels = c("All", "comm","non_comm","Bact", "Bifi", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb"), 
#                                labels = c("All genera", "Commensal", "Non-commensal","Bacteroides", "Bifidobacterium", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter"))
rural.urban.test.filter$genus = factor(rural.urban.test.filter$genus, levels = c("All", "non_comm","comm", "Esch", "Camp","Camp_MAG", "Kleb", "Entc", "Entb", "Bact", "Bifi"), 
                                labels = c("All genera", "Non-commensal", "Commensal", "Escherichia", "Campylobacter_dep", "Campylobacter", "Klebsiella", "Enterococcus", "Enterobacter", "Bacteroides", "Bifidobacterium"))

rural.urban.test.filter = rural.urban.test.filter %>% mutate(label = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*", "")))) 

hh.rates.area.filtered.CI = hh.rates.area.filtered.CI +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_text(data = rural.urban.test.filter %>% filter(tr == "Overall" & genus != "Campylobacter_dep"),
            aes(x = genus, y = pmax(rural, urban), label = label),
            vjust = -0.1, col = "red",
            size = 6,
            inherit.aes = FALSE) + 
  facet_wrap(~variable, scales = "free")

# Save into files
ggsave(file.path(here::here(),"data/figures/hh.rates.area.filtered.CI.pdf"), dpi = 300, scale = 0.3, width = 700, height = 400, units = "mm", plot = hh.rates.area.filtered.CI)
