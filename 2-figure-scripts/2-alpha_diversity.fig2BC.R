source(file.path(here::here(),"0-config.R"))
library(rstatix)
#Load R object
Kraken.phy.genus.alpha = readRDS(file.path(here::here(), "data/Kraken.phy.genus.alpha.rds"))
Kraken.phy.genus.alpha = Kraken.phy.genus.alpha %>% filter(sample_id != "stools_12093101_SM-NA1P5")

#Plot for species richness################################################################################
species_rich = ggplot(Kraken.phy.genus.alpha, aes(x=Age/30, y=Observed, fill=Age/30)) + 
  geom_point(shape=21) +
  scale_fill_viridis_c(option = "D", direction = -1) +
  facet_grid(.~area) +
  labs(x="Age (Months)", y="Observed genera", title="Comparison of species richness", fill="Age (Months)") +
  theme_bw() +
  theme(text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=16)) +
  geom_vline(xintercept = 19, linetype="dashed", color="blue", size = 0.5) +
  geom_vline(xintercept = 30, linetype="dashed", color="blue", size = 0.5)
###################################################################################################
#Plot for Shannon ################################################################################
shannon = ggplot(Kraken.phy.genus.alpha, aes(x=Age/30, y=Shannon, fill=Age/30)) + 
  geom_point(shape=21) +
  scale_fill_viridis_c(option = "D", direction = -1) +
  facet_grid(.~area) +
  labs(x="Age (Months)", y="Shannon diversity", title="Comparison of Shannon diversity", fill="Age (Months)") +
  theme_bw() +
  theme(text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=16)) +
  geom_vline(xintercept = 19, linetype="dashed", color="blue", size = 0.5) +
  geom_vline(xintercept = 30, linetype="dashed", color="blue", size = 0.5)
###################################################################################################
ggsave(file.path(here::here(),"data/figures/alpha_diverity_age.pdf"), dpi = 300, scale = 0.3, width = 800, height = 550, units = "mm", plot = ggarrange(species_rich, shannon, ncol = 1, align = c("hv")))

###################################################################################################
# Grouping by ages (control vs. water)
# load stat data
shannon.test.area.age.tr = read.csv(file.path(here::here(),"data/stats/shannon.test.area.age.tr.csv"))

stat.test <- Kraken.phy.genus.alpha[!is.na(Kraken.phy.genus.alpha$Age_group),] %>%
  group_by(area, Age_group) %>%
  wilcox_test(Shannon ~ tr) %>% add_xy_position(x = "Age_group", dodge = 0.8, scales = "fixed") %>%
  group_by(area) %>%
  adjust_pvalue(method = "BH") %>% mutate(p.adj = round(p.adj, digits = 3))

stat.test$x = c(1,2,1,2,3,4)
stat.test$xmin = stat.test$x - 0.2
stat.test$xmax = stat.test$x + 0.2

alpha.treatment = ggplot(Kraken.phy.genus.alpha[!is.na(Kraken.phy.genus.alpha$Age_group),], aes(x=Age_group, y=Shannon, col=tr)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.5) +
  scale_fill_manual(values =  c("#56B4E9", "#009E73")) +
  facet_grid(.~ area, scales = "free", space = "free") +
  labs(x="Age group", y="Shannon diversity", title="Comparison of Shannon diversity", fill="Treatment") +
  theme_bw() +
  #stat_pvalue_manual(shannon.test.area.age.tr1, x = "Age_group", group1 = "group1", group2 = "group2") + 
  stat_pvalue_manual(stat.test) + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.01, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10), text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=16))

ggsave(file.path(here::here(),"data/figures/alpha.treatment.pdf"), dpi = 300, scale = 0.3, width = 650, height = 500, units = "mm", plot = alpha.treatment)
###################################################################################################
# age comparison
compare.shannon.age = read.csv(file.path(here::here(),"data/stats/compare.shannon.age.csv"))

stat.test2 <- Kraken.phy.genus.alpha[!is.na(Kraken.phy.genus.alpha$Age_group),] %>%
  group_by(area) %>%
  wilcox_test(Shannon ~ Age_group) %>% add_xy_position(x = "Age_group", dodge = 0.8, scales = "fixed") %>%
  group_by(area) %>%
  adjust_pvalue(method = "BH") 
# adjust position
stat.test2[1, "xmin"] = 1
stat.test2[1, "xmax"] = 2

# Figure 2B
alpha.age = ggplot(Kraken.phy.genus.alpha[!is.na(Kraken.phy.genus.alpha$Age_group),], aes(x=Age_group, y=Shannon, col=Age_group)) +
  geom_boxplot() +
  #geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.5) +
  geom_jitter(alpha = 0.5) + 
  scale_fill_manual(values =  c("#56B4E9", "#009E73")) +
  facet_grid(.~ area, scales = "free", space = "free") +
  labs(x="Age group", y="Shannon diversity", title="Comparison of Shannon diversity", fill="Treatment") +
  theme_bw() +
  stat_pvalue_manual(stat.test2, hide.ns = TRUE,  tip.length = 0.01, label = "p.adj", step.increase = 0.02) + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.01, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10), text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=16))

ggsave(file.path(here::here(),"data/figures/alpha.age.pdf"), dpi = 300, scale = 0.3, width = 600, height = 500, units = "mm", plot = alpha.age)

# Grouping by ages (urban vs. rural)
#laod stat data
stat.test3 <- Kraken.phy.genus.alpha[!is.na(Kraken.phy.genus.alpha$Age_group) & Kraken.phy.genus.alpha$Age_group %in% c("19 - 30 months", "> 30 months"),] %>%
  group_by(Age_group) %>%
  wilcox_test(Shannon ~ area) %>% add_xy_position(x = "Age_group", dodge = 0.8, scales = "fixed") %>%
  adjust_pvalue(method = "BH") 

stat.test3$x = c(1,2)
stat.test3$xmin = stat.test3$x - 0.2
stat.test3$xmax = stat.test3$x + 0.2

# Figure 2C
alpha.area = ggplot(Kraken.phy.genus.alpha[!is.na(Kraken.phy.genus.alpha$Age_group) & Kraken.phy.genus.alpha$Age_group %in% c("19 - 30 months", "> 30 months"),], aes(x=Age_group, y=Shannon, col=area)) +
  geom_boxplot(width = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), alpha = 0.5) +
  scale_fill_manual(values =  c("#56B4E9", "#009E73")) +
  labs(x="Age group", y="Shannon diversity", title="Comparison of Shannon diversity", fill="Area") +
  stat_pvalue_manual(stat.test3, hide.ns = TRUE,  tip.length = 0.01, label = "p.adj", step.increase = 0.02) + 
  theme_bw() +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.01, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10), text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=16))

ggsave(file.path(here::here(),"data/figures/alpha.area.pdf"), dpi = 300, scale = 0.3, width = 450, height = 450, units = "mm", plot = alpha.area)
