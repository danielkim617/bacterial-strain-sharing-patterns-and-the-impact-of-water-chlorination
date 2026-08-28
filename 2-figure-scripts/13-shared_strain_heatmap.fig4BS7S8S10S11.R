source(file.path(here::here(),"0-config.R"))
# Load data
prop.strain_sharing = readRDS(file.path(here::here(), "data/prop.strain_sharing.rds"))
staingr.org.ML = read.table(file.path(here::here(), "data/StrainGE_ML/staingr.org.ML.tsv"), sep = "\t", header = T)
name_list1 = read.csv(file.path(here::here(), "data/StrainGE_ML/name_list1.csv"))
# Add strain name
staingr.org.ML.name = merge(staingr.org.ML,name_list1, by = "scaffold" , all.x = T) %>% filter(genus != "Camp") # Because we have Camp_MAG
# For Campylobacter MAG sequences
staingr.org.ML.name$organism[is.na(staingr.org.ML.name$organism)] <- staingr.org.ML.name$scaffold[is.na(staingr.org.ML.name$organism)]

#
urban_cols <- c(
  "strain",
  "urban_sharing_pairs",
  "urban_ref_pairs",
  "urban_all_pairs",
  "urban_sharing_rate",
  "urban_sharing_ref_rate",
  "0 - 19 months / 0 - 19 months",
  "19 - 30 months / 19 - 30 months",
  "> 30 months / > 30 months",
  "Adults > 15 years / Adults > 15 years",
  "0 - 19 months / 19 - 30 months",
  "> 30 months / 0 - 19 months",
  "0 - 19 months / Adults > 15 years",
  "> 30 months / 19 - 30 months",
  "19 - 30 months / Adults > 15 years",
  "> 30 months / Adults > 15 years",
  "0 - 19 months / 0 - 19 months_rate",
  "19 - 30 months / 19 - 30 months_rate",
  "> 30 months / > 30 months_rate",
  "Adults > 15 years / Adults > 15 years_rate",
  "0 - 19 months / 19 - 30 months_rate",
  "> 30 months / 0 - 19 months_rate",
  "0 - 19 months / Adults > 15 years_rate",
  "> 30 months / 19 - 30 months_rate",
  "19 - 30 months / Adults > 15 years_rate",
  "> 30 months / Adults > 15 years_rate",
  "0 - 19 months / 0 - 19 months_all_pairs",
  "19 - 30 months / 19 - 30 months_all_pairs",
  "> 30 months / > 30 months_all_pairs",
  "Adults > 15 years / Adults > 15 years_all_pairs",
  "0 - 19 months / 19 - 30 months_all_pairs",
  "> 30 months / 0 - 19 months_all_pairs",
  "0 - 19 months / Adults > 15 years_all_pairs",
  "> 30 months / 19 - 30 months_all_pairs",
  "19 - 30 months / Adults > 15 years_all_pairs",
  "> 30 months / Adults > 15 years_all_pairs",
  "urban_sharing_pairs_dist_median",
  "urban_all_pairs_dist_median",
  "urban_dist_p"
)

prop.strain_sharing.urban = prop.strain_sharing[, urban_cols] %>% filter(!is.na(urban_sharing_pairs_dist_median))

# prop.strain_sharing.urban = prop.strain_sharing %>% select(c("strain", "urban_sharing_pairs", "0 - 19 months / 0 - 19 months", "19 - 30 months / 19 - 30 months", "> 30 months / > 30 months", 
#                                                              "Adults > 15 years / Adults > 15 years", "0 - 19 months / 19 - 30 months", "> 30 months / 0 - 19 months", "0 - 19 months / Adults > 15 years",
#                                                              "> 30 months / 19 - 30 months", "19 - 30 months / Adults > 15 years", "> 30 months / Adults > 15 years",
#                                                              "0 - 19 months / 0 - 19 months_rate", "19 - 30 months / 19 - 30 months_rate", "> 30 months / > 30 months_rate", 
#                                                              "Adults > 15 years / Adults > 15 years_rate", "0 - 19 months / 19 - 30 months_rate", "> 30 months / 0 - 19 months_rate", "0 - 19 months / Adults > 15 years_rate",
#                                                              "> 30 months / 19 - 30 months_rate", "19 - 30 months / Adults > 15 years_rate", "> 30 months / Adults > 15 years_rate",
#                                                              "urban_all_pairs", "urban_sharing_rate","urban_sharing_pairs_dist_median", "urban_all_pairs_dist_median", "urban_dist_p"))  %>% filter(!is.na(urban_sharing_pairs_dist_median))

prop.strain_sharing.urban.m = prop.strain_sharing.urban %>% select("strain", "urban_sharing_rate","0 - 19 months / 0 - 19 months_rate", "19 - 30 months / 19 - 30 months_rate", "> 30 months / > 30 months_rate", 
                                                                   "Adults > 15 years / Adults > 15 years_rate", "0 - 19 months / 19 - 30 months_rate", "> 30 months / 0 - 19 months_rate", "0 - 19 months / Adults > 15 years_rate",
                                                                   "> 30 months / 19 - 30 months_rate", "19 - 30 months / Adults > 15 years_rate", "> 30 months / Adults > 15 years_rate") %>% melt()


rural_cols <- c(
  "strain",
  "rural_sharing_pairs",
  "rural_ref_pairs",
  "rural_all_pairs",
  "rural_sharing_rate",
  "19 - 30 months / 19 - 30 months_r",
  "> 30 months / > 30 months_r",
  "> 30 months / 19 - 30 months_r",
  "19 - 30 months / 19 - 30 months_r_rate",
  "> 30 months / > 30 months_r_rate",
  "> 30 months / 19 - 30 months_r_rate",
  "19 - 30 months / 19 - 30 months_r_all_pairs",
  "> 30 months / > 30 months_r_all_pairs",
  "> 30 months / 19 - 30 months_r_all_pairs",
  "rural_sharing_pairs_dist_median",
  "rural_all_pairs_dist_median",
  "rural_dist_p"
)

prop.strain_sharing.rural = prop.strain_sharing[, rural_cols] %>% filter(!is.na(rural_sharing_pairs_dist_median))
                                                
# prop.strain_sharing.rural = prop.strain_sharing %>% select(c("strain", "rural_sharing_pairs", "19 - 30 months / 19 - 30 months_r", "> 30 months / > 30 months_r" , "> 30 months / 19 - 30 months_r",
#                                                              "19 - 30 months / 19 - 30 months_r_rate", "> 30 months / > 30 months_r_rate" , "> 30 months / 19 - 30 months_r_rate",
#                                                              "rural_all_pairs", "rural_sharing_rate", "rural_sharing_pairs_dist_median", "rural_all_pairs_dist_median", "rural_dist_p")) %>% filter(!is.na(rural_sharing_pairs_dist_median))


prop.strain_sharing.rural.m = prop.strain_sharing.rural %>% select("strain", "rural_sharing_rate", "19 - 30 months / 19 - 30 months_r_rate", "> 30 months / > 30 months_r_rate" , "> 30 months / 19 - 30 months_r_rate") %>% melt()


# Combined data?
prop.strain_sharing.urban.m = prop.strain_sharing.urban.m %>% mutate(area = "urban")
prop.strain_sharing.rural.m = prop.strain_sharing.rural.m %>% mutate(area = "rural")
# Add genus column
prop.strain_sharing.all.m = data.frame(rbind(prop.strain_sharing.urban.m, prop.strain_sharing.rural.m))
prop.strain_sharing.all.m = merge(prop.strain_sharing.all.m, unique(staingr.org.ML.name[,c("organism", "genus")]), by.x = "strain", by.y = "organism", all.x = T)
prop.strain_sharing.all.m$variable = prop.strain_sharing.all.m$variable %>% gsub("_r_rate", "",.)
prop.strain_sharing.all.m$variable = prop.strain_sharing.all.m$variable %>% gsub("_rate", "",.)
prop.strain_sharing.all.m$variable = prop.strain_sharing.all.m$variable %>% gsub("urban_sharing", "all_sharing",.)
prop.strain_sharing.all.m$variable = prop.strain_sharing.all.m$variable %>% gsub("rural_sharing", "all_sharing",.)


# reordering age groups
prop.strain_sharing.all.m$variable = factor(prop.strain_sharing.all.m$variable, levels = c("all_sharing","0 - 19 months / 0 - 19 months","19 - 30 months / 19 - 30 months", "> 30 months / > 30 months", "Adults > 15 years / Adults > 15 years",
                                                                                           "0 - 19 months / 19 - 30 months", "> 30 months / 0 - 19 months", "0 - 19 months / Adults > 15 years",
                                                                                           "> 30 months / 19 - 30 months", "19 - 30 months / Adults > 15 years", "> 30 months / Adults > 15 years"),
                                            labels = c("All pairs","Younger child - Younger child", "Child - Child", "Older child - Older child", "Mother - Mother",
                                                       "Yonger child - Child", "Younger child - Older child", "Younger child - Mother",
                                                       "Child - Older child", "Child - Mother",
                                                       "Older child - Mother"))

### To extract list of strains to build trees for Bifidobacterium and Enterococcus ###############################################
Bifi_list = prop.strain_sharing.all.m %>% filter(genus == "Bifi") %>% select(strain) %>% unique
Bact_list = prop.strain_sharing.all.m %>% filter(genus == "Bact") %>% select(strain) %>% unique

Esch_list = prop.strain_sharing.all.m %>% filter(genus == "Esch") %>% select(strain) %>% unique
Entc_list = prop.strain_sharing.all.m %>% filter(genus == "Entc") %>% select(strain) %>% unique
Camp_MAG_list = prop.strain_sharing.all.m %>% filter(genus == "Camp_MAG") %>% select(strain) %>% unique
Entb_list = prop.strain_sharing.all.m %>% filter(genus == "Entb") %>% select(strain) %>% unique # Exclude Enterobacter for tree because there is only one strains
Kleb_list = prop.strain_sharing.all.m %>% filter(genus == "Kleb") %>% select(strain) %>% unique # Exclude Klebsiella for tree because there are only two strains

# Save into files
write.table(Bifi_list, file.path(here::here(), "data/Bifi_list.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)
write.table(Bact_list, file.path(here::here(), "data/Bact_list.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)

write.table(Esch_list, file.path(here::here(), "data/Esch_list.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)
write.table(Entc_list, file.path(here::here(), "data/Entc_list.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)
write.table(Camp_MAG_list, file.path(here::here(), "data/Camp_MAG_list.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)
##################################################################################################################################
# Rename
prop.strain_sharing.all.m$strain = gsub("Bifi_", "Bifidobacterium ", prop.strain_sharing.all.m$strain)
prop.strain_sharing.all.m$strain = gsub("Bact_", "Bacteroides ", prop.strain_sharing.all.m$strain)
prop.strain_sharing.all.m$strain = gsub("Esch_", "Escherichia ", prop.strain_sharing.all.m$strain)
prop.strain_sharing.all.m$strain = gsub("Kleb_", "Klebsiella  ", prop.strain_sharing.all.m$strain)
prop.strain_sharing.all.m <- prop.strain_sharing.all.m %>%  mutate(strain = if_else(genus == "Entc",sub("^Ente_", "Enterococcus ", strain), strain))
prop.strain_sharing.all.m$strain = gsub("_", " ", prop.strain_sharing.all.m$strain)
# Change to percent scale
prop.strain_sharing.all.m$value = prop.strain_sharing.all.m$value*100
prop.strain_sharing.all.m <- prop.strain_sharing.all.m %>% mutate(value = na_if(value, 0))
# Plot figures
library(gridExtra)  # or patchwork

# 1) compute the overall range of your 'value' column
val_range_urban <- range(prop.strain_sharing.all.m %>% filter(area == "urban") %>% pull(value), na.rm = TRUE)

# 2) define a common scale
common_fill_urban <- scale_fill_gradient2(
  low      = "white",
  high     = "red",
  limits   = val_range_urban,
  oob      = scales::squish,
  na.value = "grey"
)

# 3) build each plot *with* that same scale
p1 <- ggplot(
  prop.strain_sharing.all.m %>% filter(area == "urban", variable != "All pairs",genus %in% c("Bact", "Bifi")),
  aes(x = variable, y = strain, fill = value)
) +
  geom_tile() +
  common_fill_urban +
  theme_minimal() +
  facet_grid(rows = vars(genus), scales = "free_y", space = "free_y") +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.text.y      = element_text(angle = 0),
    panel.spacing     = unit(1, "lines")
  ) +
  labs(title = "Nairobi", x = NULL, y = NULL, fill = "Strain-sharing pairs (%)")

p2 <- ggplot(
  prop.strain_sharing.all.m %>% filter(area == "urban", variable != "All pairs", !(genus %in% c("Bact", "Bifi"))),
  aes(x = variable, y = strain, fill = value)
) +
  geom_tile() +
  common_fill_urban +
  theme_minimal() +
  facet_grid(rows = vars(genus), scales = "free_y", space = "free_y") +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.text.y      = element_text(angle = 0),
    panel.spacing     = unit(1, "lines")
  ) +
  labs(title = "Nairobi", x = NULL, y = NULL, fill = "Strain-sharing pairs (%)")

# 4) arrange them side by side with a single shared legend
urban.heatmap = grid.arrange(p1, p2, ncol = 2) # Figure S7


# Rural
val_range_rural <- range(prop.strain_sharing.all.m %>% filter(area == "rural") %>% pull(value), na.rm = TRUE)
common_fill_rural <- scale_fill_gradient2(
  low      = "white",
  high     = "red",
  limits   = val_range_rural,
  oob      = scales::squish,
  na.value = "grey"
)

p3 <- ggplot(
  prop.strain_sharing.all.m %>% filter(area == "rural", variable != "All pairs",genus %in% c("Bact", "Bifi")),
  aes(x = variable, y = strain, fill = value)
) +
  geom_tile() +
  common_fill_rural +
  theme_minimal() +
  facet_grid(rows = vars(genus), scales = "free_y", space = "free_y") +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.text.y      = element_text(angle = 0),
    panel.spacing     = unit(1, "lines")
  ) +
  labs(title = "Western Kenya", x = NULL, y = NULL, fill = "Strain-sharing pairs (%)")

p4 <- ggplot(
  prop.strain_sharing.all.m %>% filter(area == "rural", variable != "All pairs",!(genus %in% c("Bact", "Bifi"))),
  aes(x = variable, y = strain, fill = value)
) +
  geom_tile() +
  common_fill_rural +
  theme_minimal() +
  facet_grid(rows = vars(genus), scales = "free_y", space = "free_y") +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.text.y      = element_text(angle = 0),
    panel.spacing     = unit(1, "lines")
  ) +
  labs(title = "Western Kenya", x = NULL, y = NULL, fill = "Strain-sharing pairs (%)")

rural.heatmap = grid.arrange(p3, p4, ncol = 2)

# Plot together
grid.arrange(rural.heatmap, urban.heatmap, ncol = 2, widths = c(1, 1.25))

# Save into a file
ggsave(file.path(here::here(),"data/figures/strain_level_sharing.fig.pdf"), dpi = 300, scale = 0.3, width = 2500, height = 1000, units = "mm", plot = grid.arrange(rural.heatmap, urban.heatmap, ncol = 2, widths = c(1, 1.25)))

# Subset of the figure of Bifidobacterium and Enterococus
bifi_entc_urban = ggplot(prop.strain_sharing.all.m %>% filter(area == "urban",variable != "All pairs",genus %in% c("Bifi", "Entc")),aes(x = variable, y = strain, fill = value)) +
  geom_tile() +
  common_fill_urban +
  theme_minimal() +
  facet_grid(rows = vars(genus), scales = "free_y", space = "free_y") +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.text.y      = element_text(angle = 0),
    panel.spacing     = unit(1, "lines")
  ) +
  labs(title = "Nairobi", x = NULL, y = NULL, fill = "Strain-sharing pairs (%)")

bifi_entc_rural = ggplot(prop.strain_sharing.all.m %>% filter(area == "rural", variable != "All pairs",genus %in% c("Bifi", "Entc")), aes(x = variable, y = strain, fill = value)) +
  geom_tile() +
  common_fill_urban + # to have color range consistent with urban
  theme_minimal() +
  facet_grid(rows = vars(genus), scales = "free_y", space = "free_y") +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.text.y      = element_text(angle = 0),
    panel.spacing     = unit(1, "lines")
  ) +
  labs(title = "Western Kenya", x = NULL, y = NULL, fill = "Strain-sharing pairs (%)")

grid.arrange(bifi_entc_rural, bifi_entc_urban, ncol = 2, widths = c(1,1.3))
# Save into a file
ggsave(file.path(here::here(),"data/figures/Bifi_Entc_strain_level_sharing.fig.pdf"), dpi = 300, scale = 0.3, width = 1300, height = 700, units = "mm", plot = grid.arrange(bifi_entc_rural, bifi_entc_urban, ncol = 2, widths = c(1,1.3)))

### Plot heatmaps with trees
library(ggtree)

Bifi_urban = prop.strain_sharing.all.m %>% filter(area == "urban", variable != "All pairs",genus == c("Bifi"))
Bifi_rural = prop.strain_sharing.all.m %>% filter(area == "rural", variable != "All pairs",genus == c("Bifi"))

Bact_urban = prop.strain_sharing.all.m %>% filter(area == "urban", variable != "All pairs",genus == c("Bact"))
Bact_rural = prop.strain_sharing.all.m %>% filter(area == "rural", variable != "All pairs",genus == c("Bact"))

Entc_urban = prop.strain_sharing.all.m %>% filter(area == "urban", variable != "All pairs",genus == c("Entc"))
Entc_rural = prop.strain_sharing.all.m %>% filter(area == "rural", variable != "All pairs",genus == c("Entc"))

Esch_urban = prop.strain_sharing.all.m %>% filter(area == "urban", variable != "All pairs",genus == c("Esch"))
Esch_rural = prop.strain_sharing.all.m %>% filter(area == "rural", variable != "All pairs",genus == c("Esch"))

Camp_MAG_urban = prop.strain_sharing.all.m %>% filter(area == "urban", variable != "All pairs",genus == c("Camp_MAG"))
Camp_MAG_rural = prop.strain_sharing.all.m %>% filter(area == "rural", variable != "All pairs",genus == c("Camp_MAG"))
# into matrix
# Bifi
# into_matrix <- function(x){
#   x %>%
#     select(strain, variable, value) %>%
#     pivot_wider(names_from = variable, values_from = value, values_fill = 0) %>%
#     tibble::column_to_rownames("strain") %>%
#     as.matrix()
# }
into_matrix <- function(x){
  mat <- x %>%
    select(strain, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value, values_fill = 0) %>%
    tibble::column_to_rownames("strain") %>%
    as.matrix()

  mat[, levels(droplevels(x$variable)), drop = FALSE]
}

#
Bifi_urban.mat = into_matrix(Bifi_urban)
Bifi_rural.mat = into_matrix(Bifi_rural)

Bact_urban.mat = into_matrix(Bact_urban)
Bact_rural.mat = into_matrix(Bact_rural)

Entc_urban.mat = into_matrix(Entc_urban)
Entc_rural.mat = into_matrix(Entc_rural)

Esch_urban.mat = into_matrix(Esch_urban)
Esch_rural.mat = into_matrix(Esch_rural)

Camp_MAG_urban.mat = into_matrix(Camp_MAG_urban)
Camp_MAG_rural.mat = into_matrix(Camp_MAG_rural)

# ---- read a tree whose tip labels are the rownames(mat)
Bifi_tree <- ape::read.tree(file.path(here::here(),"data/tree/bifidobacterium.tre"))
Bact_tree <- ape::read.tree(file.path(here::here(),"data/tree/bacteroides.tre"))
Entc_tree <- ape::read.tree(file.path(here::here(),"data/tree/enterococcus.tre"))
Esch_tree <- ape::read.tree(file.path(here::here(),"data/tree/ecoli_new.tre"))
Camp_MAG_tree <- ape::read.tree(file.path(here::here(),"data/tree/Camp_MAG.tre"))
#
Bifi_tree$tip.label = gsub("Bifi_", "Bifidobacterium ", Bifi_tree$tip.label)
Bifi_tree$tip.label = gsub("_", " ", Bifi_tree$tip.label)

Bact_tree$tip.label = gsub("Bact_", "Bacteroides ", Bact_tree$tip.label)
Bact_tree$tip.label = gsub("_", " ", Bact_tree$tip.label)

Entc_tree$tip.label = gsub("Ente_", "Enterococcus ", Entc_tree$tip.label)
Entc_tree$tip.label = gsub("_", " ", Entc_tree$tip.label)

Esch_tree$tip.label = gsub("Esch_", "Escherichia ", Esch_tree$tip.label)
Esch_tree$tip.label = gsub("_", " ", Esch_tree$tip.label)

Camp_MAG_tree$tip.label = gsub("_", " ", Camp_MAG_tree$tip.label)

# Make plots
# Bifi
bifi_tree_p <- ggtree(Bifi_tree) + 
  geom_tiplab(size=3, align=TRUE, linesize=.5) +
  geom_treescale(x=0.5, y=25, offset=0.5)

bifi_tree_p1 <- gheatmap(
  bifi_tree_p, Bifi_rural.mat,
  offset = 0.5,
  width = 0.3,
  colnames = TRUE,
  colnames_angle = 45,
  hjust = 1,
  colnames_offset_y = 0.5,
  font.size = 4,        # smaller column label size
) + common_fill_urban + labs(fill = "Strain-sharing pairs (%)")

bifi_tree_p2 <- gheatmap(
  bifi_tree_p1, Bifi_urban.mat,
  offset = 1.0,          # first offset + first width + gap
  width = 1,
  colnames = TRUE,
  colnames_angle = 45,
  hjust = 1,
  colnames_offset_y = 0.5,
  font.size = 4,
) + common_fill_urban + labs(fill = "Strain-sharing pairs (%)")

# Entc
entc_tree_p <- ggtree(Entc_tree) + 
  geom_tiplab(size=3, align=TRUE, linesize=.5) +
  geom_treescale(x=0.5, y=12, offset=0.5) 

entc_tree_p1 <- gheatmap(
  entc_tree_p, Entc_rural.mat,
  offset = 0.5,
  width = 0.3,
  colnames = TRUE,
  colnames_angle = 45,
  hjust = 1,
  colnames_offset_y = 0.5,
  font.size = 4,        # smaller column label size
) + common_fill_urban + labs(fill = "Strain-sharing pairs (%)")

entc_tree_p2 <- gheatmap(
  entc_tree_p1, Entc_urban.mat,
  offset = 1.0 + 0.1,          # first offset + first width + gap
  width = 1,
  colnames = TRUE,
  colnames_angle = 45,
  hjust = 1,
  colnames_offset_y = 0.5,
  font.size = 4,
) + common_fill_urban + labs(fill = "Strain-sharing pairs (%)")

library(patchwork)
p_combined <- (bifi_tree_p2 / entc_tree_p2) +
  plot_layout(heights = c(2, 1), guides = "collect") &
  theme(legend.position = "right") &
  coord_cartesian(clip = "off")   # avoid clipping long labels

# Save into a file Fig 4B
ggsave(file.path(here::here(),"data/figures/Bifi_Entc_strain_level_sharing_w_tree.fig.pdf"), dpi = 300, scale = 0.3, width = 800, height = 600, units = "mm", plot = p_combined)

## additional targets with trees
tree_heatmap_fig <- function(tree_file, matrix_file){
  # 1) compute the overall range of your 'value' column
  val_range <- range(prop.strain_sharing.all.m %>% filter(area == "urban") %>% pull(value), na.rm = TRUE)
  
  # 2) define a common scale
  common_fill <- scale_fill_gradient2(
    low      = "white",
    high     = "red",
    limits   = val_range,
    oob      = scales::squish,
    na.value = "grey"
  )
  
  tree_p <- ggtree(tree_file) + 
    geom_tiplab(size=3, align=TRUE, linesize=.5) +
    geom_treescale(x=0.5, y=25, offset=0.5) +
    xlim_tree(15) 
  
  tree_p1 <- gheatmap(
    tree_p, matrix_file,
    offset = 2.5,
    width = 0.3,
    colnames = TRUE,
    colnames_angle = 45,
    hjust = 1,
    colnames_offset_y = 0.5,
    font.size = 4,        # smaller column label size
  ) + common_fill + labs(fill = "Strain-sharing pairs (%)")
  
  return(tree_p1)
}

# Western Kenya
Bifi_rural_heatmap = tree_heatmap_fig(tree_file = keep.tip(Bifi_tree, row.names(Bifi_rural.mat)), matrix_file = Bifi_rural.mat)
Bact_rural_heatmap = tree_heatmap_fig(tree_file = keep.tip(Bact_tree, row.names(Bact_rural.mat)), matrix_file = Bact_rural.mat)
Esch_rural_heatmap = tree_heatmap_fig(tree_file = keep.tip(Esch_tree, row.names(Esch_rural.mat)), matrix_file = Esch_rural.mat)
Camp_MAG_rural_heatmap = tree_heatmap_fig(tree_file = keep.tip(Camp_MAG_tree, row.names(Camp_MAG_rural.mat)), matrix_file = Camp_MAG_rural.mat)
# Nairobi
Bifi_urban_heatmap = tree_heatmap_fig(tree_file = keep.tip(Bifi_tree, row.names(Bifi_urban.mat)), matrix_file = Bifi_urban.mat)
Bact_urban_heatmap = tree_heatmap_fig(tree_file = keep.tip(Bact_tree, row.names(Bact_urban.mat)), matrix_file = Bact_urban.mat)
Esch_urban_heatmap = tree_heatmap_fig(tree_file = keep.tip(Esch_tree, row.names(Esch_urban.mat)), matrix_file = Esch_urban.mat)
Camp_MAG_urban_heatmap = tree_heatmap_fig(tree_file = keep.tip(Camp_MAG_tree, row.names(Camp_MAG_urban.mat)), matrix_file = Camp_MAG_urban.mat)

# Hmm... as plots in the main text having studies side by side by target
tree_heatmap_study_fig <- function(tree_file, matrix_file1, matrix_file2){
  # 1) compute the overall range of your 'value' column
  val_range <- range(prop.strain_sharing.all.m %>% pull(value), na.rm = TRUE)
  
  # 2) define a common scale
  common_fill <- scale_fill_gradient2(
    low      = "white",
    high     = "red",
    limits   = val_range,
    oob      = scales::squish,
    na.value = "grey"
  )
  tree_p <- ggtree(tree_file) + 
    geom_tiplab(size=3, align=TRUE, linesize=.5) +
    geom_treescale(x=0.5, y=12, offset=0.5)
  
  tree_p1 <- gheatmap(
    tree_p, matrix_file1,
    offset = 1,
    width = 0.5,
    colnames = TRUE,
    colnames_angle = 45,
    hjust = 1,
    colnames_offset_y = 0.5,
    font.size = 4,        # smaller column label size
  ) + common_fill + labs(fill = "Strain-sharing pairs (%)")
  
  tree_p2 <- gheatmap(
    tree_p1, matrix_file2,
    offset = 1.0 + 0.8,          # first offset + first width + gap
    width = 0.25,
    colnames = TRUE,
    colnames_angle = 45,
    hjust = 1,
    colnames_offset_y = 0.5,
    font.size = 4,
  ) + common_fill + labs(fill = "Strain-sharing pairs (%)")
  
  return(tree_p2)
}

# in addition to Bifidobacterium and Enterococcus
Bact_heatmap = tree_heatmap_study_fig(tree_file = Bact_tree, matrix_file1 = Bact_urban.mat, matrix_file2 = Bact_rural.mat)
Esch_urban_heatmap = tree_heatmap_study_fig(tree_file = keep.tip(Esch_tree, c(row.names(Esch_urban.mat), row.names(Esch_rural.mat)) %>% unique), matrix_file1 = Esch_urban.mat, matrix_file2 = Esch_rural.mat)
Camp_MAG_urban_heatmap = tree_heatmap_study_fig(tree_file = Camp_MAG_tree, matrix_file1 = Camp_MAG_urban.mat, matrix_file2 = Camp_MAG_rural.mat)

################################################################################################################################################
# Heatmaps for within household sharing
urban.all.pairs = readRDS(file.path(here::here(), "data/urban.all.pairs.rds"))
rural.all.pairs = readRDS(file.path(here::here(), "data/rural.all.pairs.rds"))
R21_sample_list.rural = readRDS(file.path(here::here(), "data/R21_sample_list.rural.rds"))
R21_sample_list.urban = readRDS(file.path(here::here(), "data/R21_sample_list.urban.rds"))

# Add Age_group information
urban.all.pairs[,"sample1_Age_group"] = R21_sample_list.urban[match(urban.all.pairs$sample1, R21_sample_list.urban$sample_id), c("Age_group")] 
urban.all.pairs[,"sample2_Age_group"] = R21_sample_list.urban[match(urban.all.pairs$sample2, R21_sample_list.urban$sample_id), "Age_group"]
rural.all.pairs[,"sample1_Age_group"] = R21_sample_list.rural[match(rural.all.pairs$sample1, R21_sample_list.rural$sample_id), c("Age_group")] 
rural.all.pairs[,"sample2_Age_group"] = R21_sample_list.rural[match(rural.all.pairs$sample2, R21_sample_list.rural$sample_id), "Age_group"]
# Add sorted Age group pair information
urban.all.pairs <- urban.all.pairs %>%
  rowwise() %>%
  mutate(Age_group_pair = paste(sort(c(as.character(sample1_Age_group), as.character(sample2_Age_group))), collapse = " / ")) %>%
  ungroup()
rural.all.pairs <- rural.all.pairs %>%
  rowwise() %>%
  mutate(Age_group_pair = paste(sort(c(as.character(sample1_Age_group), as.character(sample2_Age_group))), collapse = " / ")) %>%
  ungroup()

urban.all.pairs %>% filter(HH_type == "Within") %>% select(Age_group_pair) %>% unique
urban.all.pairs %>% filter(HH_type == "Within") %>% select(Age_group_pair) %>% table

rural.all.pairs %>% filter(HH_type == "Within") %>% select(Age_group_pair) %>% unique
rural.all.pairs %>% filter(HH_type == "Within") %>% select(Age_group_pair) %>% table
# read data
prop.strain_sharing_within = readRDS(file.path(here::here(), "data/prop.strain_sharing_within.rds"))

prop.strain_sharing_within.urban.m = prop.strain_sharing_within %>% filter(urban_sharing_pairs > 0) %>% select("strain", urban.all.pairs %>% filter(HH_type == "Within") %>% pull(Age_group_pair) %>% unique %>% paste0(., "_rate")) %>% melt()
                                                                   #  "0 - 19 months / 0 - 19 months_rate", "19 - 30 months / 19 - 30 months_rate", "> 30 months / > 30 months_rate", 
                                                                   # "Adults > 15 years / Adults > 15 years_rate", "0 - 19 months / 19 - 30 months_rate", "> 30 months / 0 - 19 months_rate", "0 - 19 months / Adults > 15 years_rate",
                                                                   # "> 30 months / 19 - 30 months_rate", "19 - 30 months / Adults > 15 years_rate", "> 30 months / Adults > 15 years_rate") %>% melt()

prop.strain_sharing_within.rural.m = prop.strain_sharing_within %>% filter(rural_sharing_pairs > 0) %>% select("strain", rural.all.pairs %>% filter(HH_type == "Within") %>% pull(Age_group_pair) %>% unique %>% paste0(., "_r_rate")) %>% melt()
                                                                    #"19 - 30 months / 19 - 30 months_r_rate", "> 30 months / > 30 months_r_rate" , "> 30 months / 19 - 30 months_r_rate") %>% melt()

# Combined data?
prop.strain_sharing_within.urban.m = prop.strain_sharing_within.urban.m %>% mutate(area = "urban")
prop.strain_sharing_within.rural.m = prop.strain_sharing_within.rural.m %>% mutate(area = "rural")
# Add genus column
prop.strain_sharing.all.within.m = data.frame(rbind(prop.strain_sharing_within.urban.m, prop.strain_sharing_within.rural.m))
prop.strain_sharing.all.within.m = merge(prop.strain_sharing.all.within.m, unique(staingr.org.ML.name[,c("organism", "genus")]), by.x = "strain", by.y = "organism", all.x = T)
prop.strain_sharing.all.within.m$variable = prop.strain_sharing.all.within.m$variable %>% gsub("_r_rate", "",.)
prop.strain_sharing.all.within.m$variable = prop.strain_sharing.all.within.m$variable %>% gsub("_rate", "",.)

# reordering age groups
prop.strain_sharing.all.within.m$variable = factor(prop.strain_sharing.all.within.m$variable, levels = c("0 - 19 months / 0 - 19 months","19 - 30 months / 19 - 30 months", "> 30 months / > 30 months", "Adults > 15 years / Adults > 15 years",
                                                                                           "0 - 19 months / 19 - 30 months", "> 30 months / 0 - 19 months", "0 - 19 months / Adults > 15 years",
                                                                                           "> 30 months / 19 - 30 months", "19 - 30 months / Adults > 15 years", "> 30 months / Adults > 15 years"),
                                            labels = c("Younger child - Younger child", "Child - Child", "Older child - Older child", "Mother - Mother",
                                                       "Yonger child - Child", "Younger child - Older child", "Younger child - Mother",
                                                       "Child - Older child", "Child - Mother",
                                                       "Older child - Mother"))
# Rename
prop.strain_sharing.all.within.m$strain = gsub("Bifi_", "Bifidobacterium ", prop.strain_sharing.all.within.m$strain)
prop.strain_sharing.all.within.m$strain = gsub("Bact_", "Bacteroides ", prop.strain_sharing.all.within.m$strain)
prop.strain_sharing.all.within.m$strain = gsub("Esch_", "Escherichia  ", prop.strain_sharing.all.within.m$strain)
prop.strain_sharing.all.within.m$strain = gsub("Kleb_", "Klebsiella  ", prop.strain_sharing.all.within.m$strain)
prop.strain_sharing.all.within.m <- prop.strain_sharing.all.within.m %>%  mutate(strain = if_else(genus == "Entc",sub("^Ente_", "Enterococcus ", strain), strain))
prop.strain_sharing.all.within.m$strain = gsub("_", " ", prop.strain_sharing.all.within.m$strain)

# Change to percent scale
prop.strain_sharing.all.within.m$value = prop.strain_sharing.all.within.m$value*100
prop.strain_sharing.all.within.m <- prop.strain_sharing.all.within.m %>% mutate(value = na_if(value, 0))

# Age_group_pair
# > 30 months / > 30 months        > 30 months / 0 - 19 months       > 30 months / 19 - 30 months    > 30 months / Adults > 15 years 
# 24                                 27                                 15                                 88 
# 0 - 19 months / 0 - 19 months     0 - 19 months / 19 - 30 months  0 - 19 months / Adults > 15 years 19 - 30 months / Adults > 15 years 
# 1                                  3                                 42                                 25 

# rural.all.pairs %>% filter(HH_type == "Within") %>% select(Age_group_pair) %>% table
# Age_group_pair
# > 30 months / 19 - 30 months 19 - 30 months / 19 - 30 months 
# 116                               1 

# Plot figures
# 1) compute the overall range of your 'value' column
val_range_urban <- range(prop.strain_sharing.all.within.m %>% filter(area == "urban", !(variable %in% c("Younger child - Younger child", "Yonger child - Child"))) %>% pull(value), na.rm = TRUE)

# 2) define a common scale
common_fill_urban <- scale_fill_gradient2(
  low      = "white",
  high     = "red",
  limits   = val_range_urban,
  oob      = scales::squish,
  na.value = "grey"
)

# 3) build each plot *with* that same scale
p1.w <- ggplot(
  prop.strain_sharing.all.within.m %>% filter(area == "urban", !(variable %in% c("Younger child - Younger child", "Yonger child - Child")),genus %in% c("Bact", "Bifi")),
  aes(x = variable, y = strain, fill = value)
) +
  geom_tile() +
  common_fill_urban +
  theme_minimal() +
  facet_grid(rows = vars(genus), scales = "free_y", space = "free_y") +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.text.y      = element_text(angle = 0),
    panel.spacing     = unit(1, "lines")
  ) +
  labs(title = "Nairobi", x = NULL, y = NULL, fill = "Strain-sharing pairs (%)")

p2.w <- ggplot(
  prop.strain_sharing.all.within.m %>% filter(area == "urban", !(variable %in% c("Younger child - Younger child", "Yonger child - Child")),!(genus %in% c("Bact", "Bifi"))),
  aes(x = variable, y = strain, fill = value)
) +
  geom_tile() +
  common_fill_urban +
  theme_minimal() +
  facet_grid(rows = vars(genus), scales = "free_y", space = "free_y") +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.text.y      = element_text(angle = 0),
    panel.spacing     = unit(1, "lines")
  ) +
  labs(title = "Nairobi", x = NULL, y = NULL, fill = "Strain-sharing pairs (%)")

# 4) arrange them side by side with a single shared legend
urban.w.heatmap = grid.arrange(p1.w, p2.w, ncol = 2,  widths = c(0.85, 1))

# Rural
val_range_rural <- range(prop.strain_sharing.all.within.m %>% filter(area == "rural", variable != "Child - Child") %>% pull(value), na.rm = TRUE)
common_fill_rural <- scale_fill_gradient2(
  low      = "white",
  high     = "red",
  limits   = val_range_rural,
  oob      = scales::squish,
  na.value = "grey"
)

p3.w <- ggplot(
  prop.strain_sharing.all.within.m %>% filter(area == "rural", variable != "Child - Child", genus %in% c("Bact", "Bifi")),
  aes(x = variable, y = strain, fill = value)
) +
  geom_tile() +
  common_fill_rural +
  theme_minimal() +
  facet_grid(rows = vars(genus), scales = "free_y", space = "free_y") +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.text.y      = element_text(angle = 0),
    panel.spacing     = unit(1, "lines")
  ) +
  labs(title = "Western Kenya", x = NULL, y = NULL, fill = "Strain-sharing pairs (%)")

p4.w <- ggplot(
  prop.strain_sharing.all.within.m %>% filter(area == "rural",variable != "Child - Child",!(genus %in% c("Bact", "Bifi"))),
  aes(x = variable, y = strain, fill = value)
) +
  geom_tile() +
  common_fill_rural +
  theme_minimal() +
  facet_grid(rows = vars(genus), scales = "free_y", space = "free_y") +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.text.y      = element_text(angle = 0),
    panel.spacing     = unit(1, "lines")
  ) +
  labs(title = "Western Kenya", x = NULL, y = NULL, fill = "Strain-sharing pairs (%)")

rural.w.heatmap = grid.arrange(p3.w, p4.w, ncol = 2, widths = c(1, 1))

# Plot together
grid.arrange(rural.w.heatmap, urban.w.heatmap, ncol = 2, widths = c(1, 1.5))

# Save into a file
ggsave(file.path(here::here(),"data/figures/strain_level_sharing_within.fig.pdf"), dpi = 300, scale = 0.3, width = 2500, height = 1000, units = "mm", plot = grid.arrange(rural.w.heatmap, urban.w.heatmap, ncol = 2, widths = c(1, 1.5)))
