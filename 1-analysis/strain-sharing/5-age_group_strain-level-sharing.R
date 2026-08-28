# This scripts is to calculate the proportion of strain-sharing pairs among all sample pairs for each strain identified from StrainGST. The output will be used for the heatmap figure.
source(file.path(here::here(),"0-config.R"))
# Load data
R21_sample_list = readRDS(file.path(here::here(), "data/Kraken.phy.genus.alpha.rds"))
staingr.org.ML = read.table(file.path(here::here(), "data/StrainGE_ML/staingr.org.ML.tsv"), sep = "\t", header = T)
name_list = process_nuccore_file(file.path(here::here(),"data/nuccore_result.txt"))
all = readRDS(file.path(here::here(), "data/straingst_comb/all_strainGST.tsv"))
all.ext = merge(all, R21_sample_list, by.x = "V1", by.y = "sample_id", all.x = T)

# For editing names to match one from StrainGE
# write.csv(name_list, file.path(here::here(), "data/StrainGE_ML/name_list.csv"), row.names = F)
name_list1 = read.csv(file.path(here::here(), "data/StrainGE_ML/name_list1.csv"))
# Add strain name
staingr.org.ML.name = merge(staingr.org.ML,name_list1, by = "scaffold" , all.x = T) %>% filter(genus != "Camp") # Because we have Camp_MAG
# For Campylobacter MAG sequences
staingr.org.ML.name$organism[is.na(staingr.org.ML.name$organism)] <- staingr.org.ML.name$scaffold[is.na(staingr.org.ML.name$organism)]

################################################################################################################
##. All sample pair information
R21_sample_list.rural = readRDS(file.path(here::here(), "data/R21_sample_list.rural.rds"))
R21_sample_list.urban = readRDS(file.path(here::here(), "data/R21_sample_list.urban.rds"))
urban.all.pairs = readRDS(file.path(here::here(), "data/urban.all.pairs.rds"))
rural.all.pairs = readRDS(file.path(here::here(), "data/rural.all.pairs.rds"))
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

##### strain-sharing case / number of all pairs with the same StrainGE reference
# ESO1212705120S1K.MAG.40.renamed --> ESO1212705120S1K_k141_16
# EST1361007101S1K.MAG.8.renamed --> EST1361007101S1K_k141_55499
all.ext.rename = all.ext %>%
  mutate(
    V3 = gsub(".fa.gz", "", V3),  # Your original line
    V3 = case_when(
      V3 == "ESO1212705120S1K.MAG.40.renamed" ~ "ESO1212705120S1K_k141_16",
      V3 == "EST1361007101S1K.MAG.8.renamed"  ~ "EST1361007101S1K_k141_55499",
      TRUE ~ V3  # Keep original value if no match
    )
  )
# Add age group pair category data
staingr.org.ML.name <- staingr.org.ML.name %>%
  rowwise() %>%
  mutate(Age_group_pair = paste(sort(c(as.character(sample1_Age_group), as.character(sample2_Age_group))), collapse = " / ")) %>%
  ungroup()

prop.strain_sharing = data.frame(matrix(ncol = 55))
colnames(prop.strain_sharing) <- c(
  # 1) Strain identifier
  "strain",
  
  # 2) Urban summary counts & rates
  "urban_sharing_pairs",
  "urban_ref_pairs",
  "urban_all_pairs",
  "urban_sharing_rate",
  "urban_sharing_ref_rate",
  
  # 3) Age-group counts (urban)
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
  
  # 4) Age-group rates (urban)
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
  
  # 5) Age-group all_pairs (urban)
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
  
  # 6) Urban distance summaries
  "urban_sharing_pairs_dist_median",
  "urban_all_pairs_dist_median",
  "urban_dist_p",
  
  # 7) Rural summary counts & rates
  "rural_sharing_pairs",
  "rural_ref_pairs",
  "rural_all_pairs",
  "rural_sharing_rate",
  
  # 8) Age-group counts (rural)
  "19 - 30 months / 19 - 30 months_r",
  "> 30 months / > 30 months_r",
  "> 30 months / 19 - 30 months_r",
  
  # 9) Age-group rates (rural)
  "19 - 30 months / 19 - 30 months_r_rate",
  "> 30 months / > 30 months_r_rate",
  "> 30 months / 19 - 30 months_r_rate",
  
  # 10) Age-group all_pairs (rural)
  "19 - 30 months / 19 - 30 months_r_all_pairs",
  "> 30 months / > 30 months_r_all_pairs",
  "> 30 months / 19 - 30 months_r_all_pairs",
  
  # 11) Rural distance summaries
  "rural_sharing_pairs_dist_median",
  "rural_all_pairs_dist_median",
  "rural_dist_p"
)

n=1
for (st in unique(all.ext.rename$V3)) {
  sub = all.ext.rename %>% filter(V3 == st) # This is StrainGST output
  sharing.sub = staingr.org.ML.name %>% filter(organism == st & HH_type == "Between")
  
  if(nrow(sharing.sub) == 0){ # pass if there isn't any strain-sharing associated with a strain
    next
  }
  ref.pair = combn(sub$V1, 2) %>% t() %>% as.data.frame()
  colnames(ref.pair) = c("sample1", "sample2")
  ref.pair = metadata_pair_load(ref.pair, R21_sample_list)

  prop.strain_sharing[n, "strain"] = st
  # Urban
  prop.strain_sharing[n, "urban_sharing_pairs"] = sharing.sub %>% filter(HH_type == "Between" & area == "urban") %>% nrow
  
  for(age in c("0 - 19 months / 0 - 19 months", "19 - 30 months / 19 - 30 months", "> 30 months / > 30 months", 
                 "Adults > 15 years / Adults > 15 years", "0 - 19 months / 19 - 30 months", "> 30 months / 0 - 19 months", "0 - 19 months / Adults > 15 years",
                 "> 30 months / 19 - 30 months", "19 - 30 months / Adults > 15 years", "> 30 months / Adults > 15 years")){
    prop.strain_sharing[n, age] = sharing.sub %>% filter(HH_type == "Between" & area == "urban" & Age_group_pair == age) %>% nrow
    prop.strain_sharing[n, paste0(age, "_rate")] = prop.strain_sharing[n, age]/urban.all.pairs %>% filter(HH_type == "Between" & Age_group_pair == age) %>% nrow
    prop.strain_sharing[n, paste0(age, "_all_pairs")] = urban.all.pairs %>% filter(HH_type == "Between" & Age_group_pair == age) %>% nrow
  }
  
  prop.strain_sharing[n, "urban_ref_pairs"] = ref.pair %>% filter(HH_type == "Between" & area == "urban") %>% nrow
  prop.strain_sharing[n, "urban_all_pairs"] = urban.all.pairs %>% filter(HH_type == "Between") %>% nrow
  prop.strain_sharing[n, "urban_sharing_ref_rate"] = prop.strain_sharing[n, "urban_sharing_pairs"] / prop.strain_sharing[n, "urban_ref_pairs"]
  prop.strain_sharing[n, "urban_sharing_rate"] = prop.strain_sharing[n, "urban_sharing_pairs"] / prop.strain_sharing[n, "urban_all_pairs"]
  # Get distances (Urban)
  urban.sharing.dist = sharing.sub %>% filter(HH_type == "Between" & area == "urban") %>% pull(distance)
  urban.ref.dist = ref.pair %>% filter(HH_type == "Between" & area == "urban") %>% pull(distance)
  urban.all.dist = urban.all.pairs %>% filter(HH_type == "Between") %>% pull(distance)
  # Add distance info
  if(length(urban.sharing.dist) > 0){
    prop.strain_sharing[n, "urban_sharing_pairs_dist_median"] = median(urban.sharing.dist)
    prop.strain_sharing[n, "urban_all_pairs_dist_median"] = median(urban.ref.dist)
    prop.strain_sharing[n, "urban_all_pairs_dist_median"] = median(urban.all.dist)
    # prop.strain_sharing[n, "urban_dist_p"] = wilcox.test(urban.sharing.dist, urban.ref.dist)$p.value
    prop.strain_sharing[n, "urban_dist_p"] = wilcox.test(urban.sharing.dist, urban.all.dist)$p.value
  }
  else{
    prop.strain_sharing[n, "urban_sharing_pairs_dist_median"] = NA
    prop.strain_sharing[n, "urban_all_pairs_dist_median"] = NA
    prop.strain_sharing[n, "urban_dist_p"] = NA
  }
  
  
  # Rural
  prop.strain_sharing[n, "rural_sharing_pairs"] = sharing.sub %>% filter(HH_type == "Between" & area == "rural") %>% nrow
  for(age in c("19 - 30 months / 19 - 30 months", "> 30 months / > 30 months" , "> 30 months / 19 - 30 months")){
    prop.strain_sharing[n, paste0(age, "_r")] = sharing.sub %>% filter(HH_type == "Between" & area == "rural" & Age_group_pair == age) %>% nrow
    prop.strain_sharing[n, paste0(age, "_r_rate")] = prop.strain_sharing[n, paste0(age, "_r")]/rural.all.pairs %>% filter(HH_type == "Between" & Age_group_pair == age) %>% nrow
    prop.strain_sharing[n, paste0(age, "_r_all_pairs")] = rural.all.pairs %>% filter(HH_type == "Between" & Age_group_pair == age) %>% nrow
  }
  prop.strain_sharing[n, "rural_ref_pairs"] =  ref.pair %>% filter(HH_type == "Between" & area == "rural") %>% nrow
  prop.strain_sharing[n, "rural_all_pairs"] = rural.all.pairs %>% filter(HH_type == "Between") %>% nrow
  prop.strain_sharing[n, "rural_sharing_rate"] = prop.strain_sharing[n, "rural_sharing_pairs"] / prop.strain_sharing[n, "rural_all_pairs"]
  # Get distances (rural)
  rural.sharing.dist = sharing.sub %>% filter(HH_type == "Between" & area == "rural") %>% pull(distance)
  # rural.ref.dist = ref.pair %>% filter(HH_type == "Between" & area == "rural") %>% pull(distance)
  rural.all.dist = rural.all.pairs %>% filter(HH_type == "Between") %>% pull(distance)
  # Add
  if(length(rural.sharing.dist) > 0){
    prop.strain_sharing[n, "rural_sharing_pairs_dist_median"] = median(rural.sharing.dist)
    # prop.strain_sharing[n, "rural_all_pairs_dist_median"] = median(rural.ref.dist)
    prop.strain_sharing[n, "rural_all_pairs_dist_median"] = median(rural.all.dist)
    # prop.strain_sharing[n, "rural_dist_p"] = wilcox.test(rural.sharing.dist, rural.ref.dist)$p.value
    prop.strain_sharing[n, "rural_dist_p"] = wilcox.test(rural.sharing.dist, rural.all.dist)$p.value
  }
  else{
    prop.strain_sharing[n, "rural_sharing_pairs_dist_median"] = NA
    prop.strain_sharing[n, "rural_all_pairs_dist_median"] = NA
    prop.strain_sharing[n, "rural_dist_p"] = NA
  }
  n=n+1
}

# Save into a file
write.csv(prop.strain_sharing, file.path(here::here(), "data/prop.strain_sharing.csv"), row.names = F, quote = F)
saveRDS(prop.strain_sharing, file.path(here::here(), "data/prop.strain_sharing.rds"))
prop.strain_sharing = readRDS(file.path(here::here(), "data/prop.strain_sharing.rds"))


prop.strain_sharing.urban = prop.strain_sharing %>% select(c("strain", "urban_sharing_pairs", "0 - 19 months / 0 - 19 months", "19 - 30 months / 19 - 30 months", "> 30 months / > 30 months", 
                                                             "Adults > 15 years / Adults > 15 years", "0 - 19 months / 19 - 30 months", "> 30 months / 0 - 19 months", "0 - 19 months / Adults > 15 years",
                                                             "> 30 months / 19 - 30 months", "19 - 30 months / Adults > 15 years", "> 30 months / Adults > 15 years",
                                                             "0 - 19 months / 0 - 19 months_rate", "19 - 30 months / 19 - 30 months_rate", "> 30 months / > 30 months_rate", 
                                                             "Adults > 15 years / Adults > 15 years_rate", "0 - 19 months / 19 - 30 months_rate", "> 30 months / 0 - 19 months_rate", "0 - 19 months / Adults > 15 years_rate",
                                                             "> 30 months / 19 - 30 months_rate", "19 - 30 months / Adults > 15 years_rate", "> 30 months / Adults > 15 years_rate",
                                                             "urban_all_pairs", "urban_sharing_rate","urban_sharing_pairs_dist_median", "urban_all_pairs_dist_median", "urban_dist_p"))  %>% filter(!is.na(urban_sharing_pairs_dist_median))

prop.strain_sharing.urban.m = prop.strain_sharing.urban %>% select("strain", "0 - 19 months / 0 - 19 months_rate", "19 - 30 months / 19 - 30 months_rate", "> 30 months / > 30 months_rate", 
                                                                   "Adults > 15 years / Adults > 15 years_rate", "0 - 19 months / 19 - 30 months_rate", "> 30 months / 0 - 19 months_rate", "0 - 19 months / Adults > 15 years_rate",
                                                                   "> 30 months / 19 - 30 months_rate", "19 - 30 months / Adults > 15 years_rate", "> 30 months / Adults > 15 years_rate") %>% melt()

prop.strain_sharing.rural = prop.strain_sharing %>% select(c("strain", "rural_sharing_pairs", "19 - 30 months / 19 - 30 months_r", "> 30 months / > 30 months_r" , "> 30 months / 19 - 30 months_r",
                                                             "19 - 30 months / 19 - 30 months_r_rate", "> 30 months / > 30 months_r_rate" , "> 30 months / 19 - 30 months_r_rate",
                                                             "rural_all_pairs", "rural_sharing_rate", "rural_sharing_pairs_dist_median", "rural_all_pairs_dist_median", "rural_dist_p")) %>% filter(!is.na(rural_sharing_pairs_dist_median))



prop.strain_sharing.rural.m = prop.strain_sharing.rural %>% select("strain", "19 - 30 months / 19 - 30 months_r_rate", "> 30 months / > 30 months_r_rate" , "> 30 months / 19 - 30 months_r_rate") %>% melt()

prop.strain_sharing_within = data.frame(matrix(ncol = 40))
colnames(prop.strain_sharing_within) <- c(
  # 1) Strain identifier
  "strain",
  
  # 2) Urban summary counts & rates
  "urban_sharing_pairs",
  "urban_ref_pairs",
  "urban_all_pairs",
  "urban_sharing_rate",
  "urban_sharing_ref_rate",
  
  # 3) Age-group counts (urban)
  "0 - 19 months / 0 - 19 months",
  "> 30 months / > 30 months",
  "0 - 19 months / 19 - 30 months",
  "> 30 months / 0 - 19 months",
  "0 - 19 months / Adults > 15 years",
  "> 30 months / 19 - 30 months",
  "19 - 30 months / Adults > 15 years",
  "> 30 months / Adults > 15 years",
  
  # 4) Age-group rates (urban)
  "0 - 19 months / 0 - 19 months_rate",
  "> 30 months / > 30 months_rate",
  "0 - 19 months / 19 - 30 months_rate",
  "> 30 months / 0 - 19 months_rate",
  "0 - 19 months / Adults > 15 years_rate",
  "> 30 months / 19 - 30 months_rate",
  "19 - 30 months / Adults > 15 years_rate",
  "> 30 months / Adults > 15 years_rate",
  
  # 5) Age-group all_pairs (urban)
  "0 - 19 months / 0 - 19 months_all_pairs",
  "> 30 months / > 30 months_all_pairs",
  "0 - 19 months / 19 - 30 months_all_pairs",
  "> 30 months / 0 - 19 months_all_pairs",
  "0 - 19 months / Adults > 15 years_all_pairs",
  "> 30 months / 19 - 30 months_all_pairs",
  "19 - 30 months / Adults > 15 years_all_pairs",
  "> 30 months / Adults > 15 years_all_pairs",
  
  # 7) Rural summary counts & rates
  "rural_sharing_pairs",
  "rural_ref_pairs",
  "rural_all_pairs",
  "rural_sharing_rate",
  
  # 8) Age-group counts (rural)
  "19 - 30 months / 19 - 30 months_r",
  
  "> 30 months / 19 - 30 months_r",
  
  # 9) Age-group rates (rural)
  "19 - 30 months / 19 - 30 months_r_rate",
  
  "> 30 months / 19 - 30 months_r_rate",
  
  # 10) Age-group all_pairs (rural)
  "19 - 30 months / 19 - 30 months_r_all_pairs",
  
  "> 30 months / 19 - 30 months_r_all_pairs"
  
)

n=1
for (st in unique(all.ext.rename$V3)) {
  sub = all.ext.rename %>% filter(V3 == st) # This is StrainGST output
  sharing.sub = staingr.org.ML.name %>% filter(organism == st & HH_type == "Within")
  
  if(nrow(sharing.sub) == 0){ # pass if there isn't any strain-sharing associated with a strain
    next
  }
  ref.pair = combn(sub$V1, 2) %>% t() %>% as.data.frame()
  colnames(ref.pair) = c("sample1", "sample2")
  ref.pair = metadata_pair_load(ref.pair, R21_sample_list)
  
  prop.strain_sharing_within[n, "strain"] = st
  # Urban
  prop.strain_sharing_within[n, "urban_sharing_pairs"] = sharing.sub %>% filter(HH_type == "Within" & area == "urban") %>% nrow
  #c("0 - 19 months / 0 - 19 months", "19 - 30 months / 19 - 30 months", "> 30 months / > 30 months", 
  #  "Adults > 15 years / Adults > 15 years", "0 - 19 months / 19 - 30 months", "> 30 months / 0 - 19 months", "0 - 19 months / Adults > 15 years",
  #  "> 30 months / 19 - 30 months", "19 - 30 months / Adults > 15 years", "> 30 months / Adults > 15 years")
  for(age in urban.all.pairs %>% filter(HH_type == "Within") %>% pull(Age_group_pair) %>% unique){
    prop.strain_sharing_within[n, age] = sharing.sub %>% filter(HH_type == "Within" & area == "urban" & Age_group_pair == age) %>% nrow
    prop.strain_sharing_within[n, paste0(age, "_rate")] = prop.strain_sharing_within[n, age]/urban.all.pairs %>% filter(HH_type == "Within" & Age_group_pair == age) %>% nrow
    prop.strain_sharing_within[n, paste0(age, "_all_pairs")] = urban.all.pairs %>% filter(HH_type == "Within" & Age_group_pair == age) %>% nrow
  }
  
  prop.strain_sharing_within[n, "urban_ref_pairs"] = ref.pair %>% filter(HH_type == "Within" & area == "urban") %>% nrow
  prop.strain_sharing_within[n, "urban_all_pairs"] = urban.all.pairs %>% filter(HH_type == "Within") %>% nrow
  prop.strain_sharing_within[n, "urban_sharing_ref_rate"] = prop.strain_sharing_within[n, "urban_sharing_pairs"] / prop.strain_sharing_within[n, "urban_ref_pairs"]
  prop.strain_sharing_within[n, "urban_sharing_rate"] = prop.strain_sharing_within[n, "urban_sharing_pairs"] / prop.strain_sharing_within[n, "urban_all_pairs"]  
  
  # Rural
  prop.strain_sharing_within[n, "rural_sharing_pairs"] = sharing.sub %>% filter(HH_type == "Within" & area == "rural") %>% nrow
  for(age in rural.all.pairs %>% filter(HH_type == "Within") %>% pull(Age_group_pair) %>% unique){ # c("19 - 30 months / 19 - 30 months", "> 30 months / > 30 months" , "> 30 months / 19 - 30 months")
    prop.strain_sharing_within[n, paste0(age, "_r")] = sharing.sub %>% filter(HH_type == "Within" & area == "rural" & Age_group_pair == age) %>% nrow
    prop.strain_sharing_within[n, paste0(age, "_r_rate")] = prop.strain_sharing_within[n, paste0(age, "_r")]/rural.all.pairs %>% filter(HH_type == "Within" & Age_group_pair == age) %>% nrow
    prop.strain_sharing_within[n, paste0(age, "_r_all_pairs")] = rural.all.pairs %>% filter(HH_type == "Within" & Age_group_pair == age) %>% nrow
  }
  prop.strain_sharing_within[n, "rural_ref_pairs"] =  ref.pair %>% filter(HH_type == "Within" & area == "rural") %>% nrow
  prop.strain_sharing_within[n, "rural_all_pairs"] = rural.all.pairs %>% filter(HH_type == "Within") %>% nrow
  prop.strain_sharing_within[n, "rural_sharing_rate"] = prop.strain_sharing_within[n, "rural_sharing_pairs"] / prop.strain_sharing_within[n, "rural_all_pairs"]

  n=n+1
}

urban.all.pairs %>% filter(HH_type == "Within") %>% select(Age_group_pair) %>% unique
urban.all.pairs %>% filter(HH_type == "Within") %>% select(Age_group_pair) %>% table

rural.all.pairs %>% filter(HH_type == "Within") %>% select(Age_group_pair) %>% unique
rural.all.pairs %>% filter(HH_type == "Within") %>% select(Age_group_pair) %>% table

# Save into a file
write.csv(prop.strain_sharing_within, file.path(here::here(), "data/prop.strain_sharing_within.csv"), row.names = F, quote = F)
saveRDS(prop.strain_sharing_within, file.path(here::here(), "data/prop.strain_sharing_within.rds"))
