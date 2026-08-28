# Odds ratio between household distance and treatment type on strain sharing and strain sharing between age groups
source(file.path(here::here(),"0-config.R"))
#Load data
load(file.path(here::here(), "data/all.ML.networks.RData"))
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
  mutate(Age_group_pair = paste(sort(c(sample1_Age_group, sample2_Age_group)), collapse = "—")) %>%
  ungroup()

rural.all.pairs <- rural.all.pairs %>%
  rowwise() %>%
  mutate(Age_group_pair = paste(sort(c(sample1_Age_group, sample2_Age_group)), collapse = "—")) %>%
  ungroup()

urban.all.pairs.tran = urban.all.pairs %>% filter(HH_type == "Between" & tr_type != "Different") %>% mutate(all.log = ifelse(All.share == "Yes", 1, 0), comm.log = ifelse(commensal.share == "Yes", 1, 0), noncomm.log = ifelse(non_commensal.share == "Yes", 1, 0))
rural.all.pairs.tran = rural.all.pairs %>% filter(HH_type == "Between" & tr_type != "Different") %>% mutate(all.log = ifelse(All.share == "Yes", 1, 0), comm.log = ifelse(commensal.share == "Yes", 1, 0), noncomm.log = ifelse(non_commensal.share == "Yes", 1, 0))
# scailing distances
urban.all.pairs.tran$dist_tran = urban.all.pairs.tran$distance %>% scale()
rural.all.pairs.tran$dist_tran = rural.all.pairs.tran$distance %>% scale()

# Urban
model1 <- glm(all.log ~  tr_type + dist_tran, data = urban.all.pairs.tran, family = binomial)
model2 <- glm(comm.log ~  tr_type + dist_tran, data = urban.all.pairs.tran, family = binomial)
model3 <- glm(noncomm.log ~  tr_type + dist_tran, data = urban.all.pairs.tran, family = binomial)
# Rural
model4 <- glm(all.log ~ tr_type + dist_tran, data = rural.all.pairs.tran, family = binomial)
model5 <- glm(comm.log ~ tr_type + dist_tran, data = rural.all.pairs.tran, family = binomial)
model6 <- glm(noncomm.log ~ tr_type + dist_tran, data = rural.all.pairs.tran, family = binomial)

get_Odds_ratio_conf = function(model, a1 = "urban", t1 = "All"){
  # Get confidence intervals
  confint_model <- confint(model)
  
  # Exponentiate the coefficients and confidence intervals
  odds_ratios <- exp(coef(model))
  conf_intervals <- exp(confint_model)
  
  # Combine the odds ratios and confidence intervals in a table
  results <- cbind(odds_ratios, conf_intervals)  %>% as.data.frame()
  colnames(results) <- c("Odds_Ratio", "Lower_CI", "Upper_CI")
  results$var = rownames(results)
  results = results[rownames(results) != "(Intercept)",]
  results$area = a1
  results$target= t1
  rownames(results) <- NULL # remove row names
  return(results)
}

# Get OR for all data
urban.all.OR = get_Odds_ratio_conf(model = model1, a1 = "urban", t1 = "All")
urban.comm.OR = get_Odds_ratio_conf(model = model2, a1 = "urban", t1 = "comm")
urban.noncomm.OR = get_Odds_ratio_conf(model = model3, a1 = "urban", t1 = "noncomm")

rural.all.OR = get_Odds_ratio_conf(model = model4, a1 = "rural", t1 = "All")
rural.comm.OR = get_Odds_ratio_conf(model = model5, a1 = "rural", t1 = "comm")
rural.noncomm.OR = get_Odds_ratio_conf(model = model6, a1 = "rural", t1 = "noncomm")
# Combine data
all.OR = rbind(urban.all.OR, urban.comm.OR, urban.noncomm.OR, rural.all.OR, rural.comm.OR, rural.noncomm.OR)
all.OR$var = factor(all.OR$var, levels = c("dist_tran", "tr_typeWater"), labels = c("Household distance", "Water treatment"))

# Save into R data
saveRDS(all.OR, file = file.path(here::here(), "data/all.OR.rds"))

# Genus shared more among age groups (Does it differ by treatment status?)
# Urban
# Overall
urban.all.pairs.bet.age.mat = urban.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, All.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "All") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.comm.pairs.bet.age.mat = urban.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, commensal.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "commensals") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.noncomm.pairs.bet.age.mat = urban.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, non_commensal.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "non_commensals") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Bifi.pairs.bet.age.mat = urban.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, Bifi.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Bifi") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Bact.pairs.bet.age.mat = urban.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, Bact.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Bact") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Esch.pairs.bet.age.mat = urban.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, Esch.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Esch") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Camp_MAG.pairs.bet.age.mat = urban.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, Camp_MAG.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Camp_MAG") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Entc.pairs.bet.age.mat = urban.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, Entc.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Entc") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Entb.pairs.bet.age.mat = urban.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, Entb.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Entb") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Kleb.pairs.bet.age.mat = urban.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, Kleb.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Kleb") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
# Water
urban.all.pairs.bet.age.mat.wat = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, All.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "All") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.comm.pairs.bet.age.mat.wat = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, commensal.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "commensals") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.noncomm.pairs.bet.age.mat.wat = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, non_commensal.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "non_commensals") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Bifi.pairs.bet.age.mat.wat = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, Bifi.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Bifi") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Bact.pairs.bet.age.mat.wat = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, Bact.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Bact") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Esch.pairs.bet.age.mat.wat = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, Esch.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Esch") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Camp_MAG.pairs.bet.age.mat.wat = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, Camp_MAG.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Camp_MAG") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Entc.pairs.bet.age.mat.wat = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, Entc.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Entc") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Entb.pairs.bet.age.mat.wat = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, Entb.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Entb") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—") %>% mutate(Yes = 0) # Adding columns for "Yes" because of no case
urban.Kleb.pairs.bet.age.mat.wat = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, Kleb.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Kleb") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—") %>% mutate(Yes = 0) # Adding columns for "Yes" because of no case
# Control
urban.all.pairs.bet.age.mat.con = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, All.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "All") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.comm.pairs.bet.age.mat.con = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, commensal.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "commensals") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.noncomm.pairs.bet.age.mat.con = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, non_commensal.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "non_commensals") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Bifi.pairs.bet.age.mat.con = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, Bifi.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Bifi") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Bact.pairs.bet.age.mat.con = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, Bact.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Bact") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Esch.pairs.bet.age.mat.con = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, Esch.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Esch") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Camp_MAG.pairs.bet.age.mat.con = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, Camp_MAG.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Camp_MAG") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Entc.pairs.bet.age.mat.con = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, Entc.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Entc") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
urban.Entb.pairs.bet.age.mat.con = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, Entb.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Entb") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—") %>% mutate(Yes = 0) # Adding columns for "Yes" because of no case
urban.Kleb.pairs.bet.age.mat.con = urban.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, Kleb.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Kleb") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—") %>% mutate(Yes = 0) # Adding columns for "Yes" because of no case

urban.pairs.bet.age.mat = rbind(urban.all.pairs.bet.age.mat, urban.comm.pairs.bet.age.mat, urban.noncomm.pairs.bet.age.mat, 
                                urban.Bifi.pairs.bet.age.mat, urban.Bact.pairs.bet.age.mat, urban.Esch.pairs.bet.age.mat, urban.Camp_MAG.pairs.bet.age.mat,
                                urban.Entc.pairs.bet.age.mat, urban.Entb.pairs.bet.age.mat, urban.Kleb.pairs.bet.age.mat)
urban.pairs.bet.age.mat$tr = "Overall"

urban.pairs.bet.age.mat.wat = rbind(urban.all.pairs.bet.age.mat.wat, urban.comm.pairs.bet.age.mat.wat, urban.noncomm.pairs.bet.age.mat.wat, 
                                urban.Bifi.pairs.bet.age.mat.wat, urban.Bact.pairs.bet.age.mat.wat, urban.Esch.pairs.bet.age.mat.wat, urban.Camp_MAG.pairs.bet.age.mat.wat,
                                urban.Entc.pairs.bet.age.mat.wat, urban.Entb.pairs.bet.age.mat.wat, urban.Kleb.pairs.bet.age.mat.wat)
urban.pairs.bet.age.mat.wat$tr = "Water"

urban.pairs.bet.age.mat.con = rbind(urban.all.pairs.bet.age.mat.con, urban.comm.pairs.bet.age.mat.con, urban.noncomm.pairs.bet.age.mat.con, 
                                    urban.Bifi.pairs.bet.age.mat.con, urban.Bact.pairs.bet.age.mat.con, urban.Esch.pairs.bet.age.mat.con, urban.Camp_MAG.pairs.bet.age.mat.con,
                                    urban.Entc.pairs.bet.age.mat.con, urban.Entb.pairs.bet.age.mat.con, urban.Kleb.pairs.bet.age.mat.con)
urban.pairs.bet.age.mat.con$tr = "Control"

urban.pairs.bet.age.mat = rbind(urban.pairs.bet.age.mat,urban.pairs.bet.age.mat.wat, urban.pairs.bet.age.mat.con )

urban.pairs.bet.age.mat$percent = urban.pairs.bet.age.mat$Yes / (urban.pairs.bet.age.mat$Yes + urban.pairs.bet.age.mat$No) *100
urban.pairs.bet.age.mat$cat = paste0(urban.pairs.bet.age.mat$Group1, " — ", urban.pairs.bet.age.mat$Group2)
urban.pairs.bet.age.mat$cat = factor(urban.pairs.bet.age.mat$cat, levels = c("0 - 19 months — 0 - 19 months", "19 - 30 months — 19 - 30 months", "> 30 months — > 30 months",
                                                                             "Adults > 15 years — Adults > 15 years", "0 - 19 months — 19 - 30 months", "0 - 19 months — > 30 months", "0 - 19 months — Adults > 15 years",
                                                                             "19 - 30 months — > 30 months", "19 - 30 months — Adults > 15 years", "> 30 months — Adults > 15 years"))
urban.all.pairs.bet.age.mat1 = urban.pairs.bet.age.mat %>% filter(tr == "Overall") %>% select(-c(Yes, No, cat))

# Rural
# Overall
rural.all.pairs.bet.age.mat = rural.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, All.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "All") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.comm.pairs.bet.age.mat = rural.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, commensal.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "commensals") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.noncomm.pairs.bet.age.mat = rural.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, non_commensal.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "non_commensals") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Bifi.pairs.bet.age.mat = rural.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, Bifi.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Bifi") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Bact.pairs.bet.age.mat = rural.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, Bact.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Bact") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Esch.pairs.bet.age.mat = rural.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, Esch.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Esch") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Camp_MAG.pairs.bet.age.mat = rural.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, Camp_MAG.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Camp_MAG") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Entc.pairs.bet.age.mat = rural.all.pairs %>% filter(HH_type == "Between") %>% select(Age_group_pair, Entc.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Entc") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")

# Water
rural.all.pairs.bet.age.mat.wat = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, All.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "All") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.comm.pairs.bet.age.mat.wat = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, commensal.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "commensals") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.noncomm.pairs.bet.age.mat.wat = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, non_commensal.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "non_commensals") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Bifi.pairs.bet.age.mat.wat = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, Bifi.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Bifi") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Bact.pairs.bet.age.mat.wat = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, Bact.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Bact") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Esch.pairs.bet.age.mat.wat = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, Esch.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Esch") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Camp_MAG.pairs.bet.age.mat.wat = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, Camp_MAG.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Camp_MAG") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Entc.pairs.bet.age.mat.wat = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Water") %>% select(Age_group_pair, Entc.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Entc") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")

# Control
rural.all.pairs.bet.age.mat.con = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, All.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "All") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.comm.pairs.bet.age.mat.con = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, commensal.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "commensals") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.noncomm.pairs.bet.age.mat.con = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, non_commensal.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "non_commensals") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Bifi.pairs.bet.age.mat.con = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, Bifi.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Bifi") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Bact.pairs.bet.age.mat.con = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, Bact.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Bact") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Esch.pairs.bet.age.mat.con = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, Esch.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Esch") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Camp_MAG.pairs.bet.age.mat.con = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, Camp_MAG.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Camp_MAG") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")
rural.Entc.pairs.bet.age.mat.con = rural.all.pairs %>% filter(HH_type == "Between" & tr_type == "Control") %>% select(Age_group_pair, Entc.share) %>% table() %>% as.data.frame.matrix() %>% mutate(genus = "Entc") %>% rownames_to_column(var = "Age_group_pair") %>% separate(Age_group_pair, into = c("Group1", "Group2"), sep = "—")

rural.pairs.bet.age.mat = rbind(rural.all.pairs.bet.age.mat, rural.comm.pairs.bet.age.mat, rural.noncomm.pairs.bet.age.mat, 
                                rural.Bifi.pairs.bet.age.mat, rural.Bact.pairs.bet.age.mat, rural.Esch.pairs.bet.age.mat, rural.Camp_MAG.pairs.bet.age.mat,
                                rural.Entc.pairs.bet.age.mat)
rural.pairs.bet.age.mat$tr = "Overall"

rural.pairs.bet.age.mat.wat = rbind(rural.all.pairs.bet.age.mat.wat, rural.comm.pairs.bet.age.mat.wat, rural.noncomm.pairs.bet.age.mat.wat, 
                                    rural.Bifi.pairs.bet.age.mat.wat, rural.Bact.pairs.bet.age.mat.wat, rural.Esch.pairs.bet.age.mat.wat, rural.Camp_MAG.pairs.bet.age.mat.wat,
                                    rural.Entc.pairs.bet.age.mat.wat)
rural.pairs.bet.age.mat.wat$tr = "Water"

rural.pairs.bet.age.mat.con = rbind(rural.all.pairs.bet.age.mat.con, rural.comm.pairs.bet.age.mat.con, rural.noncomm.pairs.bet.age.mat.con, 
                                    rural.Bifi.pairs.bet.age.mat.con, rural.Bact.pairs.bet.age.mat.con, rural.Esch.pairs.bet.age.mat.con, rural.Camp_MAG.pairs.bet.age.mat.con,
                                    rural.Entc.pairs.bet.age.mat.con)
rural.pairs.bet.age.mat.con$tr = "Control"

rural.pairs.bet.age.mat = rbind(rural.pairs.bet.age.mat,rural.pairs.bet.age.mat.wat, rural.pairs.bet.age.mat.con)

rural.pairs.bet.age.mat$percent = rural.pairs.bet.age.mat$Yes / (rural.pairs.bet.age.mat$Yes + rural.pairs.bet.age.mat$No) *100
rural.pairs.bet.age.mat$cat = paste0(rural.pairs.bet.age.mat$Group1, " — ", rural.pairs.bet.age.mat$Group2)
rural.pairs.bet.age.mat$cat = factor(rural.pairs.bet.age.mat$cat, levels = c("0 - 19 months — 0 - 19 months", "19 - 30 months — 19 - 30 months", "> 30 months — > 30 months",
                                                                             "Adults > 15 years — Adults > 15 years", "0 - 19 months — 19 - 30 months", "0 - 19 months — > 30 months", "0 - 19 months — Adults > 15 years",
                                                                             "19 - 30 months — > 30 months", "19 - 30 months — Adults > 15 years", "> 30 months — Adults > 15 years"))
rural.all.pairs.bet.age.mat1 = rural.pairs.bet.age.mat %>% select(-c(Yes, No, cat))

# Save into R data
save(urban.pairs.bet.age.mat, urban.all.pairs.bet.age.mat1, rural.pairs.bet.age.mat, rural.all.pairs.bet.age.mat1, file = file.path(here::here(), "data/all.pairs.bet.age.mat.RData"))

# Calculate Fisher's test p value (water vs. control by age groups)
urban.pairs.bet.age.mat.stat = calculate_percent_fisher(urban.pairs.bet.age.mat)
rural.pairs.bet.age.mat.stat = calculate_percent_fisher(rural.pairs.bet.age.mat)
# Save into files
write.table(urban.pairs.bet.age.mat.stat, file.path(here::here(), "data/urban.pairs.bet.age.mat.stat.tsv"), sep = "\t", row.names = F, quote = F)
write.table(rural.pairs.bet.age.mat.stat, file.path(here::here(), "data/rural.pairs.bet.age.mat.stat.tsv"), sep = "\t", row.names = F, quote = F)

urban.pairs.bet.age.mat.stat = read.table(file.path(here::here(), "data/urban.pairs.bet.age.mat.stat.tsv"), sep = "\t", header = T)
rural.pairs.bet.age.mat.stat = read.table(file.path(here::here(), "data/rural.pairs.bet.age.mat.stat.tsv"), sep = "\t", header = T)

urban.pairs.bet.age.mat.stat %>% filter(p_value < 0.05)
rural.pairs.bet.age.mat.stat %>% filter(p_value < 0.05)

##########################################################################################
# Grouping 0 - 19 months and 19- 30 months to 0 - 30 months in urban area (water vs. control)
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

urban.pairs.bet.age.mat.g.stat = calculate_percent_fisher(urban.pairs.bet.age.mat.g1)
urban.pairs.bet.age.mat.g.stat %>% filter(p_value < 0.05)

write.table(urban.pairs.bet.age.mat.g.stat, file.path(here::here(), "data/urban.pairs.bet.age.mat.g.stat.tsv"), sep = "\t", row.names = F, quote = F)

# Comparison of within HH strain-sharing in Urban area
urban.all.pairs.wit = urban.all.pairs %>% filter(HH_type == "Within") %>% rowwise() %>%
  mutate(within_cat = ifelse("Adults > 15 years" %in% c(sample1_Age_group, sample2_Age_group), "Mother - Child", "Child - Child"))

urban.all.pairs.wit$hhid <- sub("^\\D*(\\d+)\\d{2}(_.*)?$", "\\1", urban.all.pairs.wit$sample1)

# Mother - Child vs. Child - Child
urban.all.pairs.wit.all.tab = urban.all.pairs.wit %>% select(within_cat, All.share) %>% table %>% as.data.frame.matrix() %>% rownames_to_column(var = "Age_group_pair") %>% rowwise() %>% mutate(percent = Yes / (Yes + No)*100) %>% mutate(genus = "all", tr = "Overall")
urban.all.pairs.wit.comm.tab = urban.all.pairs.wit %>% select(within_cat, commensal.share) %>% table %>% as.data.frame.matrix() %>% rownames_to_column(var = "Age_group_pair") %>% rowwise() %>% mutate(percent = Yes / (Yes + No)*100) %>% mutate(genus = "commensals", tr = "Overall")
urban.all.pairs.wit.non_comm.tab = urban.all.pairs.wit %>% select(within_cat, non_commensal.share) %>% table %>% as.data.frame.matrix() %>% rownames_to_column(var = "Age_group_pair") %>% rowwise() %>% mutate(percent = Yes / (Yes + No)*100) %>% mutate(genus = "non_commensals", tr = "Overall")

urban.all.pairs.wit.all.tab.w = urban.all.pairs.wit %>% filter(tr_type == "Water") %>% select(within_cat, All.share) %>% table %>% as.data.frame.matrix() %>% rownames_to_column(var = "Age_group_pair") %>% rowwise() %>% mutate(percent = Yes / (Yes + No)*100) %>% mutate(genus = "all", tr = "Water")
urban.all.pairs.wit.comm.tab.w = urban.all.pairs.wit %>% filter(tr_type == "Water") %>% select(within_cat, commensal.share) %>% table %>% as.data.frame.matrix() %>% rownames_to_column(var = "Age_group_pair") %>% rowwise() %>% mutate(percent = Yes / (Yes + No)*100) %>% mutate(genus = "commensals", tr = "Water")
urban.all.pairs.wit.non_comm.tab.w = urban.all.pairs.wit %>% filter(tr_type == "Water") %>% select(within_cat, non_commensal.share) %>% table %>% as.data.frame.matrix() %>% rownames_to_column(var = "Age_group_pair") %>% rowwise() %>% mutate(percent = Yes / (Yes + No)*100) %>% mutate(genus = "non_commensals", tr = "Water")

urban.all.pairs.wit.all.tab.c = urban.all.pairs.wit %>% filter(tr_type == "Control") %>% select(within_cat, All.share) %>% table %>% as.data.frame.matrix() %>% rownames_to_column(var = "Age_group_pair") %>% rowwise() %>% mutate(percent = Yes / (Yes + No)*100) %>% mutate(genus = "all", tr = "Control")
urban.all.pairs.wit.comm.tab.c = urban.all.pairs.wit %>% filter(tr_type == "Control") %>% select(within_cat, commensal.share) %>% table %>% as.data.frame.matrix() %>% rownames_to_column(var = "Age_group_pair") %>% rowwise() %>% mutate(percent = Yes / (Yes + No)*100) %>% mutate(genus = "commensals", tr = "Control")
urban.all.pairs.wit.non_comm.tab.c = urban.all.pairs.wit %>% filter(tr_type == "Control") %>% select(within_cat, non_commensal.share) %>% table %>% as.data.frame.matrix() %>% rownames_to_column(var = "Age_group_pair") %>% rowwise() %>% mutate(percent = Yes / (Yes + No)*100) %>% mutate(genus = "non_commensals", tr = "Control")

urban.all.pairs.wit.tab = rbind(urban.all.pairs.wit.all.tab, urban.all.pairs.wit.comm.tab, urban.all.pairs.wit.non_comm.tab,
                                urban.all.pairs.wit.all.tab.w, urban.all.pairs.wit.comm.tab.w, urban.all.pairs.wit.non_comm.tab.w,
                                urban.all.pairs.wit.all.tab.c, urban.all.pairs.wit.comm.tab.c, urban.all.pairs.wit.non_comm.tab.c )
urban.all.pairs.wit.tab$tr = factor(urban.all.pairs.wit.tab$tr, levels = c("Overall", "Water", "Control"))

# Save into a file
write.csv(urban.all.pairs.wit.tab, file.path(here::here(), "data/stats/urban.all.pairs.wit.tab.csv"), row.names = F, quote = F)
