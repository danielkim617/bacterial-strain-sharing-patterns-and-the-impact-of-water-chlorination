# function for loading StrainGR output
get.pair = function(x){
  colnames(x) = c('sample', 'i', 'strain',  'gkmers', 'ikmers', 'skmers', 'cov',  'kcov', 'gcov', 'acct', 'even', 'spec', 'rapct', 'old_rapct', 'wscore', 'score')
  rural.samples = R21_sample_list %>% filter(area == "rural") %>% pull(sample_id) %>% unique()
  urban.samples = R21_sample_list %>% filter(area == "urban") %>% pull(sample_id) %>% unique()
  #Get pair of samples with the same strain identified
  strains.pair = data.frame(matrix(ncol=3))
  colnames(strains.pair) = c("strain", "sample1", "sample2")
  for (st in x$strain %>% unique()) {
    a = x %>% filter(strain == st)
    if(length(a$sample) > 1){
      b = combn(a$sample,2) %>% t() %>% as.data.frame()
      colnames(b) = c("sample1", "sample2")
      b$strain = st
      strains.pair = data.frame(rbind(strains.pair, b))
    }
  }
  strains.pair = strains.pair %>% filter(!is.na(strains.pair$strain))
  # Sort each row
  strains.pair_sorted <- t(apply(strains.pair[,2:3], 1, sort))
  # Convert back to a data frame
  strains.pair_sorted <- as.data.frame(strains.pair_sorted)
  # Remove duplicate rows
  strains.pair_sorted <- unique(strains.pair_sorted)
  ###Removing sample pairs between rural and urban area to decrease the number of comparisons for commensal organisms
  strains.pair_sorted.r = strains.pair_sorted %>% filter(sample1 %in% rural.samples & sample2 %in% rural.samples)
  strains.pair_sorted.u = strains.pair_sorted %>% filter(sample1 %in% urban.samples & sample2 %in% urban.samples)
  strains.pair_sorted = rbind(strains.pair_sorted.r, strains.pair_sorted.u)
  ###
  return(strains.pair_sorted)
}

#function for calculating strain frequency
get_strain_freq = function(x, y = R21_sample_list){
  #children = x %>% filter(!endsWith(sample,"150")) # filter out adult samples (> 15 yrs)
  colnames(x) = c('sample', 'i', 'strain',  'gkmers', 'ikmers', 'skmers', 'cov',  'kcov', 'gcov', 'acct', 'even', 'spec', 'rapct', 'old_rapct', 'wscore', 'score')
  x = merge(x, y, by.x = "sample", by.y = "sample_id", all.x =T) # Add sample information
  
  freq.table = data.frame(matrix(ncol=8))
  colnames(freq.table) = c("strain", "all", "urban", "urban_water", "urban_control", "rural", "rural_water", "rural_control")
  n=1
  for (st in x$strain %>% unique()) {
    freq.table[n,1] = st
    freq.table[n,2] = x %>% filter(strain==st) %>% pull(sample) %>% unique %>% length # All
    freq.table[n,3] = x %>% filter(strain==st & area=="urban") %>% pull(sample) %>% unique %>% length # Urban
    freq.table[n,4] = x %>% filter(strain==st & area=="urban" & tr=="Water") %>% pull(sample) %>% unique %>% length # Urban_water
    freq.table[n,5] = x %>% filter(strain==st & area=="urban" & tr=="Control") %>% pull(sample) %>% unique %>% length # Urban_control
    freq.table[n,6] = x %>% filter(strain==st & area=="rural") %>% pull(sample) %>% unique %>% length # Urban
    freq.table[n,7] = x %>% filter(strain==st & area=="rural" & tr=="Water") %>% pull(sample) %>% unique %>% length # Urban_water
    freq.table[n,8] = x %>% filter(strain==st & area=="rural" & tr=="Control") %>% pull(sample) %>% unique %>% length # Urban_control
    n=n+1
  }
  freq.table$strain = gsub(".fa.gz","", freq.table$strain)
  #Get relative frequency
  freq.table[,2] = freq.table[,2]/length(y %>% pull(sample_id) %>% unique) *100
  freq.table[,3] = freq.table[,3]/length(y %>% filter(area=="urban") %>% pull(sample_id) %>% unique)*100
  freq.table[,4] = freq.table[,4]/length(y %>% filter(area=="urban" & tr=="Water") %>% pull(sample_id) %>% unique)*100
  freq.table[,5] = freq.table[,5]/length(y %>% filter(area=="urban" & tr=="Control") %>% pull(sample_id) %>% unique)*100
  freq.table[,6] = freq.table[,6]/length(y %>% filter(area=="rural") %>% pull(sample_id) %>% unique)*100
  freq.table[,7] = freq.table[,7]/length(y %>% filter(area=="rural" & tr=="Water") %>% pull(sample_id) %>% unique)*100
  freq.table[,8] = freq.table[,8]/length(y %>% filter(area=="rural" & tr=="Control") %>% pull(sample_id) %>% unique)*100
  return(freq.table)
}

make_split_heatmap <- function(x, split_n = 4, group_vars = c("urban", "rural")) {
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(RColorBrewer)
  library(patchwork)
  
  # Melt and clean data
  x.m <- melt(x)
  x.m[x.m$value == 0 , "value"] <- NA  # replace 0 with NA
  x.m <- x.m %>% filter(variable %in% group_vars)
  # Add factor level
  x.m$variable = factor(x.m$variable, levels = c("rural", "urban"), labels = c("Western Kenya", "Nairobi"))
  
  x.m$strain <- gsub("Esch_", "Escherichia_", x.m$strain)
  x.m$strain <- gsub("Camp_", "Campylobacter_", x.m$strain) 
  x.m$strain <- gsub("Kleb_", "Klebsiella_", x.m$strain) 
  x.m$strain <- gsub("Bact_", "Bacteroides_", x.m$strain) 
  x.m$strain <- gsub("Bifi_", "Bifidobacterium_", x.m$strain) 
  
  x.m$strain <- gsub("_", " ", x.m$strain)
  
  # Shared color scale
  fill_scale <- scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    limits = range(x.m$value, na.rm = TRUE)
  )
  
  # Heatmap generator
  make_heatmap <- function(data, title) {
    ggplot(data, aes(x = variable, y = strain, fill = value)) +
      geom_tile() +
      fill_scale +
      labs(
        title = title,
        x = "",
        y = "StrainGST",
        fill = "Prevalence (%)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.position = "right"
      ) +
      coord_fixed(ratio = 1)
  }
  
  if (split_n == 1) {
    p <- make_heatmap(x.m, "")
    return(p)
  }
  
  # Split strain names
  strains <- unique(x.m$strain)
  n <- length(strains)
  splits <- split(strains, cut(seq_along(strains), split_n, labels = FALSE))
  data_splits <- lapply(splits, function(s) x.m %>% filter(strain %in% s))
  
  # Generate plots and combine
  plots <- mapply(make_heatmap, data_splits, 
                  paste("Strain Frequency (", seq_along(data_splits), " of ", split_n, ")", sep=""), 
                  SIMPLIFY = FALSE)
  combined <- wrap_plots(plots, ncol = split_n) + plot_layout(guides = "collect") & theme(legend.position = "right")
  
  return(combined)
}


#function for loading straingr output
straingr_load = function(x, y=R21_sample_list){
  colnames(x) = c("sample1", "sample2", "ref", "scaffold", "length", "common", "commonPct", "single", "singlePct", "singleAgree", "singleAgreePct", "multi", "multiPct", "sharedAlleles", "sharedAllelesPct", "variants", "variantPct", "commonVariant", "commonVariantPct", "variantExact", "variantExactPct", "AnotB", "AnotBpct", "BnotA", "BnotApct", "Agaps", "AgapPct", "Bgaps", "BgapPct", "gapJaccardSim")
  #Append sample information to distance information
  x[,c("sample1_study", "sample1_area", "sample1_hhid", "sample1_tr", "sample1_gps_lat", "sample1_gps_lon", "sample1_location")] = R21_sample_list[match(x$sample1,R21_sample_list$sample_id),c("study", "area", "hhid", "tr", "gps_hhlatitude", "gps_hhlongitude", "location")]
  x[,c("sample2_study", "sample2_area", "sample2_hhid", "sample2_tr", "sample2_gps_lat", "sample2_gps_lon", "sample2_location")] = R21_sample_list[match(x$sample2,R21_sample_list$sample_id),c("study", "area", "hhid", "tr", "gps_hhlatitude", "gps_hhlongitude", "location")] 
  #Calculate HH distance
  for (i in 1:nrow(x)) {
    x$distance[i] = distm(c(x$sample1_gps_lon[i], x$sample1_gps_lat[i]), c(x$sample2_gps_lon[i], x$sample2_gps_lat[i]))
  }
  
  #Assigning HH_type between samples
  x = x %>% filter(sample1 != "stools_12093101_SM-NA1P5" & sample2 != "stools_12093101_SM-NA1P5")
  x[x$sample1_hhid == x$sample2_hhid, "HH_type"] = "Within"
  x[x$sample1_hhid != x$sample2_hhid, "HH_type"] = "Between"
  
  #Assign treatment info
  x[x$sample1_tr == x$sample2_tr & x$sample1_tr == "Water", "tr"] = "Water"
  x[x$sample1_tr == x$sample2_tr & x$sample1_tr == "Control", "tr"] = "Control"
  
  #Assign area info
  x[x$sample1_area == x$sample2_area & x$sample1_area == "rural", "area"] = "rural"
  x[x$sample1_area == x$sample2_area & x$sample1_area == "urban", "area"] = "urban"
  
  x[x$sample1_location == x$sample2_location, "location"] = "Same_loc"
  x[x$sample1_location != x$sample2_location, "location"] = "Different_loc"
  
  x$HH_type = factor(x$HH_type, levels = c("Within", "Between"))
  return(x)
}

# Derivitive of straingr_load
metadata_pair_load = function(x, y=R21_sample_list){
  #colnames(x) = c("sample1", "sample2", "ref", "scaffold", "length", "common", "commonPct", "single", "singlePct", "singleAgree", "singleAgreePct", "multi", "multiPct", "sharedAlleles", "sharedAllelesPct", "variants", "variantPct", "commonVariant", "commonVariantPct", "variantExact", "variantExactPct", "AnotB", "AnotBpct", "BnotA", "BnotApct", "Agaps", "AgapPct", "Bgaps", "BgapPct", "gapJaccardSim")
  #Append sample information to distance information
  #x[,c("sample1_study", "sample1_area", "sample1_hhid", "sample1_tr", "sample1_gps_lat", "sample1_gps_lon", "sample1_location")] = y[match(x$sample1,y$sample_id),c("study", "area", "hhid", "tr", "gps_hhlatitude", "gps_hhlongitude", "location")]
  #x[,c("sample2_study", "sample2_area", "sample2_hhid", "sample2_tr", "sample2_gps_lat", "sample2_gps_lon", "sample2_location")] = y[match(x$sample2,y$sample_id),c("study", "area", "hhid", "tr", "gps_hhlatitude", "gps_hhlongitude", "location")]
  # update on 071125
  x[,c("sample1_study", "sample1_area", "sample1_hhid", "sample1_tr", "sample1_gps_lat", "sample1_gps_lon", "sample1_location", "sample1_Age_group")] = y[match(x$sample1,y$sample_id),c("study", "area", "hhid", "tr", "gps_hhlatitude", "gps_hhlongitude", "location", "Age_group")]
  x[,c("sample2_study", "sample2_area", "sample2_hhid", "sample2_tr", "sample2_gps_lat", "sample2_gps_lon", "sample2_location", "sample2_Age_group")] = y[match(x$sample2,y$sample_id),c("study", "area", "hhid", "tr", "gps_hhlatitude", "gps_hhlongitude", "location", "Age_group")]
  #Calculate HH distance
  for (i in 1:nrow(x)) {
    x$distance[i] = distm(c(x$sample1_gps_lon[i], x$sample1_gps_lat[i]), c(x$sample2_gps_lon[i], x$sample2_gps_lat[i]))
  }
  
  #Assigning HH_type between samples
  x = x %>% filter(sample1 != "stools_12093101_SM-NA1P5" & sample2 != "stools_12093101_SM-NA1P5")
  x[x$sample1_hhid == x$sample2_hhid, "HH_type"] = "Within"
  x[x$sample1_hhid != x$sample2_hhid, "HH_type"] = "Between"
  
  #Assign treatment info
  x[x$sample1_tr == x$sample2_tr & x$sample1_tr == "Water", "tr"] = "Water"
  x[x$sample1_tr == x$sample2_tr & x$sample1_tr == "Control", "tr"] = "Control"
  
  #Assign area info
  x[x$sample1_area == x$sample2_area & x$sample1_area == "rural", "area"] = "rural"
  x[x$sample1_area == x$sample2_area & x$sample1_area == "urban", "area"] = "urban"
  
  x[x$sample1_location == x$sample2_location, "location"] = "Same_loc"
  x[x$sample1_location != x$sample2_location, "location"] = "Different_loc"
  
  x$HH_type = factor(x$HH_type, levels = c("Within", "Between"))
  return(x)
}

#function for permutation test
permutation.test = function(a,b){
  #a = strain.share[strain.share$subcategory == overall.perm[i,1],] #Extract first category
  #b = strain.share[strain.share$subcategory == overall.perm[i,2],] #Extract second category
  #
  set.seed(1992)
  rep1 = replicate(1000, sample(c(a,b), replace = F))
  a1 = rep1[1:length(a),]
  b1 = rep1[(length(a)+1):(length(a) + length(b)),]
  
  diffs = apply(a1, 2, mean) - apply(b1, 2, mean)
  if(mean(a) == 0 & mean(b) == 0){
    p = NA
  }
  else{
    p = sum(abs(diffs) >= (abs(mean(a) - mean(b))))/1000
  }
  return(p)
}

#Calculate HH level strain-sharing rates
get_hh_rate = function(org, y= R21_sample_list, z=staingr.org.ML){
  hh.comb.ML = combn(y$hhid %>% unique, 2) %>% t()
  hh.comb.ML = rbind(hh.comb.ML, data.frame(V1=y$hhid %>% unique, V2=y$hhid %>% unique))
  colnames(hh.comb.ML) = c("HH1", "HH2")
  
  for (i in 1:nrow(hh.comb.ML)) {
    hh.list = c(hh.comb.ML[i,"HH1"], hh.comb.ML[i,"HH2"])
    if(hh.comb.ML[i,"HH1"] == hh.comb.ML[i,"HH2"]){ #Within HH
      num.hh = y[y$hhid == hh.comb.ML[i,"HH1"], "sample_id"] %>% unique %>% length
      hh.comb.ML$all.pair[i] = num.hh*(num.hh-1)/2 #number of all possible pair
      sharing.hh = z %>% filter(sample1_hhid %in% hh.list & sample2_hhid %in% hh.list)
      if(org != "All"){
        sharing.hh = sharing.hh %>% filter(genus == org)
      }
      sharing.hh = sharing.hh %>% select(., c("sample1", "sample2")) %>% unique()
      sharing.hh <- t(apply(sharing.hh, 1, sort)) %>% unique() # sort samples in pairs
      hh.comb.ML$sharing.pair[i] = sharing.hh %>% nrow() # Numb pair with the identical strains
    }
    if(hh.comb.ML[i,"HH1"] != hh.comb.ML[i,"HH2"]){ #Between households
      num.hh1 = y[y$hhid == hh.comb.ML[i,"HH1"], "sample_id"] %>% unique %>% length
      num.hh2 = y[y$hhid == hh.comb.ML[i,"HH2"], "sample_id"] %>% unique %>% length
      hh.comb.ML$all.pair[i] = num.hh1*num.hh2 #number of all possible pair
      sharing.hh = z %>% filter(sample1_hhid != sample2_hhid & sample1_hhid %in% hh.list & sample2_hhid %in% hh.list)
      if(org != "All"){
        sharing.hh = sharing.hh %>% filter(genus == org)
      }
      sharing.hh = sharing.hh %>% select(., c("sample1", "sample2")) %>% unique()
      sharing.hh <- t(apply(sharing.hh, 1, sort)) %>% unique() # sort samples in pairs
      hh.comb.ML$sharing.pair[i] = sharing.hh %>% nrow() # Numb pair with the identical strains
    }
  }
  hh.comb.ML$genus = org
  #Remove pair of households with no possible pairs (in case only one child in HH)
  hh.comb.ML = hh.comb.ML %>% filter(all.pair != 0)
  #calculate rates
  hh.comb.ML$rate = hh.comb.ML$sharing.pair/hh.comb.ML$all.pair
  #Append HH information to distance information
  hh.comb.ML[,c("HH1_study", "HH1_area", "HH1_tr", "HH1_gps_lat", "HH1_gps_lon", "HH1_location")] = y[match(hh.comb.ML$HH1, y$hhid),c("study", "area", "tr", "gps_hhlatitude", "gps_hhlongitude", "location")]
  hh.comb.ML[,c("HH2_study", "HH2_area", "HH2_tr", "HH2_gps_lat", "HH2_gps_lon", "HH2_location")] = y[match(hh.comb.ML$HH2, y$hhid),c("study", "area", "tr", "gps_hhlatitude", "gps_hhlongitude", "location")]
  #Add combination information
  hh.comb.ML[hh.comb.ML$HH1 == hh.comb.ML$HH2, "HH_type"] = "Within"
  hh.comb.ML[hh.comb.ML$HH1 != hh.comb.ML$HH2, "HH_type"] = "Between"
  hh.comb.ML[hh.comb.ML$HH1_area == hh.comb.ML$HH2_area & hh.comb.ML$HH1_area == "urban", "area"] = "urban"
  hh.comb.ML[hh.comb.ML$HH1_area == hh.comb.ML$HH2_area & hh.comb.ML$HH1_area == "rural", "area"] = "rural"
  hh.comb.ML[hh.comb.ML$HH1_area != hh.comb.ML$HH2_area, "area"] = "Between studies"
  hh.comb.ML[hh.comb.ML$HH1_tr == hh.comb.ML$HH2_tr & hh.comb.ML$HH1_tr == "Water", "tr"] = "Water"
  hh.comb.ML[hh.comb.ML$HH1_tr == hh.comb.ML$HH2_tr & hh.comb.ML$HH1_tr == "Control", "tr"] = "Control"
  hh.comb.ML[hh.comb.ML$HH1_tr != hh.comb.ML$HH2_tr, "tr"] = "Between treatment"
  return(hh.comb.ML)
}

#Comparison of HH level strain-sharing rates within vs. between
get.wit.bet.test = function(x){
  wit.bet.test = data.frame(matrix(ncol=5)) 
  colnames(wit.bet.test) = c("area", "tr", "within", "between", "p.value")
  n=1
  for (i in c("urban", "rural")) {
    for (j in c("Overall", "Water", "Control")) {
      if(j == "Overall"){
        wit.bet.test[n,1] = i
        wit.bet.test[n,2] = j
        a = x %>% filter(area == i & HH_type == "Within") %>% select(.,rate) %>% unlist()
        b = x %>% filter(area == i & HH_type == "Between") %>% select(.,rate) %>% unlist()
        wit.bet.test[n,3] = mean(a)
        wit.bet.test[n,4] = mean(b)
        wit.bet.test[n,5] = permutation.test(a,b)
      }
      else{
        wit.bet.test[n,1] = i
        wit.bet.test[n,2] = j
        a = x %>% filter(area == i & HH_type == "Within" & tr == j) %>% select(.,rate)  %>% unlist()
        b = x %>% filter(area == i & HH_type == "Between" & tr == j) %>% select(.,rate)  %>% unlist()
        wit.bet.test[n,3] = mean(a)
        wit.bet.test[n,4] = mean(b)
        wit.bet.test[n,5] = permutation.test(a,b)
      }
      n=n+1
    }
  }
  return(wit.bet.test)
}

# Get bootstrapped values
# Uses its own RNG stream (seeded, then restored) so results don't depend on
# call order relative to permutation.test(), which resets the global seed.
boot_ci <- function(x, n_boot = 1000, conf.level = 0.95, seed = 617) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
  set.seed(seed)
  boot_means <- replicate(n_boot, mean(sample(x, size = length(x), replace = TRUE)))
  if (!is.null(old_seed)) .Random.seed <<- old_seed else rm(".Random.seed", envir = .GlobalEnv)
  alpha <- 1 - conf.level
  quantile(boot_means, probs = c(alpha/2, 1 - alpha/2))
}


#Comparison of HH level strain-sharing rates within vs. between
get.wit.bet.test.CI = function(x, n_boot = 1000){
  set.seed(617) # For reproducibility  
  wit.bet.test = data.frame(matrix(ncol=9))
  colnames(wit.bet.test) = c("area", "tr", "within", "within_lower", "within_upper", "between", "between_lower", "between_upper", "p.value")
  n=1
  for (i in c("urban", "rural")) {
    for (j in c("Overall", "Water", "Control")) {
      if(j == "Overall"){
        a = x %>% filter(area == i & HH_type == "Within") %>% select(.,rate) %>% unlist()
        b = x %>% filter(area == i & HH_type == "Between") %>% select(.,rate) %>% unlist()
      } else {
        a = x %>% filter(area == i & HH_type == "Within" & tr == j) %>% select(.,rate)  %>% unlist()
        b = x %>% filter(area == i & HH_type == "Between" & tr == j) %>% select(.,rate)  %>% unlist()
      }
      a_ci = boot_ci(a, n_boot)
      b_ci = boot_ci(b, n_boot)
      wit.bet.test[n,] = c(i, j, mean(a), a_ci[1], a_ci[2], mean(b), b_ci[1], b_ci[2], permutation.test(a,b))
      n=n+1
    }
  }
  return(wit.bet.test)
}


#Compare rates (Rural vs. Urban) with CI
get.rural.urban.test.CI = function(x, n_boot = 1000){
  rural.urban.test = data.frame(matrix(ncol=9)) 
  colnames(rural.urban.test) = c("HH_type", "tr", "rural", "rural_lower", "rural_upper", "urban", "urban_lower", "urban_upper", "p.value")
  n=1
  for (i in c("Within", "Between")) {
    for (j in c("Overall", "Water", "Control")) {
      if(j == "Overall"){
        a = x %>% filter(HH_type == i & area == "rural") %>% select(.,rate) %>% unlist()
        b = x %>% filter(HH_type == i & area == "urban") %>% select(.,rate) %>% unlist()
      }
      else{
        a = x %>% filter(HH_type == i & area == "rural" & tr == j) %>% select(.,rate)  %>% unlist()
        b = x %>% filter(HH_type == i & area == "urban" & tr == j) %>% select(.,rate)  %>% unlist()
      }
      a_ci = boot_ci(a, n_boot)
      b_ci = boot_ci(b, n_boot)
      rural.urban.test[n,] = c(i, j, mean(a), a_ci[1], a_ci[2], mean(b), b_ci[1], b_ci[2], permutation.test(a,b))
      n=n+1
    }
  }
  return(rural.urban.test)
}

#Comparison of HH level strain-sharing rates within and between studies (only for bw HH strain-sharing rates) with CI
get.bw.studies.test.CI = function(x, n_boot = 1000){
  set.seed(617) # For reproducibility

  # Compute mean + bootstrap CI for each group exactly once, so a group's
  # stats don't change depending on which pairwise row it appears in
  groups = c("urban", "rural", "Between studies")
  group_stats = lapply(groups, function(g){
    v = x %>% filter(area == g & HH_type == "Between") %>% select(.,rate) %>% unlist()
    ci = boot_ci(v, n_boot)
    data.frame(group = g, mean = mean(v), lower = ci[1], upper = ci[2])
  }) %>% bind_rows()

  bw.studies.test = data.frame(matrix(ncol=10))
  colnames(bw.studies.test) = c("study1", "study2", "study1_mean", "study1_upper", "study1_lower","study2_mean", "study2_upper", "study2_lower", "p.value", "adj.p.value")

  bw.studies.test[1,c("study1", "study2")] = c("urban", "rural")
  bw.studies.test[2,c("study1", "study2")] = c("urban", "Between studies")
  bw.studies.test[3,c("study1", "study2")] = c("rural", "Between studies")

  for (i in 1:nrow(bw.studies.test)) {
    s1 = group_stats %>% filter(group == bw.studies.test[i,"study1"])
    s2 = group_stats %>% filter(group == bw.studies.test[i,"study2"])

    bw.studies.test[i,"study1_mean"] = s1$mean
    bw.studies.test[i,"study1_upper"] = s1$upper
    bw.studies.test[i,"study1_lower"] = s1$lower
    bw.studies.test[i,"study2_mean"] = s2$mean
    bw.studies.test[i,"study2_upper"] = s2$upper
    bw.studies.test[i,"study2_lower"] = s2$lower

    a = x %>% filter(area ==  bw.studies.test[i,"study1"] & HH_type == "Between") %>% select(.,rate) %>% unlist()
    b = x %>% filter(area == bw.studies.test[i,"study2"] & HH_type == "Between") %>% select(.,rate) %>% unlist()
    bw.studies.test[i,"p.value"] = permutation.test(a,b)
  }
  bw.studies.test$adj.p.value = p.adjust(bw.studies.test$p.value, method = "BH")
  return(bw.studies.test)
}

#Comparison of HH level strain-sharing rates Rural vs. Urban
get.rural.urban.test = function(x){
  rural.urban.test = data.frame(matrix(ncol=5)) 
  colnames(rural.urban.test) = c("HH_type", "tr", "rural", "urban", "p.value")
  n=1
  for (i in c("Within", "Between")) {
    for (j in c("Overall", "Water", "Control")) {
      if(j == "Overall"){
        rural.urban.test[n,1] = i
        rural.urban.test[n,2] = j
        a = x %>% filter(HH_type == i & area == "rural") %>% select(.,rate) %>% unlist()
        b = x %>% filter(HH_type == i & area == "urban") %>% select(.,rate) %>% unlist()
        rural.urban.test[n,3] = mean(a)
        rural.urban.test[n,4] = mean(b)
        rural.urban.test[n,5] = permutation.test(a,b)
      }
      else{
        rural.urban.test[n,1] = i
        rural.urban.test[n,2] = j
        a = x %>% filter(HH_type == i & area == "rural" & tr == j) %>% select(.,rate)  %>% unlist()
        b = x %>% filter(HH_type == i & area == "urban" & tr == j) %>% select(.,rate)  %>% unlist()
        rural.urban.test[n,3] = mean(a)
        rural.urban.test[n,4] = mean(b)
        rural.urban.test[n,5] = permutation.test(a,b)
      }
      n=n+1
    }
  }
  return(rural.urban.test)
}

#Comparison of HH level strain-sharing rates Water vs. Control
get.wat.cont.test = function(x){
  wat.cont.test = data.frame(matrix(ncol=5))
  colnames(wat.cont.test) = c("area", "HH_type", "water", "control", "p.value")
  n=1
  for (i in c("urban", "rural")) {
    for (j in c("Within", "Between")) {
      wat.cont.test[n,1] = i
      wat.cont.test[n,2] = j
      a = x %>% filter(area == i & HH_type == j & tr == "Water") %>% select(.,rate)  %>% unlist()
      b = x %>% filter(area == i & HH_type == j & tr == "Control") %>% select(.,rate)  %>% unlist()
      wat.cont.test[n,3] = mean(a)
      wat.cont.test[n,4] = mean(b)
      wat.cont.test[n,5] = permutation.test(a,b)
      n=n+1
    }
  }
  return(wat.cont.test)
}

#Comparison of HH level strain-sharing rates within and between studies (only for bw HH strain-sharing rates)
get.bw.studies.test = function(x){
  bw.studies.test = data.frame(matrix(ncol=6))
  colnames(bw.studies.test) = c("study1", "study2", "study1_mean", "study2_mean", "p.value", "adj.p.value")
  
  a = x %>% filter(area == "urban" & HH_type == "Between") %>% select(.,rate) %>% unlist() # urban
  b = x %>% filter(area == "rural" & HH_type == "Between") %>% select(.,rate) %>% unlist() # rural
  c = x %>% filter(area == "Between studies" & HH_type == "Between") %>% select(.,rate) %>% unlist() # bw studies
  
  bw.studies.test[1,1:2] = c("urban", "rural")
  bw.studies.test[2,1:2] = c("urban", "Between studies")
  bw.studies.test[3,1:2] = c("rural", "Between studies")
  
  for (i in 1:nrow(bw.studies.test)) {
    a = x %>% filter(area ==  bw.studies.test[i,1] & HH_type == "Between") %>% select(.,rate) %>% unlist() # urban
    b = x %>% filter(area == bw.studies.test[i,2] & HH_type == "Between") %>% select(.,rate) %>% unlist() # rural
    
    bw.studies.test[i,3] = mean(a)
    bw.studies.test[i,4] = mean(b)
    bw.studies.test[i,5] = permutation.test(a,b)
  }
  bw.studies.test$adj.p.value = p.adjust(bw.studies.test$p.value, method = "BH")
  return(bw.studies.test)
}




# For importing MASH output in phylip format
parseDistanceDF = function(phylip_file) {
  # Read the first line of the phylip file to find out how many sequences/samples it contains
  # https://forum.mothur.org/t/importing-dist-matrix-into-r/1266/6
  temp_connection = file(phylip_file, 'r')
  len = readLines(temp_connection, n=1)
  len = as.numeric(len)
  len = len +1
  close(temp_connection)
  phylip_data = read.table(phylip_file, fill=T, row.names=1, skip=1, col.names=1:len)
  colnames(phylip_data) <- row.names(phylip_data)
  return(phylip_data)
}

# function for Fisher's test on strain-sharing percentage betwee age groups by treatment status
calculate_percent_fisher = function(x){
  result = data.frame(matrix(ncol=5))
  colnames(result) = c("cat", "genus", "water_percent", "control_percent", "p_value")
  n=1
  for (i in unique(x$cat)) {
    for (j in unique(x$genus)) {
      water = x %>% filter(cat == i & genus == j & tr == "Water") %>% select("Yes", "No")
      control = x %>% filter(cat == i & genus == j & tr == "Control") %>% select("Yes", "No")
      # Extract Yes/No counts dynamically
      water_counts <- as.numeric(water[1, ])    # Extract first row as numeric
      control_counts <- as.numeric(control[1, ]) # Extract first row as numeric
      # Create a contingency table
      contingency_table <- matrix(c(water_counts, control_counts),
                                  nrow = 2, byrow = TRUE,
                                  dimnames = list(Treatment = c("Water", "Control"),
                                                  Outcome = c("Yes", "No")))
      fisher_result <- fisher.test(contingency_table)
      # add results
      result[n,"cat"] = i
      result[n,"genus"] = j
      result[n,"water_percent"] = water_counts[1]/(water_counts[1] + water_counts[2])*100
      result[n,"control_percent"] = control_counts[1]/(control_counts[1] + control_counts[2])*100
      result[n,"p_value"] = fisher_result$p.value
      n=n+1
    }
  }
  return(result)
}



process_nuccore_file <- function(file) {
  # Read all lines from the file
  lines <- readLines(file)
  
  # Split lines into records using blank lines as separators
  records <- list()
  current <- character()
  for (line in lines) {
    if (nchar(trimws(line)) == 0) {
      if (length(current) > 0) {
        records[[length(records) + 1]] <- current
        current <- character()
      }
    } else {
      current <- c(current, line)
    }
  }
  if (length(current) > 0) {
    records[[length(records) + 1]] <- current
  }
  
  # Prepare vectors for scaffold IDs and organism names
  scaffold <- character()
  organism <- character()
  
  for (rec in records) {
    ## Extract scaffold ID from the accession line (e.g. NZ_CP020055.1)
    acc_line <- rec[grepl("^(NZ|NC)_", rec)]
    if (length(acc_line) > 0) {
      acc <- strsplit(acc_line[1], "\\s+")[[1]][1]
    } else {
      acc <- NA
    }
    
    ## Process the organism name from the first line
    first_line <- rec[1]
    # Remove the leading record number and dot, then take text before the comma
    title <- sub("^[0-9]+\\.\\s*", "", first_line)
    title <- strsplit(title, ",")[[1]][1]
    
    # Split the title into tokens
    tokens <- strsplit(title, "\\s+")[[1]]
    
    ## Determine the strain token:
    ## If the third token is "strain", use the fourth token as the strain identifier.
    if (length(tokens) >= 4 && tokens[3] == "strain") {
      strain_token <- tokens[4]
    } else if (length(tokens) >= 3) {
      strain_token <- tokens[3]
    } else {
      strain_token <- NA
    }
    
    # Abbreviate the genus (first token) to its first 4 letters and get species from the second token
    genus_abbrev <- ifelse(length(tokens) >= 1, substr(tokens[1], 1, 4), NA)
    species <- ifelse(length(tokens) >= 2, tokens[2], NA)
    
    ## Special case: If the organism is Campylobacter upsaliensis RM3195, append "_map"
    if (tokens[1] == "Campylobacter" && tokens[2] == "upsaliensis" && strain_token == "RM3195") {
      org_name <- "Camp_upsaliensis_RM3195_map"
    } else {
      org_name <- paste0(genus_abbrev, "_", species, "_", strain_token)
    }
    
    scaffold <- c(scaffold, acc)
    organism <- c(organism, org_name)
  }
  
  # Return a data frame with scaffold IDs and processed organism names
  df <- data.frame(scaffold = scaffold, organism = organism, stringsAsFactors = FALSE)
  return(df)
}




# Getting proportion of strain-sharing sample pairs within a range of HH distances
get_sharing_dist = function(pairs, min, max){
  dist.bin = data.frame()
  for (g in c("All", "Bact", "Bifi", "Esch", "Camp_MAG", "Entc", "Entb", "Kleb","commensal", "non_commensal")) {
    for (d in seq(min, max, by = 10)) {
      i = paste0(g, ".share")
      if(pairs %>% filter(HH_type == "Between" & distance <= d) %>% pull(all_of(i)) %>% unique() %>% length() == 1){
        new.bin = pairs %>% filter(HH_type == "Between" & distance <= d) %>% select(all_of(i)) %>% table %>% as.data.frame 
        new.bin[,i] = as.character(new.bin[,i]) #Change class to charactor to add a new row
        new.bin[2,] = c("Yes", 0) # Manually add "Yes"
        new.bin[,i] = as.factor(new.bin[,i]) #Change back to factor level
        new.bin[,"Freq"] = as.numeric(new.bin[,"Freq"])
        new.bin = new.bin %>% pivot_wider(names_from = all_of(i), values_from = Freq) %>% mutate(percent = (Yes / (Yes + No)) * 100, genus = g, distance = d)
      }
      else if(pairs %>% filter(HH_type == "Between" & distance <= d) %>% nrow == 0){
        next
      }
      else{
        new.bin = pairs %>% filter(HH_type == "Between" & distance <= d) %>% select(all_of(i)) %>% table %>% as.data.frame %>% pivot_wider(names_from = all_of(i), values_from = Freq) %>% mutate(percent = (Yes / (Yes + No)) * 100, genus = g, distance = d)
      }
      dist.bin = rbind(dist.bin, new.bin)
    }
  }
  return(dist.bin)
}

