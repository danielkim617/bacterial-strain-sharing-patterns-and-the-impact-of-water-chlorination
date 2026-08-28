# this script is to get sample pairs with the same strain identified from straingst results. The output will be used for strain-sharing analysis.
source(file.path(here::here(),"0-config.R"))
#Load R object
R21_sample_list = readRDS(file.path(here::here(), "data/Kraken.phy.genus.alpha.rds"))

#Loading straingst results
ecoli.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.ecoli.tsv"), sep = "\t")
camp.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Campylobacter.tsv"), sep = "\t")
camp_all.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Campylobacter_all.tsv"), sep = "\t")
camp_MAG.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Campylobacter_MAG.tsv"), sep = "\t")
enterob.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Enterobacter.tsv"), sep = "\t")
enterco.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Enterococcus.tsv") , sep = "\t")
kleb.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Klebsiella.tsv") , sep = "\t")
staph.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Staphylococcus.tsv") , sep = "\t")
#From Colin's DB
bacteroides.strains = read.table(file.path(here::here(), "data/colin_DB_straingst/bacteroides_straingst.tsv"), sep = "\t")
bifidobacterium.strains = read.table(file.path(here::here(),"data/colin_DB_straingst/bifidobacterium_straingst.tsv"), sep = "\t")

#Get pair of samples with the same strain identified
ecoli.strains.pair_sorted =  get.pair(ecoli.strains)
camp.strains.pair_sorted =  get.pair(camp.strains)
camp_all.strains.pair_sorted =  get.pair(camp_all.strains)
camp_MAG.strains.pair_sorted =  get.pair(camp_MAG.strains)
enterob.strains.pair_sorted =  get.pair(enterob.strains)
enterco.strains.pair_sorted =  get.pair(enterco.strains)
kleb.strains.pair_sorted =  get.pair(kleb.strains)
#staph.strains.pair_sorted =  get.pair(staph.strains) # No pair exist as only detected in one sample
bacteroides.strains.pair_sorted =  get.pair(bacteroides.strains) #From this point, sample pairs between rural and urban area will not be included.
bifidobacterium.strains.pair_sorted =  get.pair(bifidobacterium.strains)

#Save into files
write.table(ecoli.strains.pair_sorted, file.path(here::here(),"data/comb_ecoli.strains.pair_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(camp.strains.pair_sorted, file.path(here::here(),"data/camp.strains.pair_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(camp_all.strains.pair_sorted, file.path(here::here(),"data/camp_all.strains.pair_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(camp_MAG.strains.pair_sorted, file.path(here::here(),"data/camp_MAG.strains.pair_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(enterob.strains.pair_sorted, file.path(here::here(),"data/enterob.strains.pair_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(enterco.strains.pair_sorted, file.path(here::here(),"data/enterco.strains.pair_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(kleb.strains.pair_sorted, file.path(here::here(),"data/kleb.strains.pair_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(bacteroides.strains.pair_sorted, file.path(here::here(),"data/bacteroides.strains.pair_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(bifidobacterium.strains.pair_sorted, file.path(here::here(),"data/bifidobacterium.strains.pair_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)

## To get sample pairs between villages for commensals to determine baseline strain-sharing rates
# function for loading StrainGR output
get.pair_bw = function(x){
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
  strains.pair_sorted <- t(apply(strains.pair[,c("sample1", "sample2")], 1, sort))
  # Convert back to a data frame
  strains.pair_sorted <- as.data.frame(strains.pair_sorted)
  # Remove duplicate rows
  strains.pair_sorted <- unique(strains.pair_sorted)
  ### retain pairs only between urban and rural areas
  strains.pair_sorted.1 = strains.pair_sorted %>% filter(sample1 %in% rural.samples & sample2 %in% urban.samples)
  strains.pair_sorted.2 = strains.pair_sorted %>% filter(sample1 %in% urban.samples & sample2 %in% rural.samples)
  strains.pair_sorted = rbind(strains.pair_sorted.1, strains.pair_sorted.2)
    ###
  return(strains.pair_sorted)
}

# sample pairs between areas
# Commensals
bacteroides.strains.pair_bw_sorted =  get.pair_bw(bacteroides.strains)
bifidobacterium.strains.pair_bw_sorted =  get.pair_bw(bifidobacterium.strains)
# Poten. pathogenic
ecoli.strains.pair_bw_sorted =  get.pair_bw(ecoli.strains)
camp.strains.pair_bw_sorted =  get.pair_bw(camp.strains)
camp_all.strains.pair_bw_sorted =  get.pair_bw(camp_all.strains)
camp_MAG.strains.pair_bw_sorted =  get.pair_bw(camp_MAG.strains)
enterob.strains.pair_bw_sorted =  get.pair_bw(enterob.strains)
enterco.strains.pair_bw_sorted =  get.pair_bw(enterco.strains)
kleb.strains.pair_bw_sorted =  get.pair_bw(kleb.strains)

#Save into files
write.table(bacteroides.strains.pair_bw_sorted, file.path(here::here(),"data/bacteroides.strains.pair_bw_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(bifidobacterium.strains.pair_bw_sorted, file.path(here::here(),"data/bifidobacterium.strains.pair_bw_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)

write.table(ecoli.strains.pair_bw_sorted, file.path(here::here(),"data/ecoli.strains.pair_bw_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(camp.strains.pair_bw_sorted, file.path(here::here(),"data/camp.strains.pair_bw_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(camp_all.strains.pair_bw_sorted, file.path(here::here(),"data/camp_all.strains.pair_bw_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(camp_MAG.strains.pair_bw_sorted, file.path(here::here(),"data/camp_MAG.strains.pair_bw_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(enterob.strains.pair_bw_sorted, file.path(here::here(),"data/enterob.strains.pair_bw_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(enterco.strains.pair_bw_sorted, file.path(here::here(),"data/enterco.strains.pair_bw_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(kleb.strains.pair_bw_sorted, file.path(here::here(),"data/kleb.strains.pair_bw_sorted.tsv"), sep = "\t", row.names = F, quote = F, col.names = F)
