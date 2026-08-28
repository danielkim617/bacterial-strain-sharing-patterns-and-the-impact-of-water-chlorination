# this script is to load StrainGR outputs and add species information for the random forest model. The output will be used for random forest-based strain sharing calls.
source(file.path(here::here(),"0-config.R"))
#Load R object
R21_sample_list = readRDS(file.path(here::here(), "data/Kraken.phy.genus.alpha.rds"))

#Load StrainGR outputs and add species information for the random forest model
#E. coli
staingr.ecoli = read.table(file.path(here::here(), "data/straingr_comb/ecoli_new.straingr.compare.combined.tsv") , header = F, sep = "\t")
staingr.ecoli = straingr_load(staingr.ecoli)
ecoli.ref = read.table(file.path(here::here(), "data/StrainGE_ML/e.coli.ref.tsv"), sep = "\t", header = F)
colnames(ecoli.ref) = c("scaffold", "species")
staingr.ecoli$species =  ecoli.ref[match(staingr.ecoli$scaffold, ecoli.ref$scaffold), "species"] 

#Bacteroides
bacteroides.straingr = read.table(file.path(here::here(), "data/StrainGE_ML/bacteroides_combined_strangr.tsv") , sep = "\t", header = F)
bacteroides.straingr = straingr_load(bacteroides.straingr)
bacteroides.ref = read.table(file.path(here::here(), "data/StrainGE_ML/bacteroides.ref.tsv"), sep = "\t", header = F)
colnames(bacteroides.ref) = c("scaffold", "species")
bacteroides.ref[bacteroides.ref$species == "Bact_sp", "species"] = "Bact_sp."
bacteroides.straingr$species =  bacteroides.ref[match(bacteroides.straingr$scaffold, bacteroides.ref$scaffold), "species"] 

#Bifidobacterium
bifidobacterium.straingr = read.table(file.path(here::here(), "data/StrainGE_ML/bifidobacterium_combined_strangr.tsv"), sep = "\t", header = F)
bifidobacterium.straingr = straingr_load(bifidobacterium.straingr)
bifidobacterium.ref = read.table(file.path(here::here(), "data/StrainGE_ML/bifidobacterium.ref.tsv"), sep = "\t", header = F)
colnames(bifidobacterium.ref) = c("scaffold", "species")
bifidobacterium.ref[bifidobacterium.ref$species == "Bifi_sp", "species"] = "Bifi_sp."
bifidobacterium.straingr$species =  bifidobacterium.ref[match(bifidobacterium.straingr$scaffold, bifidobacterium.ref$scaffold), "species"] 

#Campylobacter
camp.all.straingr = read.table(file.path(here::here(), "data/StrainGE_ML/Campylobacter_all.straingr.compare.combined.tsv"), sep = "\t", header = F)
camp.all.straingr = straingr_load(camp.all.straingr)
camp.all.ref = read.table(file.path(here::here(), "data/StrainGE_ML/camp.all.ref.tsv"), sep = "\t", header = F)
colnames(camp.all.ref) = c("scaffold", "species")
camp.all.ref[camp.all.ref$species == "Camp_sp", "species"] = "Camp_sp."
camp.all.straingr$species =  camp.all.ref[match(camp.all.straingr$scaffold, camp.all.ref$scaffold), "species"] 

#Campylobacter (concateated MAG)
camp.mag.straingr = read.table(file.path(here::here(), "data/StrainGE_ML/Campylobacter_MAG.straingr.compare.combined.tsv"), sep = "\t", header = F)
camp.mag.straingr = straingr_load(camp.mag.straingr)
#As there are only few number of reference genomes detected...(I will manually add those)
camp.mag.ref = data.frame(scaffold = camp.mag.straingr$scaffold %>% unique)
camp.mag.ref$species = c("Camp_infans", "Camp_sp.", "Camp_infans", "Camp_coli", "Camp_jejuni", "Camp_jejuni", "Camp_jejuni")
camp.mag.straingr$species =  camp.mag.ref[match(camp.mag.straingr$scaffold, camp.mag.ref$scaffold), "species"] 

#Enterococcus
entc.all.straingr = read.table(file.path(here::here(), "data/StrainGE_ML/Enterococcus.straingr.compare.combined.tsv"), sep = "\t", header = F)
entc.all.straingr = straingr_load(entc.all.straingr)
entc.all.ref = read.table(file.path(here::here(), "data/StrainGE_ML/enterococcus.ref.tsv"), sep = "\t", header = F)
colnames(entc.all.ref) = c("scaffold", "species")
#entc.all.ref[entc.all.ref$species == "Camp_sp", "species"] = "Camp_sp."
entc.all.straingr$species =  entc.all.ref[match(entc.all.straingr$scaffold, entc.all.ref$scaffold), "species"] 

#Klebsiella
kleb.straingr = read.table(file.path(here::here(), "data/StrainGE_ML/R21_strainGR_combined_kleb_summary.tsv"), sep = "\t", header = F)
kleb.straingr = kleb.straingr %>% filter(V1 != "sample1") #remove rows in the middle containing header information
kleb.straingr = straingr_load(kleb.straingr)
klebsiella.ref = read.table(file.path(here::here(), "data/StrainGE_ML/klebsiella.ref.tsv"), sep = "\t", header = F)
colnames(klebsiella.ref) = c("scaffold", "species")
klebsiella.ref[klebsiella.ref$species == "Kleb_sp", "species"] = "Kleb_sp."
kleb.straingr$species =  klebsiella.ref[match(kleb.straingr$scaffold, klebsiella.ref$scaffold), "species"] 

#Enterobacter
entb.straingr = read.table(file.path(here::here(), "data/StrainGE_ML/Enterobacter.straingr.compare.combined.tsv"), sep = "\t", header = F)
entb.straingr = straingr_load(entb.straingr)
entb.ref = read.table(file.path(here::here(), "data/StrainGE_ML/enterobacter.ref.tsv"), sep = "\t", header = F)
colnames(entb.ref) = c("scaffold", "species")
entb.ref[entb.ref$species == "Entb_sp", "species"] = "Entb_sp."
entb.straingr$species =  entb.ref[match(entb.straingr$scaffold, entb.ref$scaffold), "species"] 

#Write into files
write.table(staingr.ecoli, file.path(here::here(), "data/StrainGE_ML/e.coli.ML.input.tsv"), sep = "\t", row.names = F, quote = F)
write.table(bifidobacterium.straingr, file.path(here::here(), "data/StrainGE_ML/bifidobacterium.ML.input.tsv"), sep = "\t", row.names = F, quote = F)
write.table(bacteroides.straingr, file.path(here::here(), "data/StrainGE_ML/bacteroides.ML.input.tsv"), sep = "\t", row.names = F, quote = F)
write.table(camp.all.straingr, file.path(here::here(), "data/StrainGE_ML/camp.all.ML.input.tsv"), sep = "\t", row.names = F, quote = F)
write.table(camp.mag.straingr, file.path(here::here(), "data/StrainGE_ML/camp.mag.ML.input.tsv"), sep = "\t", row.names = F, quote = F)
write.table(entc.all.straingr, file.path(here::here(), "data/StrainGE_ML/entc.ML.input.tsv") , sep = "\t", row.names = F, quote = F)
write.table(entb.straingr, file.path(here::here(), "data/StrainGE_ML/entb.ML.input.tsv") , sep = "\t", row.names = F, quote = F)
write.table(kleb.straingr, file.path(here::here(), "data/StrainGE_ML/kleb.ML.input.tsv") , sep = "\t", row.names = F, quote = F)

################################################################################################################################################
################################################################################################################################################

# For strain comparisons between studies
#E. coli
bw_staingr.ecoli = read.table(file.path(here::here(), "data/straingr_comb/bw_ecoli_new.straingr.compare.combined.tsv") , header = F, sep = "\t")
bw_staingr.ecoli = straingr_load(bw_staingr.ecoli)
ecoli.ref = read.table(file.path(here::here(), "data/StrainGE_ML/e.coli.ref.tsv"), sep = "\t", header = F)
colnames(ecoli.ref) = c("scaffold", "species")
bw_staingr.ecoli$species =  ecoli.ref[match(bw_staingr.ecoli$scaffold, ecoli.ref$scaffold), "species"] 

#Bacteroides
bw_bacteroides.straingr = read.table(file.path(here::here(), "data/straingr_comb/bw_bact.straingr.compare.combined.tsv") , sep = "\t", header = F)
bw_bacteroides.straingr = straingr_load(bw_bacteroides.straingr)
bacteroides.ref = read.table(file.path(here::here(), "data/StrainGE_ML/bacteroides.ref.tsv"), sep = "\t", header = F)
colnames(bacteroides.ref) = c("scaffold", "species")
bacteroides.ref[bacteroides.ref$species == "Bact_sp", "species"] = "Bact_sp."
bw_bacteroides.straingr$species =  bacteroides.ref[match(bw_bacteroides.straingr$scaffold, bacteroides.ref$scaffold), "species"] 

#Bifidobacterium
bw_bifidobacterium.straingr = read.table(file.path(here::here(), "data/straingr_comb/bw_bifi.straingr.compare.combined.tsv"), sep = "\t", header = F)
bw_bifidobacterium.straingr = straingr_load(bw_bifidobacterium.straingr)
bifidobacterium.ref = read.table(file.path(here::here(), "data/StrainGE_ML/bifidobacterium.ref.tsv"), sep = "\t", header = F)
colnames(bifidobacterium.ref) = c("scaffold", "species")
bifidobacterium.ref[bifidobacterium.ref$species == "Bifi_sp", "species"] = "Bifi_sp."
bw_bifidobacterium.straingr$species =  bifidobacterium.ref[match(bw_bifidobacterium.straingr$scaffold, bifidobacterium.ref$scaffold), "species"] 

#Campylobacter (concateated MAG)
bw_camp.mag.straingr = read.table(file.path(here::here(), "data/straingr_comb/bw_camp_MAG.straingr.compare.combined.tsv"), sep = "\t", header = F)
bw_camp.mag.straingr = straingr_load(bw_camp.mag.straingr)
#As there are only few number of reference genomes detected...(I will manually add those)
camp.mag.ref = data.frame(scaffold = bw_camp.mag.straingr$scaffold %>% unique)
camp.mag.ref$species = c("Camp_infans", "Camp_sp.", "Camp_infans","Camp_jejuni", "Camp_upsaliensis")
bw_camp.mag.straingr$species =  camp.mag.ref[match(bw_camp.mag.straingr$scaffold, camp.mag.ref$scaffold), "species"] 

#Enterococcus
bw_entc.all.straingr = read.table(file.path(here::here(), "data/straingr_comb/bw_entc.straingr.compare.combined.tsv"), sep = "\t", header = F)
bw_entc.all.straingr = straingr_load(bw_entc.all.straingr)
entc.all.ref = read.table(file.path(here::here(), "data/StrainGE_ML/enterococcus.ref.tsv"), sep = "\t", header = F)
colnames(entc.all.ref) = c("scaffold", "species")
bw_entc.all.straingr$species =  entc.all.ref[match(bw_entc.all.straingr$scaffold, entc.all.ref$scaffold), "species"] 

#Klebsiella
bw_kleb.straingr = read.table(file.path(here::here(), "data/straingr_comb/bw_kleb.straingr.compare.combined.tsv"), sep = "\t", header = F)
bw_kleb.straingr = bw_kleb.straingr %>% filter(V1 != "sample1") #remove rows in the middle containing header information
bw_kleb.straingr = straingr_load(bw_kleb.straingr)
klebsiella.ref = read.table(file.path(here::here(), "data/StrainGE_ML/klebsiella.ref.tsv"), sep = "\t", header = F)
colnames(klebsiella.ref) = c("scaffold", "species")
klebsiella.ref[klebsiella.ref$species == "Kleb_sp", "species"] = "Kleb_sp."
bw_kleb.straingr$species =  klebsiella.ref[match(bw_kleb.straingr$scaffold, klebsiella.ref$scaffold), "species"] 

#Enterobacter
bw_entb.straingr = read.table(file.path(here::here(), "data/straingr_comb/bw_entb.straingr.compare.combined.tsv"), sep = "\t", header = F)
bw_entb.straingr = straingr_load(bw_entb.straingr)
entb.ref = read.table(file.path(here::here(), "data/StrainGE_ML/enterobacter.ref.tsv"), sep = "\t", header = F)
colnames(entb.ref) = c("scaffold", "species")
entb.ref[entb.ref$species == "Entb_sp", "species"] = "Entb_sp."
bw_entb.straingr$species =  entb.ref[match(bw_entb.straingr$scaffold, entb.ref$scaffold), "species"] 

#Write into files (comparisons bw studies)
write.table(bw_staingr.ecoli, file.path(here::here(), "data/StrainGE_ML/bw_e.coli.ML.input.tsv"), sep = "\t", row.names = F, quote = F)
write.table(bw_bifidobacterium.straingr, file.path(here::here(), "data/StrainGE_ML/bw_bifidobacterium.ML.input.tsv"), sep = "\t", row.names = F, quote = F)
write.table(bw_bacteroides.straingr, file.path(here::here(), "data/StrainGE_ML/bw_bacteroides.ML.input.tsv"), sep = "\t", row.names = F, quote = F)
write.table(bw_camp.mag.straingr, file.path(here::here(), "data/StrainGE_ML/bw_camp.mag.ML.input.tsv"), sep = "\t", row.names = F, quote = F)
write.table(bw_entc.all.straingr, file.path(here::here(), "data/StrainGE_ML/bw_entc.ML.input.tsv") , sep = "\t", row.names = F, quote = F)
write.table(bw_entb.straingr, file.path(here::here(), "data/StrainGE_ML/bw_entb.ML.input.tsv") , sep = "\t", row.names = F, quote = F)
write.table(bw_kleb.straingr, file.path(here::here(), "data/StrainGE_ML/bw_kleb.ML.input.tsv") , sep = "\t", row.names = F, quote = F)






