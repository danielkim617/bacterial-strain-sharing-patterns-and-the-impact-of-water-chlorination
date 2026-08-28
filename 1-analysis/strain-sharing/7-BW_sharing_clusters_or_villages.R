# this script is to compare strain sharing rates within vs. between clusters or villages at the household level using the output from the random forest model.
source(file.path(here::here(),"0-config.R"))
#Load data
load(file.path(here::here(), "data/HH.rates.RData"))
###########################################################################################################################
#Comparison of HH level strain-sharing rates within vs. between
get.cluster.test = function(x, n_boot = 1000){
  wit.bet.test = data.frame(matrix(ncol=8)) 
  colnames(wit.bet.test) = c("area", "within_cluster", "within_lower", "within_upper", "between_cluster", "between_lower", "between_upper", "p.value")
  n=1
  for (i in c("urban", "rural")) {
      
    wit.bet.test[n,1] = i
    if(i == "rural"){
      a = x %>% filter(area == i, HH_type == "Between", HH1_location == HH2_location) %>% pull(rate)
      b = x %>% filter(area == i, HH_type == "Between", HH1_location != HH2_location) %>% pull(rate)
    }
    if(i == "urban"){ # In case of urban village category is same as tr
      a = x %>% filter(area == i, HH_type == "Between", HH1_tr == HH2_tr) %>% pull(rate)
      b = x %>% filter(area == i, HH_type == "Between", HH1_tr != HH2_tr) %>% pull(rate)
    }
    a_ci = boot_ci(a, n_boot)
    b_ci = boot_ci(b, n_boot)
    wit.bet.test[n,2] = mean(a)
    wit.bet.test[n,3] = a_ci[1]
    wit.bet.test[n,4] = a_ci[2]
    wit.bet.test[n,5] = mean(b)
    wit.bet.test[n,6] = b_ci[1]
    wit.bet.test[n,7] = b_ci[2]
    wit.bet.test[n,8] = permutation.test(a,b)
    n=n+1
  }
  return(wit.bet.test)
}

#Compare rates (Within vs. Between clusters)
wit.bet.clstr.test.All = get.cluster.test(hh.comb.ML.All)
wit.bet.clstr.test.comm = get.cluster.test(hh.comb.ML.comm)
wit.bet.clstr.test.non_comm = get.cluster.test(hh.comb.ML.non_comm)
wit.bet.clstr.test.Esch = get.cluster.test(hh.comb.ML.Esch)
wit.bet.clstr.test.Bact = get.cluster.test(hh.comb.ML.Bact)
wit.bet.clstr.test.Camp = get.cluster.test(hh.comb.ML.Camp)
wit.bet.clstr.test.Camp.mag = get.cluster.test(hh.comb.ML.Camp.mag)
wit.bet.clstr.test.Bifi = get.cluster.test(hh.comb.ML.Bifi)
wit.bet.clstr.test.Entc = get.cluster.test(hh.comb.ML.Entc)
wit.bet.clstr.test.Entb = get.cluster.test(hh.comb.ML.Entb)
wit.bet.clstr.test.Kleb = get.cluster.test(hh.comb.ML.Kleb)

wit.bet.clstr.test.All$genus = "All"
wit.bet.clstr.test.comm$genus = "comm"
wit.bet.clstr.test.non_comm$genus = "non_comm"
wit.bet.clstr.test.Esch$genus = "Esch"
wit.bet.clstr.test.Bact$genus = "Bact"
wit.bet.clstr.test.Camp$genus = "Camp"
wit.bet.clstr.test.Camp.mag$genus = "Camp_MAG"
wit.bet.clstr.test.Bifi$genus = "Bifi"
wit.bet.clstr.test.Entc$genus = "Entc"
wit.bet.clstr.test.Entb$genus = "Entb"
wit.bet.clstr.test.Kleb$genus = "Kleb"

#stat test (Within vs. Between)
wit.bet.clstr.test = data.frame(rbind(wit.bet.clstr.test.All,wit.bet.clstr.test.comm, wit.bet.clstr.test.non_comm,wit.bet.clstr.test.Esch, wit.bet.clstr.test.Bact, wit.bet.clstr.test.Camp, wit.bet.clstr.test.Camp.mag, wit.bet.clstr.test.Bifi,wit.bet.clstr.test.Entc,wit.bet.clstr.test.Entb, wit.bet.clstr.test.Kleb))
#Save into a file
write.csv(wit.bet.clstr.test, file.path(here::here(), "data/stats/wit.bet.clstr.test.csv"), row.names = F, quote = F)





