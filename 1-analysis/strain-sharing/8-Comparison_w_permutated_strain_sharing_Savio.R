# This script is to permutate household distances of sample pairs sharing strains to examine the distance-dependent strain-sharing patterns. 
#Run on Savio
source(file.path(here::here(),"0-config.R"))
#Load data
urban.all.pairs = readRDS(file.path(here::here(), "data/urban.all.pairs.rds"))
rural.all.pairs = readRDS(file.path(here::here(), "data/rural.all.pairs.rds"))

# Permutation 1000 iteration - To examine how the graph looks like with the permutated data
# Strain-sharing
library(foreach)
library(doParallel)
# Define the permutation function
permutation_strain <- function(x, min, max) {
  # Create a cluster with the number of cores available
  num_cores <- detectCores() - 1  # Save one core for the OS
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Export necessary objects to the cluster
  clusterExport(cl, varlist = c("x", "min", "max", "required_packages"), envir = environment())
  
  # Parallel execution using foreach
  b <- foreach(i = 1:1000, .combine = rbind, .packages = "foreach") %dopar% {
    # Load necessary packages within each worker
    for (pkg in required_packages) {
      if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
      }
    }
    
    # Source the configuration script within each worker
    source(file.path(here::here(), "0-functions.R"))
    
    # Perform the permutation
    a <- x
    set.seed(617 + i)
    for (g in c("Bact.share", "Bifi.share", "Esch.share", "Camp_MAG.share", "Entc.share", "Entb.share", "Kleb.share", "All.share", "commensal.share", "non_commensal.share")) {
      a[, all_of(g)] <- sample(a[, all_of(g)], replace = FALSE)
    }
    perm_a <- get_sharing_dist(a, min, max)
    perm_a$iter <- i
    perm_a$tr = "All"
    
    perm_w <- get_sharing_dist(a %>% filter(tr_type == "Water"), min, max)
    perm_w$iter <- i
    perm_w$tr = "Water"
    
    perm_c <- get_sharing_dist(a %>% filter(tr_type == "Control"), min, max)
    perm_c$iter <- i
    perm_c$tr = "Control"
    
    perm = rbind(perm_a, perm_w, perm_c)
    perm
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # Return the combined result
  return(data.frame(b))
}

# Run permutation
urban_strain_perm = permutation_strain(urban.all.pairs, 50, 1500)

rural_strain_perm = permutation_strain(rural.all.pairs, 50, 1500)

#Save into files
saveRDS(urban_strain_perm, file.path(here::here(), "data/perm_output_all/urban_strain_perm.rds"))

saveRDS(rural_strain_perm, file.path(here::here(), "data/perm_output_all/rural_strain_perm.rds"))






