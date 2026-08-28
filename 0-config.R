#config file
rm(list = ls()) # Clear the environment

# Install missing packages
packages <- c("reshape2", "parallel", "readxl", "ggplot2", "ggpubr", "ggthemes", "ape",
              "stringr", "RColorBrewer", "viridis", "tidyr", "unpivotr", "tidyverse", 
              "vegan", "mclust", "geosphere", "Rtsne", "corncob", "corrplot", "magrittr", 
              "biomformat", "tibble", "igraph", "tnet", "haven", "minpack.lm", "pairwiseAdonis",
              "foreach", "doParallel", "uwot", "ggforce", "concaveman", "dplyr", "binom", "patchwork")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

if (!requireNamespace("phyloseq", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("phyloseq")
    remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}

if (!requireNamespace("ggtree", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("ggtree")
}

#Load libraries
library(reshape2)
library(parallel)
library(readxl)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(stringr)
library(RColorBrewer)
library(viridis)
library(tidyr)
library(unpivotr)
library(tidyverse)
library(vegan)
library(mclust)
library(geosphere)
library(phyloseq)
library(Rtsne)
library(pairwiseAdonis)
library(corncob)
library(corrplot)
#library(microbiome)
#library(DirichletMultinomial)
library(magrittr)
library(parallel)
library(biomformat)
library(tibble)
#library(mgcv)
library(igraph)
library(tnet)
library(haven)
library(minpack.lm)
library(foreach)
library(doParallel)
library(uwot)
library(ggforce)
library(concaveman)
library(dplyr)
library(binom)
library(patchwork)
library(ggtree)
library(ape)

set.seed(617)

#Loading the script where functions are located
source(file.path(here::here(), "0-functions.R"))
