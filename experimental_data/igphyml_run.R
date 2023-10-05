# fname <- "/scratch/projects/tribal/experimental_data/GCB_NP_1/test_igphy/input_data_igphyml-pass.tab"
library(alakazam)
library(dowser)
library(tidyverse)
# library(igraph)
# library(ape)

inputfile <- snakemake@input[['airr']]
# pth <- "/scratch/projects/tribal/experimental_data"
# dataset <- "day_14"
# inputfile <- file.path(pth, dataset, "igphyml/input_data.tsv")
dat <- read_tsv(inputfile)
clones <- formatClones(dat, collapse=FALSE)


clones_out <- snakemake@output[['clones']]
saveRDS(clones, clones_out)

nproc <- snakemake@params[['nproc']]
trees <- getTrees(clones, build="igphyml",optimize ='tlr',
 exec="/scratch/projects/tribal/igphyml/src/igphyml", nproc=7, collapse=TRUE)

print("saving trees....")

trees_out <- snakemake@output[['trees']]
saveRDS(trees, trees_out)


