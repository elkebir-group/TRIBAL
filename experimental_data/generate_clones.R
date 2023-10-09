library(alakazam)
library(dowser)
library(tidyverse)

inputfile <- snakemake@input[['airr']]
mappingfile <- snakemake@input[['mapping']]
clones_out <- snakemake@output[['clones']]

# pth <- "/scratch/projects/tribal/experimental_data"
# dataset <- "day_14"
# inputfile <- file.path(pth, dataset, "igphyml/input_data.tsv")
# mappingfile <- file.path(pth, dataset, "igphyml/name_isotype_mapping.csv")

in_dat <- read_tsv(inputfile)

mapping <- read_csv(mappingfile) %>% 
    select(sequence_id, clone_id)
dat <- inner_join(in_dat, mapping)

clones <- formatClones(dat, collapse=FALSE, verbose=TRUE)

saveRDS(clones, clones_out)
