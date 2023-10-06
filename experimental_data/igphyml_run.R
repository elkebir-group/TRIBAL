# fname <- "/scratch/projects/tribal/experimental_data/GCB_NP_1/test_igphy/input_data_igphyml-pass.tab"
library(alakazam)
library(dowser)


inputfile <- snakemake@input[['airr']]
nproc <- snakemake@params[['nproc']]
trees_out <- snakemake@output[['trees']]

clones <- readRDS(inputfile)

trees <- getTrees(clones, build="igphyml",optimize ='tlr',
 exec="/scratch/projects/tribal/igphyml/src/igphyml", nproc=nproc, collapse=TRUE)

print("saving trees....")


saveRDS(trees, trees_out)

