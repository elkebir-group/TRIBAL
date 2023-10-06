library(dowser)
library(alakazam)
library(ape)

clones <- readRDS(snakemake@input[['clones']])

clones <- readRDS("/scratch/projects/tribal/experimental_data/GCB_NP_1/igphyml/clones.rds")

dat <- snakemake@wildcards[['dataset']]
script <- snakemake@wildcards[['script']]
mode <- snakemake@wildcards[['mode']]
nproc <- snakemake@params[['nproc']]


# dat <- "day_14"
# script <- "marginal"
# mode <- "refine_ilp"
# nproc <- 10

pth <- sprintf("/scratch/projects/tribal/experimental_data/%s/tribal_recomb/%s/%s/newick", dat, script, mode)

trees <- list()

for(cl in clones$clone_id){
    suffix <- sprintf("%s.nwk.csv", cl)
    fname <- file.path(pth, suffix)
    all_strings <- read.csv(fname)
    newick <- all_strings[1,"newick"]
    tree  <- ape::read.tree(text = newick)
    trees[[cl]] <- tree
}

# clones$trees <- trees
# clones2 <- clones %>% filter(seqs == ntips)
# t <- trees[['B_147_6_76_148_1_41']]
# dat <- clones2$data[[1]]
# seq_ids <- dat@data$sequence_id
# myseqs <-t[['tip.label']]
# setdiff(myseqs, seq_ids)
# setdiff(seq_ids, myseqs)

# ntips <- unlist(lapply(trees, function(x) length(x[['tip.label']])))



# clones$ntips <- ntips



# length(t[['tip.label']])



outtrees <- getTrees(clones, 
            build="igphyml",
            optimize ='lr',
            fixtrees=TRUE,
            exec="/scratch/projects/tribal/igphyml/src/igphyml", 
            nproc=nproc, 
            collapse=FALSE)

saveRDS(outtrees, snakemake@output[['outtrees']])


params <- outtrees$parameters

lhood <- numeric(0)
i <- 1
for (l in params) {

    lhood[i] <- l$lhood 
    i <- i + 1

}

res <- data.frame(clone_id = names(params), likelihood = lhood)
write.csv(res, snakemake@output[['likelihoods']], row.names = FALSE)