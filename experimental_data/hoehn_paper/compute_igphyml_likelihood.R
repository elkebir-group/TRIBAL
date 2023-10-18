library(dowser)
library(alakazam)
library(ape)
library(glue)
library(tidyverse)
clones <- readRDS(snakemake@input[['clones']])




mode <- snakemake@wildcards[['mode']]
nproc <- snakemake@params[['nproc']]
mouse <-  snakemake@wildcards[['mouse']]
nmin <- snakemake@params[['nmin']]


# clones <- readRDS("/scratch/projects/tribal/experimental_data/hoehn_paper/clone_input.rds")

clones <- filter(clones, mouseident==glue("Mouse {mouse}"), seqs >= nmin)

# dat <- "day_14"
# script <- "marginal"
# mode <- "refine_ilp"
# nproc <- 7



trees <- list()

for(cl in clones$clone_id){
    print(cl)

    fname <- glue("Mouse_{mouse}/tribal_ml/{mode}/{cl}/nwk.csv")
    print(fname)
    # fname <- file.path(pth, suffix)
    all_strings <- read.csv(fname, stringsAsFactors=FALSE)

    newick <- all_strings[1,"newick"]

    tree  <- ape::read.tree(text = newick)
    # resolved_tree <-ape:: multi2di(tree)
    trees[[cl]] <- tree
}

clones$trees <- trees


# clones_test <- clones[1,]
# png(file="test.png")
# plot(t,no.margin=TRUE,edge.width=2)
# dev.off()
# t <- trees[['B_79_1_8_64_1_56']]
# dat <- clones$data[[1]]
# seq_ids <- dat@data$sequence_id
# myseqs <-t[['tip.label']]
# setdiff(myseqs, seq_ids)
# setdiff(seq_ids, myseqs)

# ntips <- unlist(lapply(trees, function(x) length(x[['tip.label']])))

# clones2 <- clones 
# clones2$ntips <- ntips
# clones$ntips <- ntips
# filter(clones2, (ntips-1) != seqs)


# length(t[['tip.label']])

# for(cl in bad_list){
#     print(i)
#     clones_temp <- filter(clones, clone_id == cl)
#     temp <-clones_temp$data[[1]]@data
#     clone_seq_ids <- temp$sequence_id
#     tree <- trees[[i]]
#     tip_labels <- tree[['tip.label']]
#     tip_labels <- tip_labels[!str_detect( tip_labels,"GERM")]
#     diff1 <- setdiff(clone_seq_ids, tip_labels)
#     if(length(diff1) > 0){
#         print(diff1)
#     }
#     diff2 <- setdiff(tip_labels, clone_seq_ids)
#     if(length(diff2) > 0){
#         print(diff2)
#     }

#     # if(length(diff1) > 0 | length(diff2) > 0){
#     #     print(clones[i, "clone_id"])
#     # }


# }


# bad_list <- c("B_120_2_8_210_1_17", "B_150_3_6_142_1_6","B_156_1_7_148_1_46","B_82_9_8_148_1_41")
# clones2 <- clones %>% filter(!( clone_id  %in% bad_list) )
outtrees <- getTrees(clones,
                build="igphyml",
                optimize ='lr',
                fixtrees=TRUE,
                exec="/scratch/projects/tribal/igphyml/src/igphyml", 
                nproc=nproc, 
                collapse=FALSE)


# out_fname <- "day_14/tribal_recomb/marginal/refine_ilp/igphyml.trees.rds"
# likelihoods <- "day_14/tribal_recomb/marginal/refine_ilp/likelihoods.csv"
# saveRDS(outtrees, out_fname)
#12
# outtrees <- list()
# for( i in 90:nrow(clones)){
#     print(i)
#     print(clones$clone_id[i])
#     outtrees[[i]] <- getTrees(clones[i, ], 
#                 build="igphyml",
#                 optimize ='lr',
#                 fixtrees=TRUE,
#                 exec="/scratch/projects/tribal/igphyml/src/igphyml", 
#                 nproc=nproc, 
#                 collapse=FALSE)
# }


saveRDS(outtrees, snakemake@output[['outtrees']])



params <- outtrees$parameters

lhood <- numeric(0)
i <- 1
for (l in params) {

    lhood[i] <- l$lhood 
    i <- i + 1

}


res <- data.frame(clone_id = names(params), likelihood = lhood)
# write.csv(res, likelihoods, row.names = FALSE)
write.csv(res, snakemake@output[['likelihoods']], row.names = FALSE)