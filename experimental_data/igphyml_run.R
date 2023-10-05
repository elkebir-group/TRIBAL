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


# trees <- readRDS(file.path(pth, dataset, "igphyml/trees.rds"))
results_path <- snakemake@output[['results']]

if(!dir.exists(results_path)){
    dir.create(results_path)
}


# db <- readIgphyml(fname, collapse=TRUE) #, #format="phylo")
# dbphylo <- readIgphyml(fname, collapse=TRUE, format="phylo")
# library(phylobase)
# nodeData(t1)
# tree <- db$trees[[1]]
# t1 <- dbphylo$trees[[1]]



map_nodes <- function(nodes){

  label <- character()
  id <- 1:length(nodes)
  seq <- character()
  for(i in id){
    label[i] <- nodes[[i]]$id
    seq[i] <- nodes[[i]]$sequence
  }
  return(data.frame(id=id, label=label, sequence=seq))
  
}

map_edges <- function(edge_mat, node_map){
  edge_df <- as.data.frame(edge_mat)
  colnames(edge_df) <- c("parent", "child")
  edges <-left_join(edge_df, node_map, by=c("parent"="id")) %>% rename(parent_label=label) %>%
    left_join(node_map, by=c("child"="id")) %>% rename(child_label = label) %>%
    select(-parent, -child) %>%
    rename(parent=parent_label, child=child_label)
  return(edges)
  
}




save_tree <- function(i,tib, pth="igphyml_trees"){
  row <- tib[i,]
  
  tree <- row$trees[[1]]

  data <- row$data[[1]]@data
  node_map <- map_nodes(tree$nodes)
  tree_name <- tree$name
#   mouse <- str_replace(row$mouseident[1], " ", "_")
  write_csv(select(node_map, label, sequence), file.path(pth,  sprintf("%s.sequence.csv", tree_name)), col_names =F)

  tree_edges <- map_edges(tree$edge, node_map %>% select(-sequence))
  write_csv(tree_edges, file.path(pth, sprintf("%s.tree.csv", tree_name)), col_names = F)


#   write_csv(iso.df, file.path(pth, mouse, sprintf("%s.isotypes.csv", tree_name)), col_names = F)
#   write_csv(anno.df, file.path(pth, mouse, sprintf("%s.annotations.csv", tree_name)), col_names = F)
  return(list(tree= tree_edges, sequences =select(node_map, label, sequence )))
}

# outpath <- "/scratch/projects/tribal/experimental_data/GCB_NP_1/test_igphy/results"
# save_tree(1, trees, outpath)

lapply(1:nrow(trees), save_tree, trees, pth =results_path)
