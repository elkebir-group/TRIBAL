library(alakazam)
library(dowser)
library(tidyverse)

map_edges <- function(edge_mat, node_map){
  edge_df <- as.data.frame(edge_mat)
  colnames(edge_df) <- c("parent", "child")
  edges <-left_join(edge_df, node_map, by=c("parent"="id")) %>% rename(parent_label=label) %>%
    left_join(node_map, by=c("child"="id")) %>% rename(child_label = label) %>%
    select(-parent, -child) %>%
    rename(parent=parent_label, child=child_label)
  return(edges)
  
}

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


trees <- readRDS(snakemake@input[['trees']])
results_path <- snakemake@params[['outdir']]
if(!dir.exists(results_path)){
    dir.create(results_path)
}
res <- lapply(1:nrow(trees), save_tree, trees, pth =results_path)
saveRDS(res, snakemake@output[['all_files']])


# save_tree(1, trees, outpath)


