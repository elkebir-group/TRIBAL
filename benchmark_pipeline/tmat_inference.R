library(tidyverse)
reps <- c(1, 4:6)
cells <- 35
size <- 25
fname <- "inference_error.txt"

pth <- sprintf("/scratch/projects/tribal/bcr-phylo-benchmark/sim_data/replicates/cells%d/size%d", cells, size)

get_errors <- function(r, pth){
 
  error_fname <- file.path(pth, sprintf("rep%d",r),"2.0/0.365/tribal", fname)
  print(error_fname)
  df <- read.csv(error_fname,header=T)

  df$rep <- r

  return(df)
}

error.df <- bind_rows(lapply(reps, get_errors, pth))
error.df$cells <- cells 
error.df$size <- size
error.df <- mutate(error.df, cells =factor(cells), size=factor(size))

error.df$RMSE <- sqrt(error.df$MSE)
head(error.df)

ggplot(error.df, aes(x=cells, y=RMSE)) + geom_boxplot() + facet_wrap(~size, labeller="label_both")
ggplot(error.df, aes(x=cells, y=MAE)) + geom_boxplot() + facet_wrap(~size, labeller="label_both")
