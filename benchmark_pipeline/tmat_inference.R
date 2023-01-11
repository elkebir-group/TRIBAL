source("/scratch/projects/tribal/init.R")
reps <- 1:7
cells <- 35
size <- 25
fname <- "inference_error.txt"

my_theme <- get_theme(base_size)


pth <- sprintf("/scratch/projects/tribal/bcr-phylo-benchmark/sim_data/replicates/cells%d/size%d", cells, size)

get_errors <- function(r, pth){
 
  error_fname <- file.path(pth, sprintf("rep%d",r),"2.0/0.365/tribal", fname)
  print(error_fname)
  df <- read.csv(error_fname,header=T)


  return(df)
}

error.df <- bind_rows(lapply(reps, get_errors, pth))
error.df$cells <- cells 
error.df$size <- size
error.df <- mutate(error.df, cells =factor(cells), size=factor(size))

error.df$RMSE <- sqrt(error.df$MSE)
head(error.df)


entropy.df <- leaf.iso %>% group_by(rep) %>% summarize(entropy = entropy(isotype))
error.df <- inner_join(entropy.df, error.df)
entropy.sum <-entropy_per_clono %>% group_by(rep) %>% summarize(mean_entropy=mean(entropy))

cor(error.df$entropy, error.df$MAE)

ent_corr_plot <-ggplot(error.df, aes(x=entropy, y=MAE)) + geom_point(size=2.5) + geom_smooth(method="lm") +
  xlab("entropy of simulated observed isotypes") + my_theme +
  annotate("text", x = 2.3, y = 0.061, label = "rho == -0.96", size=8,
           parse = TRUE) +  ylab("isotype transition probability MAE")

ent_corr_plot
plot_save(file.path(plot_path, "mae_entropy_corr.pdf"), ent_corr_plot, width=width)

error.df %>% group_by(cells, size) %>% summarize(median_mae = median(MAE), iqr_mae = IQR(MAE))

mae_plot <- ggplot(error.df, aes(x="", y=MAE)) + geom_boxplot() + xlab("experiment") + 
  geom_point(size=2.5)   + my_theme +
  scale_x_discrete(labels="k=25,n=35") +  ylab("isotype transition probability MAE")
plot_save(file.path(plot_path, "mae_plot.pdf"), mae_plot, width=width/3) 



#compare inferred versus simulated state probablities
cell.counts <- iso.sim %>% group_by(rep) %>% count()
iso.sim.dist <- iso.sim  %>% group_by(rep, isotype) %>%  summarize(count= n())
iso.sim.dist <- inner_join(iso.sim.dist, cell.counts)

iso.prop <- iso.sim.dist %>% mutate(prop =count/n)


get_states <- function(r, pth, fname="state_probs.txt"){
  
  state_fname <- file.path(pth, sprintf("rep%d",r),"2.0/0.365/tribal", fname)

  df <- read.csv(state_fname,header=F, col.names=c("prob")) %>% mutate(isotype=row_number()-1)
  
  df$rep <- r
  
  return(df)
}


total.leaf.counts <- leaf.iso %>% group_by(rep) %>% count()
leaf.counts <- leaf.iso %>% group_by(rep, isotype) %>% summarize(count = n())
gt_states <-bind_rows(lapply(reps, get_states, pth))
leaf.prop <- inner_join(leaf.counts, total.leaf.counts) %>% mutate(leaf_prop=count/n) %>% select(-count,-n)
comp.states <- inner_join(gt_states, iso.prop) %>%left_join(leaf.prop)

comp.states.df <- comp.states %>% pivot_longer(cols=contains("prop")) %>%
  mutate(AE=abs(prob-value))


comp.states.sum <- comp.states.df %>% group_by(rep, name) %>% summarize(MAE = sum(AE)/n_distinct(isotype))
iso_dist_comp <- ggplot(comp.states.sum, aes(x=name, y=MAE)) + geom_boxplot() +
  scale_x_discrete(labels=c("observed isotype\nproportions", "inferred isotype\nproportions")) +
  xlab("isotype distribution estimate") + my_theme + ylab("isotype proportion MAE")

comp.states.sum %>% group_by(name) %>% summarize(med = median(MAE))

plot_save(file.path(plot_path, "isotype_dist_comp.pdf"), iso_dist_comp, width=width*2/3)


                                                      