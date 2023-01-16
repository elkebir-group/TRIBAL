source("/scratch/projects/tribal/init.R")
reps <- 1:7
cells <- 35
size <- 25
fname <- "inference_error.txt"

my_theme <- get_theme(base_size)



get_errors <- function(r, pth, fname){
 
  error_fname <- file.path(pth, sprintf("rep%d",r),"2.0/0.365/tribal", fname)

  if(file.exists(error_fname)){
    df <- read.csv(error_fname,header=T)
    df$rep <- r
  }else{
    print(sprintf("%s does not exist!!",error_fname))
    df <- data.frame()
  }



  return(df)
}

error.df <- data.frame()
size <- 25
for(c in cell_values){
  pth <- sprintf("/scratch/projects/tribal/bcr-phylo-benchmark/sim_data/replicates/cells%d/size%d", c, size)
  print(pth)
  e.df <- bind_rows(lapply(reps, get_errors, pth, "inference_error.txt"))
  if(nrow(e.df) >0){
    e.df$cells <- c 
    e.df$size <- size
    error.df <- bind_rows(error.df, e.df)
  }

}


error.df <- mutate(error.df,
                   RMSE = sqrt(MSE))





entropy.df <- leaf.iso %>% group_by(cells, rep) %>% summarize(entropy = entropy(isotype)) %>% ungroup()
error.df <- inner_join(entropy.df, error.df, by=c("cells", "rep")) %>% inner_join(combos_to_use)


error.df %>% group_by(cells) %>% summarize(cor = cor(entropy, MAE))
cor(error.df$entropy, error.df$MAE)

ent_corr_plot <-ggplot(error.df, aes(x=entropy, y=MAE, color=factor(cells))) + 
  geom_point(size=3) + geom_smooth(method="lm", se=F, linetype="dashed") +
  xlab("entropy of simulated observed isotypes") + my_theme +
  scale_color_discrete(name="cells") #+
  # annotate("text", x = 2.0, y = 0.06, label = "rho == -0.95", size=8,
  #          parse = TRUE) +  ylab("isotype transition probability MAE") +
  # theme(legend.position=c(.775,.86))

ent_corr_plot
plot_save(file.path(plot_path, "mae_entropy_corr.pdf"), ent_corr_plot, width=width*2/3)

error.df %>% group_by(cells, size) %>% summarize(median_mae = median(MAE), iqr_mae = IQR(MAE))

mae_plot <- ggplot(error.df, aes(x=factor(cells), y=MAE)) + geom_boxplot() + xlab("experiment") + 
  geom_point(size=2.5)   + my_theme +
  scale_x_discrete(name="cells") +  ylab("isotype transition probability MAE")
mae_plot
plot_save(file.path(plot_path, "mae_plot.pdf"), mae_plot, width=width/2) 



#compare inferred versus simulated state probablities
cell.counts <- iso.sim %>% group_by(rep) %>% count()
iso.sim.dist <- iso.sim  %>% group_by(rep, isotype) %>%  summarize(count= n()) 
iso.sim.dist <- inner_join(iso.sim.dist, cell.counts)

iso.prop <- iso.sim.dist %>% mutate(prop =count/n)


get_states <- function(r, pth, cells, fname="state_probs.txt"){
  pth <- sprintf("/scratch/projects/tribal/bcr-phylo-benchmark/sim_data/replicates/cells%d/size25", cells)
  state_fname <- file.path(pth, sprintf("rep%d",r),"2.0/0.365/tribal", fname)

  df <- read.csv(state_fname,header=F, col.names=c("prob")) %>% mutate(isotype=row_number()-1)
  
  df$rep <- r
  df$cells <- cells
  
  return(df)
}


total.leaf.counts <- leaf.iso %>% group_by(cells, rep) %>% count()
leaf.counts <- leaf.iso %>% group_by(cells, rep, isotype) %>% summarize(count = n())
gt_states.35 <-bind_rows(lapply(reps, get_states, pth, 35))
gt_states.65 <-bind_rows(lapply(c(1,2,4,5), get_states, pth, 65))
gt_states <- bind_rows(gt_states.35, gt_states.65)
leaf.prop <- inner_join(leaf.counts, total.leaf.counts) %>% mutate(leaf_prop=count/n) %>% select(-count,-n)
comp.states <- inner_join(gt_states, iso.prop) %>%left_join(leaf.prop)

comp.states.df <- comp.states %>% pivot_longer(cols=contains("prop")) %>%
  mutate(AE=abs(prob-value))


comp.states.sum <- comp.states.df %>% group_by(cells,rep, name) %>% 
  summarize(MAE = sum(AE)/n_distinct(isotype)) %>%
  inner_join(combos_to_use)
iso_dist_comp <- ggplot(comp.states.sum, aes(x=factor(cells), y=MAE, fill=name)) + geom_boxplot() +
  scale_fill_discrete(name="", labels=c("obs.cells", "TRIBAL")) +
  xlab("cells") + my_theme + 
  # theme(legend.position=c(.7,.86))+
  ylab("isotype proportion MAE")# +
  # facet_wrap(~cells)
iso_dist_comp
comp.states.sum %>% group_by(name) %>% summarize(med = median(MAE))

plot_save(file.path(plot_path, "isotype_dist_comp.pdf"), iso_dist_comp, width=width*0.5)


                                                      