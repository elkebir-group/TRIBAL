library(tidyverse)
reps <- 1:6
cells <- 35
size <- 25
methods <- c("tribal_score", "tribal_score_gt", "tribal_refine", "tribal_refine_gt", "validaggreg")

get_results <- function(r, methods ){
  rep_res <- data.frame()
  for(m in methods){
    pth <- sprintf("/scratch/projects/tribal/bcr-phylo-benchmark/sim_data/replicates/cells%d/size%d/rep%d/2.0/0.365", cells, size, r)
    fname <- sprintf('%s.tsv', m)
    print(fname)
    df <- read.table(file.path(pth,fname), header=T)
    if('rep' %in% colnames(df)){
      df <- df %>% rename(clonotype=rep)
    }else{
      df <- df %>% group_by(method) %>% mutate(clonotype=row_number()) %>% ungroup()
    }
    
    rep_res <- bind_rows(rep_res, df)

  }
  rep_res$rep <- r
  return(rep_res)
  
    
}

get_isotypes <- function(r, size){
  rep_res <- data.frame()
  for(i in 1:size){
    pth <- sprintf("/scratch/projects/tribal/bcr-phylo-benchmark/sim_data/replicates/cells%d/size%d/rep%d/2.0/0.365/%d", cells, size, r,i)
    fname <- "GCsim.isotypes"

    df <- read.csv(file.path(pth,fname), header=T, col.names = c("cell", "isotype"))
    df$clonotype <- i

    
    rep_res <- bind_rows(rep_res, df)
    
  }
  rep_res$rep <- r
  return(rep_res)
  
  
}

entropy <- function(x,states=0:6){
  n <- length(x)

  
  iso_prob <- data.frame(iso=x) %>% group_by(iso) %>% summarize(p_x = n()/n)

  for(s in states){
    if(!(s %in% x)){
      iso_prob <- bind_rows(data.frame(iso=s, p_x=0.001),iso_prob)
    }
  }
  iso_probs<- iso_prob %>% group_by(iso) %>% summarize(p_x = p_x/sum(p_x))

  
  entropy_val <- -1*sum(iso_prob$p_x*log2(iso_prob$p_x))
  return(entropy_val)
  
  
}

entropy(test)

arrange(entropy_per_clono,entropy)

leaf.iso <- bind_rows(lapply(reps, get_isotypes, size))

filter(leaf.iso, clonotype==10, rep==6)

entropy_per_clono <- leaf.iso %>% group_by(rep, clonotype) %>% summarize(entropy = entropy(isotype))
ggplot(entropy_per_rep, aes(x=factor(rep), y=entropy)) + geom_boxplot() +
  xlab("Replication") + ylab("Entropy per clonotype")

tree.results <- bind_rows(lapply(reps, get_results, methods)) %>%
  filter(RF < 999) %>%
  mutate(clonotype= factor(clonotype), rep=factor(rep),
         alpha= ifelse(is.na(alpha), 1.0, alpha))

  
head(tree.results)

cand.df <- read.csv("/scratch/projects/tribal/bcr-phylo-benchmark/sim_data/replicates/cells35/size25/candidate_counts.csv") %>%
  mutate(rep=factor(rep), clonotype=factor(clonotype))
cand.df %>% group_by(rep) %>% summarize(med= median(ncand), max = max(ncand), quant3 =quantile(ncand, 0.75)) 
head(cand.df)
ggplot(cand.df, aes(x=factor(rep), y=ncand)) + geom_boxplot() + xlab("Replication") +
  scale_y_log10()
ggplot(tree.results, aes(x=method, fill=factor(alpha), y=RF)) + 
  geom_boxplot() + 
  # facet_wrap(~rep, labeller="label_both") +
  vtext +
  scale_fill_discrete(name="lambda")




ggplot(tree.results, aes(x=method, fill=factor(alpha), y=log(MRCA))) + 
  geom_boxplot() + 
  # facet_wrap(~rep, labeller="label_both") +
  vtext +
  ylab("log MRCA") +
  scale_fill_discrete(name="lambda")

ggplot(tree.results %>% filter(str_detect(method, "tribal")), aes(x=method, fill=factor(alpha), y=MRCA_ISO)) + 
  geom_boxplot() + 
  facet_wrap(~rep, labeller="label_both") +
  vtext +
  ylab("MRCA_ISO") +
  scale_fill_discrete(name="lambda")

head(tree.results)

filter(tree.results, str_detect( method,'tribal_search')) %>% ggplot(aes(x=MRCA_ISO, y=RF)) +
  geom_point() +
  facet_grid(alpha ~rep, labeller="label_both") +
  geom_smooth(method="lm") +
  ggtitle("tribal_search")
filter(tree.results, str_detect( method,'tribal_search')) %>% ggplot(aes(x=MRCA_ISO, y=log(MRCA))) +
  geom_point() +
  facet_grid(alpha ~rep, labeller="label_both") +
  geom_smooth(method="lm") +
  ggtitle("tribal_search")
  

head(tree.results)
tree.wide.rf <- select(tree.results, clonotype, method, alpha, RF, rep, MRCA_ISO) %>% 
  pivot_wider(id_cols =c("clonotype","rep"), names_from=c("method", "alpha"), values_from=c("RF", "MRCA_ISO"))


ggplot(tree.wide.rf, aes(y=RF_dnapars_1-RF_tribal_refine_0.75 , x=MRCA_ISO_tribal_refine_0.75)) +
  geom_point() + 
  facet_wrap(~rep, labeller="label_both") + geom_smooth(method="lm")

ggplot(tree.wide.rf, aes(x=RF_tribal_refine_0.75, y=RF_dnapars_1, color=(RF_tribal_refine_0.75 <= RF_dnapars_1))) +
  geom_point() + scale_color_discrete(name="tribal better") #+
  # facet_wrap(~rep, labeller="label_both")

ggplot(tree.wide.rf, aes(x=log(MRCA_tribal_refine_0.75), y=log(MRCA_dnapars_1), 
                         color=(RF_tribal_refine_0.75 <= RF_dnapars_1))) +
  geom_point() + scale_color_discrete(name="tribal better")

# ggplot(tree.wide.rf, aes(x=tribal_refine_0.75, y=dnapars_1, color=(tribal_refine_0.75 <= dnapars_1))) +
#   geom_point() + scale_color_discrete(name="tribal better") +
#   facet_wrap(~rep, labeller="label_both")


ggplot(tree.wide.rf, aes(x=tribal_search_0.95, y=IgPhyML_1, color=(tribal_search_0.95 <= IgPhyML_1))) +
  geom_point() + scale_color_discrete(name="tribal<=") +
  facet_wrap(~rep, labeller="label_both")



entropy_per_clono <- entropy_per_clono %>% ungroup() %>% mutate(rep = factor(rep), clonotype=factor(clonotype))
tree.wide.rf <- inner_join(tree.wide.rf, entropy_per_clono)

tree.wide.rf <- tree.wide.rf %>% mutate(dpar_ts0.95_diff = dnapars_1-tribal_search_0.95, 
                                        igphy_ts0.95_diff=IgPhyML_1- tribal_search_0.95)

tree.wide.rf <- inner_join(tree.wide.rf, cand.df)

ggplot(tree.wide.rf, aes(x=ncand, y=dpar_ts0.95_diff)) + geom_point(aes(color=rep)) +
  geom_smooth(method="lm", aes(color=rep, group=rep), se=F) + #scale_x_log10() +
  facet_wrap(~rep, scales="free")

ggplot(tree.wide.rf, aes(x=ncand, y=dnapars_1-tribal_refine_0.75)) + geom_point(aes(color=rep)) +
  geom_smooth(method="lm", aes(color=rep, group=rep), se=F) + scale_x_log10() +
  facet_wrap(~rep, scales="free")

ggplot(tree.wide.rf, aes(x=entropy, y=dpar_ts0.95_diff)) + geom_point(aes(color=rep)) +
  geom_smooth(method="lm", aes(color=rep, group=rep), se=F)


ggplot(tree.wide.mrca, aes(x=log(tribal_search_0.95), y=log(dnapars_1), color=(tribal_search_0.95 <= dnapars_1))) +
  geom_point() + scale_color_discrete(name="tribal better") +
  facet_wrap(~rep, labeller="label_both") +
  ggtitle("MRCA")




ggplot(tree.wide.mrca, aes(x=log(tribal_search_0.95), y=log(IgPhyML_1), color=(tribal_search_0.95 <= IgPhyML_1))) +
  geom_point() + scale_color_discrete(name="tribal <=") +
  facet_wrap(~rep, labeller="label_both", scales="free") +
  ggtitle("MRCA")
