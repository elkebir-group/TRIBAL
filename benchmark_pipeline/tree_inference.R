




source("/scratch/projects/tribal/init.R")

#get benchmarking methods
iso_res <- data.frame()
for(c in cell_values){
  for(m in other_methods){
    fname <- sprintf("/scratch/projects/tribal/benchmark_pipeline/sim_data/replicates/cells%d/size25/%s.tribal.tsv", c, m)
    print(fname) 
    if(file.exists(fname)){
      df <- read.table(fname, header=T)
      df$cells <- c
      iso_res <- bind_rows(iso_res, df)
    }

     }

}

cand.counts.35 <- read.csv(sprintf("/scratch/projects/tribal/benchmark_pipeline/sim_data/replicates/cells%d/size25/candidate_counts.csv", 35, 25)) %>%
  mutate(cells = 35)
cand.counts.65 <- read.csv(sprintf("/scratch/projects/tribal/benchmark_pipeline/sim_data/replicates/cells%d/size25/candidate_counts.csv", 65, 25)) %>%
  mutate(cells = 65)
cand.counts <- bind_rows(cand.counts.35, cand.counts.65) %>% inner_join(combos_to_use)

cand.counts %>% group_by(cells) %>% summarise(median=median(ncand), iqr= IQR(ncand), mean=mean(ncand), sd=sd(ncand))
head(cand.counts)

cand.tree.dist <- ggplot(cand.counts %>% inner_join(combos_to_use), aes(x=ncand)) +
  geom_histogram(bins=50) + facet_wrap(~cells, labeller="label_both", scales="free") +
  xlab("parsimony forest size") + my_theme

plot_save(file.path(plot_path, "forest_size.pdf"), cand.tree.dist)

#ggplot(cand.counts, aes(x=factor(cells), y=ncand)) +  geom_boxplot() + scale_y_log10() +methods <- c("tribal_score", "tribal_score_gt", "tribal_search",  "tribal_search_gt", "tribal_refine", "tribal_refine_gt", "validaggreg")
#methods <- c("tribal_score", "tribal_search", "tribal_refine")
methods <- c("tribal_search")

get_results <- function(r, methods, cells, size ){
  rep_res <- data.frame()
  for(m in methods){
    pth <- sprintf("/scratch/projects/tribal/bcr-phylo-benchmark/sim_data/replicates/cells%d/size%d/rep%d/2.0/0.365", cells, size, r)
    fname <- sprintf('%s.tsv', m)
    if(file.exists(file.path(pth,fname))){
      print(file.path(pth,fname))
      df <- read.table(file.path(pth,fname), header=T)
      if('rep' %in% colnames(df)){
        df <- df %>% rename(clonotype=rep)
      }else{
        df <- df %>% group_by(method) %>% mutate(clonotype=row_number()) %>% ungroup()
      }
      
      
      rep_res <- bind_rows(rep_res, df)
    }
    
  }
      if(nrow(rep_res) > 0){
        rep_res$rep <- r
        rep_res$cells <- cells
        rep_res$size <- size
        
      }
  
      return(rep_res)
      

  
    
}




tree.results <- data.frame()

for(c in cell_values){
  size <- 25
  res <- bind_rows(lapply(reps, get_results, methods,c, size)) 
  tree.results <- bind_rows(tree.results, res)
  
}

#t.res <- filter(tree.results,  alpha==0.75, method=="tribal_search") %>% mutate(method="TRIBAL")
t.res <- tree.results

tree.results <- bind_rows(t.res, iso_res) %>% 
#  select(-alpha) %>%
  mutate(method =as.character(method)) %>% inner_join(combos_to_use)

# dnaml1 <- filter(tree.results, method=="dnaml") %>% select(-RF, -MRCA_ISO,  -N_taxa) 
# dnaml2 <- filter(tree.results, method=="dnaml_tribal") %>% select(-MRCA,  -N_taxa) %>%mutate(method="dnaml", size=25)
# dnaml <- inner_join(dnaml2, dnaml1)
# 
# igphy1 <- filter(tree.results, method=="igphyml") %>% select(-RF, -MRCA_ISO,  -N_taxa) 
# igphy2 <- filter(tree.results, method=="igphyml_tribal") %>% select(-MRCA,  -N_taxa) %>%mutate(method="dnaml", size=25)
# igphy<- inner_join(igphy1, igphy2)
# 
# tree.results <- filter(tree.results, !(str_detect(method, "dnaml"))) %>% bind_rows(dnaml)

tree.results.long <- tree.results %>% 
  pivot_longer( cols=c("RF", "MRCA_ISO", "MRCA")) %>%
  filter(value < 99)


res.sum <- tree.results.long %>% group_by(cells, method, name) %>%
  summarize(mean=mean(value), sd=sd(value), median=median(value), iqr=IQR(value))
View(res.sum %>% ungroup() %>% arrange(cells,name, median))

tree.results.long <- tree.results.long %>%
  mutate(
    metric = factor(name, levels=c("RF", "MRCA", "MRCA_ISO"),
                         labels = c("RF\ndistance", "sequence\nMRCA","isotype\nMRCA"),
                         ordered=T),
    method_labs = factor(method, levels=c("tribal_search","igphyml", "dnapars","gctree_isotype", "dnaml" ),
                    labels = c("TRIBAL",  "IgPhyML", "dnapars", "GCTree","dnaml" )))

cells.65.plot <- ggplot(tree.results.long %>% filter(cells==65), aes(x=method_labs, y=value, fill=method_labs)) + geom_boxplot() + 
  facet_wrap(~metric, scales="free_y") + vertx_theme + ylab("")  + xlab("method") +
  scale_fill_manual(values=method_cols) + theme(legend.position="none")

plot_save(file.path(plot_path, "cells.65.tree.metrics.pdf"), cells.65.plot, width=0.73*width, height=height)


cells.65.plot

 ggplot(tree.results.long %>% filter(cells==65), aes(x=method, y=value, fill=factor(alpha))) + geom_boxplot() + 
  facet_wrap(~name, scales="free_y") + vertx_theme + ylab("")  + xlab("method") +
   ggtitle("65 cells") + get_theme_vertx(base_size) + scale_fill_discrete(name="lambda")
 
 ggplot(tree.results.long %>% filter(cells==35), aes(x=method, y=value,fill=factor(alpha))) + geom_boxplot() + 
   facet_wrap(~name, scales="free_y") + vertx_theme + ylab("")  + xlab("method") +
   ggtitle("35 cells") + get_theme_vertx(base_size) + scale_fill_discrete(name="lambda")

cells.65.large <-ggplot(tree.results.long %>% filter(cells==65) %>%
         inner_join(large_forest),
       aes(x=method_labs, y=value, fill=method_labs)) + geom_boxplot() + 
  facet_wrap(~metric, scales="free_y") + vertx_theme + ylab("")  + xlab("method") +
  scale_fill_manual(values=method_cols) +  theme(legend.position="none")

plot_save(file.path(plot_path, "cells.65.tree.large.metrics.pdf"), cells.65.large, width=width, height=height)

cells.35.plot <- ggplot(tree.results.long %>% filter(cells==35), aes(x=method_labs, y=value)) + geom_boxplot() + 
  facet_wrap(~metric, scales="free_y") + vertx_theme + ylab("")  + xlab("method")

plot_save(file.path(plot_path, "cells.35.tree.metrics.pdf"), cells.35.plot)



cells.35.large <-ggplot(tree.results.long %>% filter(cells==35) %>%
                          inner_join(large_forest),
                        aes(x=method_labs, y=value)) + geom_boxplot() + 
  facet_wrap(~metric, scales="free_y") + vertx_theme + ylab("")  + xlab("method")
plot_save(file.path(plot_path, "cells.35.tree.large.metrics.pdf"), cells.35.large)

#large_forest <- cand.counts %>% group_by(cells) %>% filter( ncand >= quantile(ncand, 0.75) )
#large_forest <- cand.counts %>% group_by(cells) %>% filter( ncand >=4 )
large_forest <- cand.counts %>% group_by(cells) %>% filter( ncand >= quantile(ncand, 0.5) )

ggplot(tree.results.long %>% filter(cells==65) %>% inner_join(large_forest), aes(x=method_labs, y=value)) + geom_boxplot() + 
  facet_wrap(~metric, scales="free_y") + vertx_theme + ylab("")  + xlab("method")

sim.entropy <- iso.sim %>% inner_join(combos_to_use) %>% group_by(cells, rep, clonotype) %>%
  summarize(entropy = entropy(isotype))
sim.entropy <- leaf.iso %>% inner_join(combos_to_use) %>% group_by(cells, rep, clonotype) %>%
  summarize(entropy = entropy(isotype))

ggplot(sim.entropy, aes(x=entropy)) + geom_histogram() + facet_wrap(~cells) 
large.ent <- sim.entropy %>% group_by(cells) %>% filter(entropy >quantile(entropy, 0.5))
ggplot(tree.results.long %>% filter(cells==65, method != "tribal_score") %>%
         inner_join(large_forest),
       aes(x=method, y=value)) + geom_boxplot() + 
  facet_wrap(~metric, scales="free_y") + vtext + ylab("") +
  ggtitle("65 cells per clonotype, filtered clonotypes")


larg.sum <- tree.results.long %>%  inner_join(large_forest) %>% group_by(cells, method,metric) %>%
  summarize(median = median(value))
View(larg.sum)
ggplot(tree.results.long %>% filter(cells==65) %>%
         inner_join(large_forest),
       aes(x=method, y=value)) + geom_boxplot() + 
  facet_wrap(~metric, scales="free_y") + vtext + ylab("") +
  ggtitle("65 cells per clonotype, filtered clonotypes")

ggplot(tree.results.long %>% filter(cells==35, method != "tribal_score") %>%
         inner_join(large_forest) %>% filter(ncand > 5),
       aes(x=method, y=value)) + geom_boxplot() + 
  facet_wrap(~metric, scales="free_y") + vtext + ylab("") 



foo <- tree.results.long %>% filter(cells==65, method != "tribal_score") %>%
  inner_join(large.ent) %>% inner_join(large_forest)
ggplot(tree.results.long %>% filter(cells==65, method != "tribal_score") %>%
         inner_join(large.ent) %>% 
         inner_join(large_forest),
       aes(x=method, y=value)) + geom_boxplot() + 
  facet_wrap(~metric, scales="free_y") + vtext + ylab("") #+
  # ggtitle("35 cells per clonotype, filtered clonotypes")

ggplot(tree.results.long %>% inner_join(large_forest), 
       aes(x=method, y=value)) + geom_boxplot() + 
  facet_grid(metric~cells, scales="free_y") + vtext 

# +
#   ggtitle("35 cells per clonotype, filtered clonotypes")

# ggplot(tree.results.long %>% filter(cells==35, rep %in% c(7,6,5,4,3)), aes(x=method, y=value)) + geom_boxplot() + 
#   facet_wrap(~name, scales="free_y") + vtext

ggplot(tree.results %>% filter(MRCA < 99), aes(x=method, y=MRCA)) + 
  geom_boxplot() + facet_wrap(~cells) + vtext 

ggplot(tree.results %>% filter(RF < 1000), aes(x=method, y=RF)) + 
  geom_boxplot() + facet_wrap(~cells) + vtext

ggplot(tree.results %>% filter(MRCA_ISO < 1000), aes(x=method, y=MRCA_ISO)) + 
  geom_boxplot() + facet_wrap(~cells, labeller="label_both") + vtext

# 
# dnapars_iso <- filter(tree.results, method=="tribal_score_gt") %>% select(-MRCA,-RF) %>% mutate(method="dnapars") %>%
#   pivot_longer(cols=c("MRCA_ISO"))
# 
# 
# tree.result.long <- tree.results %>% filter(method=="tribal_search"| !str_detect(method, "tribal")) %>% 
#   pivot_longer(cols=c("MRCA", "MRCA_ISO", "RF")) 
# 
#   
# tree.result.long <- tree.result.long %>% mutate(alpha= ifelse(lambda==1, 0.95,lambda),
#                                                 method=ifelse(method=="tribal_search", "TRIBAL", as.character(method))) %>%
#   filter(lambda==0.95) %>%
#  # bind_rows(dnapars_iso %>% select(-alpha)) %>%
#tree.results.long <-mutate(tree.results.long ,metric=ifelse(name=="MRCA_ISO", "isotype\nMRCA", ifelse(name=="MRCA", "sequence\nMRCA", "RF")))
metric_box_plot <- ggplot(tree.results.long, aes(x=method, y=value, fill=factor(lambda))) + geom_boxplot() +
  facet_grid(metric~cells, scales="free_y")  + xlab("method") + ylab("") +
  theme(legend.position="top") + vtext + scale_fill_discrete(name="lambda")


metric_box_plot
plot_save(file.path(plot_path,"metric_box_plot.pdf"), metric_box_plot)
tree.result.long %>% group_by(method, name) %>% summarize(mean= mean(value, na.rm=T), sd= sd(value, na.rm=T))








head(tree.results)
tree.wide.rf <- tree.results %>%
  filter(RF < 1000, MRCA < 1000, MRCA_ISO < 1000) %>%
  select( clonotype, cells, method, RF, rep, MRCA_ISO, MRCA) %>% 
  pivot_wider(id_cols =c("clonotype","rep", "cells"), names_from=c("method"), values_from=c("RF", "MRCA", "MRCA_ISO"))


# ggplot(tree.wide.rf, aes(y=RF_dnapars_1-RF_tribal_refine_0.75 , x=MRCA_ISO_tribal_refine_0.75)) +
#   geom_point() + 
#   facet_wrap(~rep, labeller="label_both") + geom_smooth(method="lm")
# 
# ggplot(tree.wide.rf, aes(x=RF_tribal_refine_0.75, y=RF_dnapars_1, color=(RF_tribal_refine_0.75 <= RF_dnapars_1))) +
#   geom_point() + scale_color_discrete(name="tribal better") #+
#   # facet_wrap(~rep, labeller="label_both")

ggplot(tree.wide.rf, aes(x=RF_TRIBAL, y=RF_igphyml, color=(RF_TRIBAL <=RF_igphyml))) +
  geom_point() + scale_color_discrete(name="tribal better") + facet_wrap(~cells, scales="free_y") +
  ylim(c(0,25)) + my_theme + xlab("TRIBAL RF") + ylab("IgPhyML RF")

ggplot(tree.wide.rf, aes(x=MRCA_TRIBAL, y=MRCA_igphyml, color=(MRCA_TRIBAL <=MRCA_igphyml))) +
  geom_point() + scale_color_discrete(name="tribal better") + facet_wrap(~cells, scales="free_y") +
  my_theme + xlab("TRIBAL MRCA") + ylab("IgPhyML MRCA")


ggplot(tree.wide.rf, aes(x=RF_TRIBAL, y=RF_gctree_isotype, color=(RF_TRIBAL <=RF_gctree_isotype))) +
  geom_point() + scale_color_discrete(name="tribal better") + facet_wrap(~cells, scales="free_y") +
  ylim(c(0,25)) + my_theme + xlab("TRIBAL RF") + ylab("GCTree RF")

ggplot(tree.wide.rf, aes(x=MRCA_TRIBAL, y=MRCA_gctree, color=(MRCA_TRIBAL <=MRCA_igphyml))) +
  geom_point() + scale_color_discrete(name="tribal better") + facet_wrap(~cells, scales="free_y") +
  my_theme + xlab("TRIBAL MRCA") + ylab("IgPhyML MRCA")


ggplot(tree.wide.rf, aes(x=MRCA_TRIBAL, y=MRCA_dnapars_1, 
                         color=(MRCA_TRIBAL <= MRCA_dnapars_1))) +
   geom_point() + scale_color_discrete(name="TRIBAL better") + 
  facet_wrap(~cells, scales="free_y", nrow=2)

ggplot(tree.wide.rf, aes(x=RF_tribal_search_0.75, y=RF_dnapars_1, 
                         color=(RF_tribal_search_0.75 <= RF_dnapars_1))) +
  geom_point() + scale_color_discrete(name="TRIBAL better") + 
  facet_wrap(~cells, scales="free_y", nrow=2)

tree.wide.rf <- mutate(tree.wide.rf, RF_dnapars_tribal = RF_dnapars_1 - RF_tribal_search_0.95)
tree.wide.rf <- mutate(tree.wide.rf, MRCA_dnapars_tribal = MRCA_dnapars_1 - MRCA_tribal_search_0.95)
ggplot(tree.wide.rf, aes(x=factor(cells), y=RF_dnapars_tribal)) + geom_jitter()
ggplot(tree.wide.rf, aes(x=factor(cells), y=MRCA_dnapars_tribal)) + geom_jitter()
ggplot(tree.wide.rf, aes(x=RF_tribal_refine_0.75, y=RF_dnapars_1, 
                         color=(RF_tribal_refine_0.75 <= RF_dnapars_1))) +
  geom_point() + scale_color_discrete(name="TRIBAL better") + 
  facet_wrap(~cells, scales="free_y")

ggplot(tree.wide.rf, aes(x=RF_tribal_search_0.95, y=RF_dnapars_1, 
                         color=(RF_tribal_search_0.95 <= RF_dnapars_1))) +
  geom_point() + scale_color_discrete(name="TRIBAL better") + 
  facet_wrap(~cells, scales="free_y")

igphyml_mrca_comp <- ggplot(tree.wide.rf, aes(x=MRCA_tribal_search_0.75, y=MRCA_igphyml_1, 
                         color=(MRCA_tribal_search_0.75 <= MRCA_igphyml_1))) +
  geom_point(size=3) + scale_color_discrete(name="TRIBAL <= IgPhyML") +
  facet_wrap(~cells, scales="free_y", nrow=2) + ylim(c(0, 0.0016))
  xlab("TRIBAL MRCA") + ylab("IgPhyML MRCA") 
  igphyml_mrca_comp
  
  igphyml_rf_comp <- ggplot(tree.wide.rf, aes(x=RF_tribal_score_0.75, y=RF_igphyml_1, 
                                                color=(RF_tribal_score_0.75 <= RF_igphyml_1))) +
    geom_point(size=2) + scale_color_discrete(name="TRIBAL <= IgPhyML") +
    facet_wrap(~cells, scales="free_y", nrow=2)  + ylim(c(0,25)) +
  
  xlab("TRIBAL RF") + ylab("IgPhyML RF") 
  
 ggplot(tree.wide.rf, aes(x=log(MRCA_tribal_score_0.75), y=log(MRCA_IgPhyML_1), 
                                              color=(MRCA_tribal_score_0.75 <= MRCA_IgPhyML_1))) +
  geom_point(size=3) + scale_color_discrete(name="TRIBAL <= IgPhyML") +
  facet_wrap(~cells, scales="free_y") +
  xlab("TRIBAL MRCA") + ylab("IgPhyML MRCA") 

igphyml_rf_comp <- ggplot(tree.wide.rf, aes(x=RF_tribal_search_0.95, y=RF_IgPhyML_1, 
                                              color=(RF_tribal_search_0.95 <= RF_IgPhyML_1))) +
  geom_point(size=3) + scale_color_discrete(name="TRIBAL <= IgPhyML") +
  facet_wrap(~cells, scales="free_y") +
  xlab("TRIBAL MRCA") + ylab("IgPhyML MRCA") 
igphyml_rf_comp
+ get_theme(base_size)+ #+xlim(c(0, 0.005)) +
theme(legend.position=c(.775,.86))



plot_save(file.path(plot_path, "igphyml_mrca_comp.pdf"), igphyml_mrca_comp)
dnapars_mrca_comp <-ggplot(tree.wide.rf, aes(x=MRCA_tribal_search_0.95, y=MRCA_dnapars_1, 
                         color=(MRCA_tribal_search_0.95 <= MRCA_dnapars_1))) +
  geom_point(size=2) + scale_color_discrete(name="TRIBAL <= dnapars") +
  xlab("TRIBAL MRCA") + ylab("dnapars MRCA") + my_theme +ylim(c(0, 0.0018))

dnaml_mrca_comp <- ggplot(tree.wide.rf, aes(x=MRCA_tribal_search_0.95, y=MRCA_dnaml_1, 
                         color=(MRCA_tribal_search_0.95 <= MRCA_dnaml_1))) +
  geom_point(size=2) + scale_color_discrete(name="TRIBAL <= dnaml") +
  xlab("TRIBAL MRCA") + ylab("dnaml MRCA") + my_theme +ylim(c(0, 0.0018))


   # scale_x_log10() + scale_y_log10() +
  #annotation_logticks() #+ ylim(c(log(1e-7), log(1e-2))) + xlim(c(log(1e-7), log(1e-2)))
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
