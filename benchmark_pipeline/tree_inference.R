
reps <-c(1:7)
cells <- 35
size <- 25

other_methods <- c("dnaml", "dnapars", "igphyml")
# other_methods <- c("igphyml")
source("/scratch/projects/tribal/init.R")
iso_res <- data.frame()
cols <- c("rep", "clonotype", "method", "tree_id", "N_taxa", "RF", "MRCA", "MRCA_ISO")
for(c in cell_values){
  for(m in other_methods){
    fname <- sprintf("/scratch/projects/tribal/benchmark_pipeline/sim_data/replicates/cells%d/size25/%s.tribal.tsv", c, m)
    print(fname) 
    if(file.exists(fname)){
      df <- read.table(fname, header=F, col.names = cols)
      df$cells <- c
      iso_res <- bind_rows(iso_res, df)
    }

     }

}


plot_path <- "/scratch/projects/tribal/benchmark_pipeline/figures"
#methods <- c("tribal_score", "tribal_score_gt", "tribal_search",  "tribal_search_gt", "tribal_refine", "tribal_refine_gt", "validaggreg")
methods <- c("tribal_score", "tribal_search", "tribal_refine")

vertx_theme <- get_theme_vertx(base_size)


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

tree.results <- bind_rows(tree.results, iso_res) %>% mutate(lambda= factor(ifelse(is.na(alpha),1, alpha))) %>%
mutate(method =as.character(method))

# dnaml1 <- filter(tree.results, method=="dnaml") %>% select(-RF, -MRCA_ISO,  -N_taxa) 
# dnaml2 <- filter(tree.results, method=="dnaml_tribal") %>% select(-MRCA,  -N_taxa) %>%mutate(method="dnaml", size=25)
# dnaml <- inner_join(dnaml2, dnaml1)
# 
# igphy1 <- filter(tree.results, method=="igphyml") %>% select(-RF, -MRCA_ISO,  -N_taxa) 
# igphy2 <- filter(tree.results, method=="igphyml_tribal") %>% select(-MRCA,  -N_taxa) %>%mutate(method="dnaml", size=25)
# igphy<- inner_join(igphy1, igphy2)
# 
# tree.results <- filter(tree.results, !(str_detect(method, "dnaml"))) %>% bind_rows(dnaml)

tree.results.long <- tree.results %>% pivot_longer( cols=c("RF", "MRCA_ISO", "MRCA")) %>%
  #filter(method != "dnapars") %>% 
  filter(value < 99) %>%
  #mutate(method = ifelse(method=="dnapars_tribal", "dnapars", method)) %>%
  filter(!str_detect(method, "gt"))

tree.results.long %>% group_by(cells, method,lambda, name) %>%
  summarize(mean=mean(value), sd=sd(value)) %>% View()

ggplot(tree.results.long, aes(x=method,fill=lambda, y=value)) + geom_boxplot() + 
  facet_grid(name~cells, scales="free_y") + vtext




ggplot(tree.results, aes(x=method, y=log(MRCA), fill=lambda)) + 
  geom_boxplot() + facet_wrap(~cells) + vtext

ggplot(tree.results, aes(x=method, y=RF, fill=lambda)) + 
  geom_boxplot() + facet_wrap(~cells) + vtext

ggplot(tree.results, aes(x=method, y=MRCA_ISO, fill=lambda)) + 
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
#   mutate(metric=ifelse(name=="MRCA_ISO", "isotype\nMRCA", ifelse(name=="MRCA", "sequence\nMRCA", "RF")))
tree.results.long <-mutate(tree.results.long ,metric=ifelse(name=="MRCA_ISO", "isotype\nMRCA", ifelse(name=="MRCA", "sequence\nMRCA", "RF")))
metric_box_plot <- ggplot(tree.results.long %>% filter(cells==35), aes(x=method, y=value, fill=factor(lambda))) + geom_boxplot() +
  facet_grid(metric~cells, scales="free_y")  + xlab("method") + ylab("") +
  theme(legend.position="top") + vtext + scale_fill_discrete(name="lambda")


metric_box_plot
plot_save(file.path(plot_path,"metric_box_plot.pdf"), metric_box_plot)
tree.result.long %>% group_by(method, name) %>% summarize(mean= mean(value, na.rm=T), sd= sd(value, na.rm=T))

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

ggplot(tree.results %>% filter(str_detect(method, "tribal")), 
       aes(x=method, fill=factor(alpha), y=MRCA_ISO)) + 
  geom_boxplot() + 
  # facet_wrap(~rep, labeller="label_both") +
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
tree.wide.rf <- tree.results %>%
  filter(RF < 99) %>%
  select( clonotype, cells, method, lambda, RF, rep, MRCA_ISO, MRCA) %>% 
  pivot_wider(id_cols =c("clonotype","rep", "cells"), names_from=c("method", "lambda"), values_from=c("RF", "MRCA", "MRCA_ISO"))


# ggplot(tree.wide.rf, aes(y=RF_dnapars_1-RF_tribal_refine_0.75 , x=MRCA_ISO_tribal_refine_0.75)) +
#   geom_point() + 
#   facet_wrap(~rep, labeller="label_both") + geom_smooth(method="lm")
# 
# ggplot(tree.wide.rf, aes(x=RF_tribal_refine_0.75, y=RF_dnapars_1, color=(RF_tribal_refine_0.75 <= RF_dnapars_1))) +
#   geom_point() + scale_color_discrete(name="tribal better") #+
#   # facet_wrap(~rep, labeller="label_both")

ggplot(tree.wide.rf, aes(x=RF_tribal_search_0.75, y=RF_IgPhyML_1, color=(RF_tribal_search_0.75 <=RF_IgPhyML_1))) +
  geom_point() + scale_color_discrete(name="tribal better") + facet_wrap(~cells, scales="free_y") +
  ylim(c(0,25))

ggplot(tree.wide.rf, aes(x=MRCA_tribal_search_0.95, y=MRCA_dnapars_1, 
                         color=(MRCA_tribal_search_0.95 <= MRCA_dnapars_1))) +
   geom_point() + scale_color_discrete(name="TRIBAL better") + 
  facet_wrap(~cells, scales="free_y")


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

igphyml_mrca_comp <- ggplot(tree.wide.rf, aes(x=MRCA_tribal_search_0.75, y=MRCA_IgPhyML_1, 
                         color=(MRCA_tribal_search_0.75 <= MRCA_IgPhyML_1))) +
  geom_point(size=3) + scale_color_discrete(name="TRIBAL <= IgPhyML") +
  facet_wrap(~cells, scales="free_y") + ylim(c(0, 0.0016))
  xlab("TRIBAL MRCA") + ylab("IgPhyML MRCA") 

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
