pth <- "/scratch/projects/tribal/experimental_data"
table_fname <- "day_14_dandelion_table.tsv"
root_fname <- "day_14_root_sequences.csv"
# 
run_dir <- "day_14"
# 
run_path <- file.path(pth, run_dir)
dat_fname = sprintf("%s/%s_dandelion_table.tsv", run_path, run_dir)
root_fname = sprintf("%s/%s_root_sequences.csv", run_path, run_dir)
# 
# min_size = 5
# out_path = sprintf("/scratch/projects/tribal/real_data/%s/input", run_dir)
# 
#  outfile = sprintf("%s/clonotype_summary.csv",run_path)
# id_mapping = sprintf("%s/barcode_id_mapping.csv", run_path)
min_size <- 5
library(tidyverse)
###########################
dat_fname <- snakemake@input[['data_fname']]
root_fname <- snakemake@input[['root_fname']]
min_size <- as.numeric(snakemake@params[['min_size']])
out_path <- snakemake@params[['pth']]
outfile <- snakemake@output[["summary"]]
id_mapping <-  snakemake@output[["id_mapping"]]
##########################
cols <- c("barcode" , "heavy_constant_seq", 	"heavy_isotype", 	
          "light_constant_seq", "light_isotype", 	"heavy_v_seq", 
          "heavy_v_allele", "heavy_d_allele",	"heavy_j_allele",	
          "light_v_seq",	"light_v_allele",	"light_j_allele",	"clonotype",	
          "heavy_isotype_exp")
dat <- read.table(dat_fname, header=F, skip=1, col.names=cols) %>%
      select(barcode, clonotype, heavy_isotype_exp,  heavy_v_seq, light_v_seq)
print("read data")

good_clonos <- dat %>% group_by(clonotype) %>% count() %>% filter(n >= min_size) %>% pull(clonotype)
n_distinct(dat$clonotype)
dat %>% group_by(clonotype) %>% count()  %>% summary()


#dat.filt %>% group_by(clonotype) %>% count() %>% ungroup() %>% summarize(med= median(n), max=max(n))
root.cols <- c("clonotype", "light_v_seq", "heavy_v_seq")
root.dat <- read.csv(root_fname,header=F, skip=1, col.names= root.cols) %>%
  filter(clonotype %in% good_clonos) %>%
  mutate(barcode="naive", heavy_isotype_exp="Ighm")


dat.filt <- filter(dat, clonotype %in% good_clonos)
dat.filt %>% group_by(clonotype) %>% count() %>% summary()
dat.filt <- bind_rows(dat.filt, root.dat)
nrow(dat.filt)
n_distinct(dat.filt$clonotype)


#dat.filt <- dat.filt %>% mutate(heavy_len = str_length(heavy_v_seq), light_len = str_length(light_v_seq))

#dat.len <- dat.filt %>% group_by(clonotype) %>% summarize(dist_light_len = n_distinct(light_len), dist_heavy_len = n_distinct(heavy_len))

dat.filt <- dat.filt %>% group_by(clonotype) %>%  
  mutate(id = ifelse(barcode != "naive", paste("seq", row_number(), sep=""), barcode))

write.fasta <- function(ids, seqs, fname){
  # fileConn<-file(fname, open="w")
  sink(fname)
  for(i in 1:length(ids)){
    writeLines(sprintf(">%s\n%s", ids[i], seqs[i]) )

  }
  sink(file=NULL)
  # close(fileConn)
}


# out_path <- file.path(pth, "input")
if(!dir.exists(out_path)){
  dir.create(out_path)
}

write_all_fastas <- function(clono, dat.f, out_path){
  clono.dat <- filter(dat.f, clonotype==clono)
  clono_out_path <- file.path(out_path, clono)
  if(!dir.exists(clono_out_path)){
    dir.create(clono_out_path)
  }
  ids <- clono.dat$id
  light_fname <- file.path(clono_out_path, "light.fasta")
  
  write.fasta(ids, as.character(clono.dat$light_v_seq), light_fname)

  heavy_fname <- file.path(clono_out_path,"heavy.fasta")

  write.fasta(ids, as.character(clono.dat$heavy_v_seq), heavy_fname)
  iso_fname <- file.path(clono_out_path, "isotype.fasta")
  write.fasta(ids, as.character(clono.dat$heavy_isotype_exp), iso_fname)

  
}

lapply(good_clonos,write_all_fastas, dat.filt, out_path)

dat.filt %>% select(clonotype, barcode, id) %>% 
  write.csv(id_mapping, row.names=F, quote=F)
dat.filt %>% group_by(clonotype) %>% count() %>%
  write.csv(outfile, row.names=F, quote=F)

