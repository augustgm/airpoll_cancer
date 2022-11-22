library(data.table)
library(ggplot2)
library(vcfR)
library(cowplot)
library(dplyr)
library(reshape2)
library(tidyverse)
library(deconstructSigs)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg19)
source("genome_utils.R")


comparisons_df = data.frame(
  comparison=names(vcf_files),
  category=factor(c(rep("PBS",5),rep("Pollution",5)),levels=c("PBS","Pollution"))
)


load("mutect_filtered.RData")

mutect_filtered_df = lapply(names(mutect_filtered),function(f_name){
  x = mutect_filtered[[f_name]]
  x$comparison = f_name
  x
})
mutect_filtered_df = bind_rows(mutect_filtered_df)

mutation_counts = mutect_filtered_df %>%
  group_by(mutation_id) %>%
  summarise(count=n(),
            unique_samples = length(unique(comparison)))

unique_mutations = mutation_counts %>%
  filter(unique_samples==1) %>%
  pull(mutation_id)

non_unique_mutations = mutation_counts %>%
  filter(unique_samples>1) %>%
  pull(mutation_id)

mutect_filtered_df %>%
  mutate(AF=as.numeric(AF)) %>%
  ggplot(aes(x=AF)) +
  geom_density() +
  facet_wrap(~comparison,nrow=2)



  

sigs.inputs = lapply(names(mutect_filtered),function(x){
  print(x)
  cdf = mutect_filtered[[x]]
  
  mutect_filtered_anno_DCS_df = cdf %>%
    filter(mutation_id%in%unique_mutations) %>%
    mutate(var_count=as.numeric(AF)*as.numeric(DP)) %>%
    filter(var_count>=3) %>%
    mutate(Sample=x,
           chr=CHROM,
           pos=POS,
           ref=REF,
           alt=ALT) %>%
    mutate(pos=as.numeric(as.character(pos))) %>%
    select(Sample,chr,pos,ref,alt)
  
  sigs.input <- mut.to.sigs.input(mut.ref = mutect_filtered_anno_DCS_df, 
                                  sample.id = "Sample", bsg=BSgenome.Mmusculus.UCSC.mm10,
                                  chr = "chr", 
                                  pos = "pos", 
                                  ref = "ref", 
                                  alt = "alt")
  
  sigs.input %>%
    melt() %>%
    mutate(nt_change=substr(variable,3,5)) %>%
    mutate(comparison=x)
    
})
names(sigs.inputs) = names(mutect_filtered)
sigs.inputs_df = bind_rows(sigs.inputs)

sig_plots = lapply(names(mutect_filtered),function(x){
  print(x)
  plotTriNucleotideContext(sigs.inputs[[x]],title_text =x)
})
plot_grid(plot_grid(plotlist = sig_plots[1:5],ncol=1),plot_grid(plotlist = sig_plots[6:10],ncol=1),ncol=2)

color_bar = c("#28bdee","#000000","#e62a27","#cbcacb","#a2c966","#ecc8c5")
sigs.inputs_df %>%
  mutate(variable=substr(variable,3,5)) %>%
  group_by(comparison,variable) %>%
  summarise(mutation_count=sum(value)) %>%
  left_join(.,comparisons_df,by="comparison") %>%
  ggplot(aes(comparison,y=mutation_count,fill=variable)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=color_bar) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~category,scales="free_x")


all_mutations_plot = sigs.inputs_df %>%
  mutate(variable=substr(variable,3,5)) %>%
  group_by(comparison) %>%
  summarise(mutation_count=sum(value)) %>%
  left_join(.,comparisons_df,by="comparison") %>%
  ggplot(aes(x=category,y=mutation_count)) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(aes(color=comparison)) +
  theme_cowplot() +
  theme(
    legend.position = "none"
  ) +
  ggtitle("All Mutations")

nt_change_counts_plot = sigs.inputs_df %>%
  mutate(variable=substr(variable,3,5)) %>%
  group_by(comparison,variable) %>%
  summarise(mutation_count=sum(value)) %>%
  left_join(.,comparisons_df,by="comparison") %>%
  ggplot(aes(x=category,y=mutation_count)) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(aes(color=comparison)) +
  facet_wrap(~variable,nrow=1) +
  theme_cowplot() 

plot_grid(all_mutations_plot,nt_change_counts_plot,rel_widths=c(1,4))

sigs.inputs_df %>%
  mutate(variable=substr(variable,3,5)) %>%
  group_by(comparison) %>%
  summarise(mutation_count=sum(value)) %>%
  left_join(.,comparisons_df,by="comparison") %>%
  group_by(category) %>%
  summarise(mutation_count=list(mutation_count)) %>%
  spread(category,mutation_count) %>%
  mutate(p_value = t.test(unlist(PBS), unlist(Pollution))$p.value,
         t_value = t.test(unlist(PBS), unlist(Pollution))$statistic)



sigs.inputs_df %>%
  mutate(variable=substr(variable,3,5)) %>%
  group_by(comparison,variable) %>%
  summarise(mutation_count=sum(value)) %>%
  left_join(.,comparisons_df,by="comparison") %>%
  group_by(variable,category) %>%
  summarise(mutation_count=list(mutation_count)) %>%
  spread(category,mutation_count) %>%
  group_by(variable)%>% 
  mutate(p_value = t.test(unlist(PBS), unlist(Pollution))$p.value,
         t_value = t.test(unlist(PBS), unlist(Pollution))$statistic)
  


tri_context_df = getMouseHumanGenomeContext()

tri_context_df %>%
  ggplot(aes(x=freq.mm10,y=freq.hg19)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label=tri_context_name)) +
  theme_cowplot() +
  geom_abline()

sigs.inputs_adj_df = sigs.inputs_df %>%
  mutate(tri_context_name=paste0(substr(variable,1,1),substr(variable,3,3),substr(variable,7,7))) %>%
  left_join(.,tri_context_df,by="tri_context_name") %>%
  mutate(adj_value=value*freq.hg19/freq.mm10) %>%
  group_by(comparison) %>%
  mutate(norm_adj_value=prop.table(adj_value))

sigs.inputs_adj_mm = sigs.inputs_adj_df %>%
  select(variable,comparison,norm_adj_value) %>%
  spread(variable,norm_adj_value,fill=0)
sigs.inputs_adj_m = sigs.inputs_adj_mm[,2:ncol(sigs.inputs_adj_mm)]
rownames(sigs.inputs_adj_m) = sigs.inputs_adj_mm$comparison


load("signatures.genome.cosmic.v3.may2019.rda")

sigs_to_consider = c("SBS1","SBS4","SBS5","SBS2","SBS13","SBS40","SBS92","SBS17a","SBS17b","SBS18")
which_sig_res = lapply(names(mutect_filtered),function(x){
  print(x)
  sig_res = whichSignatures(tumor.ref = sigs.inputs_adj_m, 
                            signatures.ref = signatures.genome.cosmic.v3.may2019, 
                            sample.id = x,
                            signature.cutoff=0.06,
                            associated=sigs_to_consider)
  df = sig_res$weights
  df$comparison = x
  df
})

which_sig_res = bind_rows(which_sig_res)



which_sig_res_unknown = lapply(names(mutect_filtered),function(x){
  print(x)
  sig_res = whichSignatures(tumor.ref = sigs.inputs_adj_m, 
                            signatures.ref = signatures.genome.cosmic.v3.may2019, 
                            sample.id = x,
                            signature.cutoff=0.06,
                            associated=sigs_to_consider)
  df = data.frame(unknown=sig_res$unknown,
                  comparison=x)
})

which_sig_res_unknown = bind_rows(which_sig_res_unknown)



plot_sigs = which_sig_res %>%
  melt(id.vars=c("comparison")) %>%
  filter(variable%in%sigs_to_consider) %>%
  left_join(.,comparisons_df,by=c("comparison")) %>%
  ggplot(aes(x=comparison,y=variable,fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red3") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~category,scales="free_x")


plot_unknown = which_sig_res_unknown %>%
  left_join(.,comparisons_df,by=c("comparison")) %>%
  ggplot(aes(x=comparison,y=unknown)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~category,scales="free_x") +
  ggtitle("Proportion Unassigned")

plot_grid(plot_sigs,plot_unknown,rel_widths=c(2,1))

which_sig_res %>%
  melt(id.vars=c("comparison")) %>%
  filter(variable%in%sigs_to_consider) %>%
  left_join(.,comparisons_df,by=c("comparison")) %>%
  group_by(category) %>%
  summarise(value_list=list(value)) %>%
  spread(category,value_list) %>%
  mutate(p_val=t.test(unlist(Pollution),unlist(PBS))$p.value)
