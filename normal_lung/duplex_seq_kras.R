library(data.table)
library(ggplot2)
library(janitor)
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

duplex_anno_df = fread("/camp/project/proj-tracerx-lung/tctProjects/duplex/data/tx_peace_jdg/process_sr_ltx/results/PEACE_BDRE_duplex_anno_df.txt")
  
duplex_anno_plot_df = duplex_anno_df %>%  
  filter(SYMBOL=="KRAS") %>%
  filter((!is.na(isDriver)&grepl("Strong",isDriver))|!is.na(cosmic70)) %>%
  mutate(isDriver=ifelse(is.na(isDriver),"Other",isDriver)) %>%
  mutate(exon_id=ifelse(POS>25390000,"Exon2","Exon3")) %>%
  mutate(mutation_aa=ifelse(is.na(protein_coordinate),paste(exon_id,ALT),mutation_aa)) %>%
  mutate(Protein_position=ifelse(is.na(Protein_position),ALT,Protein_position)) %>%
  arrange(protein_coordinate,status) %>%
  group_by(protein_coordinate) %>%
  mutate(row_index=row_number())

# Gene Model Plot
vaf_plot = duplex_anno_plot_df %>%
  mutate(mutation_aa=factor(mutation_aa,levels=unique(duplex_anno_plot_df$mutation_aa))) %>%
  mutate(protein_coordinate=factor(protein_coordinate)) %>%
  ggplot(aes(x=mutation_aa,y=log10(vaf))) +
  ggbeeswarm::geom_quasirandom(aes(color=isDriver,shape=status)) +
  scale_color_manual(values=c("red3","royalblue3","black")) +
  ylim(-4.5,-1.5) +
  facet_wrap(~sample_type,ncol=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

exon_df = data.frame(
  exon_names = c(2,3),
  start = c(1,36),
  end = c(35,96)
)
exon_df$bottom = 0.9
exon_df$top = 1.1

exon_plot = exon_df %>%
  ggplot() +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=bottom, ymax=top),fill="grey90") +
  theme_cowplot() +
  geom_point(data=duplex_anno_plot_df,aes(x=protein_coordinate,y=row_index+1,color=isDriver,shape=status)) +
  scale_color_manual(values=c("red3","royalblue3","black")) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )

plot_grid(exon_plot,vaf_plot,ncol=1,rel_heights = c(3,5))


# KRAS By Cancer Type
duplex_anno_plot_df %>%
  filter(sample_type=="Autopsy") %>%
  mutate(mutation_aa=factor(mutation_aa,levels=unique(duplex_anno_plot_df$mutation_aa))) %>%
  mutate(protein_coordinate=factor(protein_coordinate)) %>%
  mutate(smoking_status=ifelse(smoking_status=="","Unknown",smoking_status)) %>%
  ggplot(aes(x=cancer_type_group,y=log10(vaf))) +
  ggbeeswarm::geom_quasirandom(aes(color=BRAF_inhibitor,shape=smoking_status),size=2) +
  scale_color_manual(values=c("darkorchid1","forestgreen","black")) +
  ylim(-4.5,-1.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", color = "black", width = 0.7, lwd = 0.2)