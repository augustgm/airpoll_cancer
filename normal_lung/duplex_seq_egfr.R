library(data.table)
library(ggplot2)
library(janitor)
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

duplex_anno_df = fread("/camp/project/proj-tracerx-lung/tctProjects/duplex/data/tx_peace_jdg/process_sr_ltx/results/PEACE_BDRE_duplex_anno_df.txt")

duplex_anno_plot_df = duplex_anno_df %>%  
  filter(SYMBOL=="EGFR") %>%
  filter((!is.na(isDriver)&grepl("Strong",isDriver))|!is.na(cosmic70)) %>%
  mutate(isDriver=ifelse(is.na(isDriver),"Other",isDriver)) %>%
  arrange(protein_coordinate,status) %>%
  group_by(protein_coordinate) %>%
  mutate(row_index=row_number())

# Gene Model Plot
vaf_plot = duplex_anno_plot_df %>%
  mutate(mutation_aa=factor(mutation_aa,levels=unique(duplex_anno_plot_df$mutation_aa))) %>%
  mutate(protein_coordinate=factor(protein_coordinate)) %>%
  ggplot(aes(x=mutation_aa,y=log10(vaf))) +
  ggbeeswarm::geom_quasirandom(aes(color=isDriver,shape=status)) +
  scale_color_manual(values=c("salmon","red3","royalblue3","black")) +
  ylim(-4.5,-1.5) +
  facet_wrap(~sample_type,ncol=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

exon_df = data.frame(
  exon_names = c(18,19,20,21),
  start = c(688,729,762,824),
  end = c(728,761,823,875)
)
exon_df$bottom = 0.9
exon_df$top = 1.1

exon_plot = exon_df %>%
  ggplot() +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=bottom, ymax=top),fill="grey90")+
  theme_cowplot() +
  geom_point(data=duplex_anno_plot_df,aes(x=protein_coordinate,y=row_index+1,color=isDriver,shape=status)) +
  scale_color_manual(values=c("salmon","red3","royalblue3","black")) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )


plot_grid(exon_plot,vaf_plot,ncol=1,rel_heights = c(1,5))

