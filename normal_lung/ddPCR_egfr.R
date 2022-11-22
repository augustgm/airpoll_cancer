library(data.table)
library(ggplot2)
library(janitor)
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)
library(patchwork)

saga_results_anno_ltx_peace = fread("/camp/project/proj-tracerx-lung/tctProjects/lime/projects/egfr_normal/data/SAGA_Report_20210702/saga_results_anno_ltxpea_df_20220404_TND.txt")

plot_freq_count = function(col_oi,saga_results_anno_ltx_include){
  color_bar = c("grey","#1b9e77","#d95f02","#4a1486","#7570b3","#7570b3","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")
  names(color_bar) = c("Absent","exon19dels","g719s","l858r","g719s;l858r","l858r;l861q","l858r;s768i","l861q","l861q;s768i","s768i")
  
  plot_df = saga_results_anno_ltx_include %>%
    data.frame() %>%
    melt(id.vars=c("patient_id","Normal_EGFRm_Summary_TND")) %>%
    filter(variable==col_oi) 
  
  
  count_plot = plot_df %>%
    group_by(value,Normal_EGFRm_Summary_TND) %>%
    tally() %>%
    ggplot(aes(x=value,y=n,fill=Normal_EGFRm_Summary_TND)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=color_bar) +
    theme_cowplot() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle(paste0("Exclusively Normal Mutations - ",col_oi))
  
  freq_plot = plot_df %>%
    group_by(value,Normal_EGFRm_Summary_TND) %>%
    tally() %>%
    mutate(freq=prop.table(n)) %>%
    ggplot(aes(x=value,y=freq,fill=Normal_EGFRm_Summary_TND)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=color_bar) +
    theme_cowplot()     +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
  
  freq_present_plot = plot_df %>%
    filter(Normal_EGFRm_Summary_TND!="Absent") %>%
    group_by(value,Normal_EGFRm_Summary_TND) %>%
    tally() %>%
    mutate(freq=prop.table(n)) %>%
    ggplot(aes(x=value,y=freq,fill=Normal_EGFRm_Summary_TND)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=color_bar) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
  
  plot_grid(count_plot,freq_plot,freq_present_plot,ncol=1,align="v")
  
}

## TX PEACE
saga_results_anno_ltx_include = saga_results_anno_ltx_peace 

color_bar = c("grey","#1b9e77","#d95f02","#4a1486","#7570b3","#7570b3","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")
names(color_bar) = c("Absent","exon19dels","g719s","l858r","g719s;l858r","l858r;l861q","l858r;s768i","l861q","l861q;s768i","s768i")

TX_PEACE_count_plot = saga_results_anno_ltx_include %>%
  data.frame() %>%
  melt(id.vars=c("patient_id","Normal_EGFRm_Summary_TND")) %>%
  filter(variable=="TX_PEACE") %>%
  group_by(value,Normal_EGFRm_Summary_TND) %>%
  tally() %>%
  mutate(Normal_EGFRm_Summary_TND=factor(Normal_EGFRm_Summary_TND,levels=names(color_bar))) %>%
  mutate(freq=prop.table(n)) %>%
  ggplot(aes(x=value,y=n,fill=Normal_EGFRm_Summary_TND)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=color_bar) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ylab("Sample Count")

TX_PEACE_freq_plot = saga_results_anno_ltx_include %>%
  data.frame() %>%
  melt(id.vars=c("patient_id","Normal_EGFRm_Summary_TND")) %>%
  filter(variable=="TX_PEACE") %>%
  group_by(value,Normal_EGFRm_Summary_TND) %>%
  tally() %>%
  mutate(Normal_EGFRm_Summary_TND=factor(Normal_EGFRm_Summary_TND,levels=names(color_bar))) %>%
  mutate(freq=prop.table(n)) %>%
  ggplot(aes(x=value,y=freq,fill=Normal_EGFRm_Summary_TND)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=color_bar) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


TX_PEACE_count_freq_plot = (TX_PEACE_count_plot / TX_PEACE_freq_plot) +
    plot_layout(guides = "collect")


## Carbon Review
saga_results_anno_ltx_include = saga_results_anno_ltx_peace %>%
  filter(TX_PEACE=="TRACERx")

color_bar = c("grey","#1b9e77","#d95f02","#4a1486","#7570b3","#7570b3","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")
names(color_bar) = c("Absent","exon19dels","g719s","l858r","g719s;l858r","l858r;l861q","l858r;s768i","l861q","l861q;s768i","s768i")

Carbon_count_plot = saga_results_anno_ltx_include %>%
    data.frame() %>%
    melt(id.vars=c("patient_id","Normal_EGFRm_Summary_TND")) %>%
    filter(variable=="parenchyma_carbon") %>%
    filter(!is.na(value)) %>%
    group_by(value,Normal_EGFRm_Summary_TND) %>%
    tally() %>%
    mutate(Normal_EGFRm_Summary_TND=factor(Normal_EGFRm_Summary_TND,levels=names(color_bar))) %>%
    mutate(freq=prop.table(n)) %>%
    ggplot(aes(x=value,y=n,fill=Normal_EGFRm_Summary_TND)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=color_bar) +
    theme_cowplot() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank()) +
    ylab("Sample Count")

  
  Carbon_freq_plot = saga_results_anno_ltx_include %>%
    data.frame() %>%
    melt(id.vars=c("patient_id","Normal_EGFRm_Summary_TND")) %>%
    filter(variable=="parenchyma_carbon") %>%
    filter(!is.na(value)) %>%
    group_by(value,Normal_EGFRm_Summary_TND) %>%
    tally() %>%
    mutate(Normal_EGFRm_Summary_TND=factor(Normal_EGFRm_Summary_TND,levels=names(color_bar))) %>%
    mutate(freq=prop.table(n)) %>%
    ggplot(aes(x=value,y=freq,fill=Normal_EGFRm_Summary_TND)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=color_bar) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab("parenchyma_carbon")

Carbon_count_freq_plot = (Carbon_count_plot / Carbon_freq_plot) +
  plot_layout(guides="collect")

  

 # VAFs
status_to_include = c("Norm_Only","Tum_Norm_Distinct")
multi_mut_exclusion_cases = c("LTX0162","LTX0934","LTX0946","LTX0566")

vaf_to_plot = saga_results_ltx_peace %>%
  mutate(VAF=(as.numeric(gsub("%","",VAF))/100)) %>%
  dplyr::rename(mutation_match_status=match_status) %>%
  filter(mutation_match_status%in%status_to_include) %>%
  filter(!patient_id%in%multi_mut_exclusion_cases) %>%
  left_join(.,saga_results_anno_ltx_peace,by=c("patient_id"="patient_id")) %>%
  filter(!is.na(passed_mutation),VAF>0) %>%
  group_by(TX_PEACE,external_patient_id,variable,VAF,RawVAF,sex,smoking_status,histology,age_group,tumour_EGFRm,mean_pm25_2007_2016,
           parenchyma_carbon,ln_carbon) %>%
  summarise(mutation_match_status=paste(unique(mutation_match_status[!is.na(mutation_match_status)]),collapse=";")) %>%
  mutate(tumour_EGFRm_short=ifelse(tumour_EGFRm%in%c("Absent"),tumour_EGFRm,"EGFRm")) %>%
  mutate(PM2.5_group=ifelse(mean_pm25_2007_2016>12.271,"High >12.271","Int")) %>%
  mutate(PM2.5_group=ifelse(mean_pm25_2007_2016<=10.579,"Low <10.579",PM2.5_group))


color_bar = c("#1b9e77","#916906","#7570b3","#e7298a","#e4ab02")
names(color_bar) = c("exon19dels","g719s","l858r","l861q","s768i")


TX_plot = vaf_to_plot %>%
  filter(TX_PEACE=="TRACERx") %>%
  filter(!is.na(parenchyma_carbon)) %>%
  ggplot(aes(x= parenchyma_carbon  ,y=log10(VAF))) +
  ggbeeswarm::geom_quasirandom(aes(shape= smoking_status ,color=PM2.5_group),size=2) +
  ylab("log10(VAF)") +
    theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +

  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", color = "black", width = 0.7, lwd = 0.2) +
  ylim(-5,-1) +
  stat_compare_means(method="t.test")
