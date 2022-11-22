library(ggpubr)
library(data.table)
library(ggplot2)
library(janitor)
library(dplyr)
library(cowplot)
library(tidyr)
library(patchwork)

duplex_anno_df = fread("PEACE_BDRE_duplex_anno_df.txt")

mut_count_df = duplex_anno_df %>%
  left_join(.,carbon_review_df,by=c("sample"="patient_id")) %>%
  mutate(smoking_status=ifelse(smoking_status=="","Unknown",smoking_status)) %>%
  group_by(cohort,sample,smoking_status,age,parenchyma_carbon,sex) %>%
  tally()

cos_mut_count_df = duplex_anno_df %>%
  left_join(.,carbon_review_df,by=c("sample"="patient_id")) %>%
  filter(!is.na(cosmic70)|mutation_aa%in%strong_driver) %>%
  mutate(smoking_status=ifelse(smoking_status=="","Unknown",smoking_status)) %>%
  group_by(cohort,sample,smoking_status,age,parenchyma_carbon,sex) %>%
  tally()

cosmic_mut_plot = cos_mut_count_df %>%
  filter(smoking_status=="Never Smoked") %>%
  filter(cohort=="PEACE") %>%
  ggplot(aes(x=age,y=n)) +
  geom_point(size=3) +
  stat_cor(method = "spearman")  +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~cohort+smoking_status,scales="free") +
  theme_cowplot() +
  ylab("Mutation Count") +
  xlab("Age") +
  ggtitle("Driver Mutations Identified Using Duplex-seq")