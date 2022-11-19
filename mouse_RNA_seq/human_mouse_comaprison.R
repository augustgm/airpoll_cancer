# Load libraries
library(biomaRt)
library(data.table)
library(dplyr)
library(glue)
library(ggplot2)
library(ggrepel)

# Get human and mouse orthologues 
human_ensembl <-  useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes  <-  c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene")
orth.mouse <- getBM(attributes,
                    filters="with_mmusculus_homolog",
                    values=TRUE,
                    mart = human_ensembl,
                    uniqueRows = TRUE)
                    
# Load DEG results table
# Human COPA study
carlsten_de_never_smoker <- read.table("NeverSmoker differential expression.csv", header=T, fill=T, sep =",")
# Mouse study
bjorn_ta_tc = read.table("exp1_TA_TC.csv", header=T, fill=T, sep=",")

# Compare DEG results between human and mouse
df = carlsten_de_never_smoker %>%
  rename(ensembl_gene_id = geneID) %>%
  right_join(orth.mouse) %>%
  dplyr::select(baseMean, log2FoldChange, pvalue, padj, geneSymbol, ensembl_gene_id, mmusculus_homolog_ensembl_gene) %>%
  rename(human_mean=baseMean, human_log2FoldChange=log2FoldChange, human_pvalue=pvalue, human_padj=padj, human_geneSymbol=geneSymbol) %>%
  right_join(bjorn_ta_tc %>%
               rename(mmusculus_homolog_ensembl_gene = ensembl_gene_id) %>%
               dplyr::select(mmusculus_homolog_ensembl_gene, baseMean, log2FoldChange, pvalue, padj, geneSymbol, mmusculus_homolog_ensembl_gene) %>%
               rename(mouse_mean=baseMean, mouse_log2FoldChange=log2FoldChange, mouse_pvalue=pvalue, mouse_padj=padj, mouse_geneSymbol=geneSymbol)) %>%
  filter(human_geneSymbol %ni% c(NA)) %>%
  mutate(geneSymbol = ifelse(human_geneSymbol != mouse_geneSymbol, 
                             glue("{human_geneSymbol} ({mouse_geneSymbol})"), human_geneSymbol)) %>%
  mutate(concordance = ifelse(human_log2FoldChange > 0 & mouse_log2FoldChange > 0, "concordant",
                              ifelse(human_log2FoldChange < 0 & mouse_log2FoldChange < 0, "concordant", "discordant")),
         label = ifelse(abs(`mouse_log2FoldChange`)>1.5, "label",
                        ifelse(abs(human_log2FoldChange)>1.5,"label", "no label")))
                        
# Plot human and mouse genes comparison
ggplot(df,
       aes(human_log2FoldChange, mouse_log2FoldChange)) +
  geom_rect(data = df[1,], aes(ymin = -Inf, ymax = 0, xmin = -Inf, xmax = 0), fill = "navy", alpha = 0.5) +
  geom_rect(data = df[1,], aes(ymin = 0, ymax = Inf, xmin = 0, xmax = Inf), fill = "firebrick", alpha = 0.5) +
  geom_point(aes(colour = concordance)) +
  xlab("log2FC (COPA human data)") + ylab("log2FC (Mouse data)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(data = df[df$label == "label",], aes(label = geneSymbol, colour = concordance, size = concordance), max.overlaps = 100, segment.size = 0.05, force = 10) +
  #geom_text_repel(data = df[df$concordance == "concordant",], aes(label = geneSymbol), colour = "black", size = 2.5, max.overlaps = Inf, segment.size = 0.1, force = 10) +
  scale_size_manual(values = c(2.5, 1.5)) +
  scale_colour_manual(name = "Human-mouse concordance", values = c("black", "#bababa")) +
  scale_fill_manual(name = "Human - mouse concordance", values = c("blue", "red"))+
  theme_classic() +
  theme(text = element_text(size=14, family = "Arial"), 
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none") +
  xlab("log2FC (COPA human data)") + ylab("log2FC (Mouse data)")


