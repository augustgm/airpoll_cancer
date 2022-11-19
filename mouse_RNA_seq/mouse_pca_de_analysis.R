# load packages
library(DESeq2)
library(tidyverse)

# Load count table
data <- read.table("RNAseqEATA/RSEM_RawReadCounts.txt", header = TRUE, sep = "\t")
data <- data[, -2]
data <- data %>%
  column_to_rownames("gene_id") %>%
  as.matrix() %>%
  round(digits = 0)

# Define custom variables
projectTitle <- "exp1"

# set column data
coldata <- data.frame(group = c(rep("EA", times = 5), 
                                rep("EC", times = 5), 
                                rep("TA", times = 5), 
                                rep("TC", times = 5)),
                      row.names = colnames(data))

# create a coldata frame and initiate DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = data, 
                              colData = coldata, 
                              design = ~ group)

# filter genes with >1 cpm in 2 or more samples: 14582
keep <- rowSums(edgeR::cpm(data) > 1) >= 2
dds <- dds[keep, ]
table(keep)

# estimate size factors, dispersion and perform negative binomial Wald test
dds <- dds %>%
  estimateSizeFactors() %>%
  estimateDispersions() %>%
  nbinomWaldTest()

# variance stabilized for PCA
vsd <- vst(dds, blind = FALSE)

# get genes for PC1 and PC2
pca <- assay(vsd) %>%
  t() %>%
  prcomp()
pca_contribution <- pca$rotation %>%
  abs() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  select(gene_id, PC1, PC2)
pc1_genes <- pca_contribution %>%
  arrange(-PC1) %>%
  select(-PC2)
pc2_genes <- pca_contribution %>%
  arrange(-PC2) %>%
  select(-PC1)
save(pc1_genes, file = "pc1_genes.RData")
save(pc2_genes, file = "pc2_genes.RData")

# get raw PCA data and plot
pca_data_raw <- plotPCA(vsd, intgroup = c("group"), returnData = TRUE)
variance <- attr(pca_data_raw, "percentVar")
pca_data_raw %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(col = group), size = 3) +
  ggrepel::geom_text_repel(aes(label = name)) +
  labs(x = paste0("PC1: ", round(100*variance[1], 2), "% variance"),
       y = paste0("PC2: ", round(100*variance[2], 2),"% variance"),
       col = "Group") +
  theme_bw() +
  theme(aspect.ratio = 1)
ggsave(file = paste0(projectTitle, "_pca_plot.pdf"))

# do a pairwise comparison for contrasts of interest
contrasts <- list(c("EA", "EC"), 
                  c("TA", "TC"),
                  c("EA", "TA"),
                  c("EC", "TC")) # set contrasts
res <- list()
for(i in seq_along(contrasts)) {
  res[[i]] <- results(dds, alpha = 0.05, contrast = c("group", contrasts[[i]]))
}
names(res) <- map_chr(contrasts, function(x) {paste0(x, collapse = "_")})

# annotate genes
mart <- biomaRt::useEnsembl(biomart = "ensembl", 
                            dataset = "mmusculus_gene_ensembl",
                            host = "uswest.ensembl.org")
ids <- c("ensembl_gene_id", 'external_gene_name', 'chromosome_name', 
         "band", 'start_position', 'end_position')

# get annotation for the filtered genes
ensembl_gene_id_anno <- biomaRt::getBM(attributes = ids, 
                                       mart = mart, 
                                       values = rownames(data)[keep])

# extract results and annotate genes
dds_res <- map(res, function(x) {
  x %>%
    as.data.frame() %>% 
    tibble::rownames_to_column("ensembl_gene_id") %>%
    as_tibble() %>%
    arrange(padj, pvalue) %>%
    left_join(ensembl_gene_id_anno, by = "ensembl_gene_id")
})

# plot MA plots for both contrasts
dds_res_total <- bind_rows(dds_res, .id = "contrast")

dds_res_total %>%
  filter(!is.na(padj)) %>%
  mutate(signif = ifelse(padj < 0.05, "Padj < 0.05", "N.S")) %>%
  ggplot(aes(x = baseMean, y = log2FoldChange, col = signif, alpha = signif)) +
  geom_point(shape = 1) +
  scale_x_log10() +
  scale_color_manual(values = c("Padj < 0.05" = "red", "N.S." = "black")) +
  scale_alpha_manual(values = c(0.25, 1)) +
  facet_wrap(~contrast, scales = "free_y") +
  labs(x = "Mean abundance (TPM)", y = "log2FC", col = "") +
  theme_bw() +
  theme(aspect.ratio = 1)
ggsave(file = paste0(projectTitle, "_ma_plot.pdf"))

write.csv(x = dds_res$EA_EC, file = paste0(projectTitle, "_EA_EC.csv"), quote = FALSE, row.names = FALSE)
write.csv(x = dds_res$TA_TC, file = paste0(projectTitle, "_TA_TC.csv"), quote = FALSE, row.names = FALSE)
write.csv(x = dds_res$EA_TA, file = paste0(projectTitle, "_EA_TA.csv"), quote = FALSE, row.names = FALSE)
write.csv(x = dds_res$EC_TC, file = paste0(projectTitle, "_EC_TC.csv"), quote = FALSE, row.names = FALSE)


