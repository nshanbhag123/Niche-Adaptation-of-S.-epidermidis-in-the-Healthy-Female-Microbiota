### 1. Load accessory matrix
acc <- read.table("/Users/sjablonska/Desktop/staph_matrix-PRESENCE-ABSENCE.txt",
                  header=TRUE,
                  row.names=1,
                  sep="\t",
                  check.names=FALSE)

### 2. Keep numeric only
acc <- acc[, sapply(acc, is.numeric)]

### 3. Keep accessory-only genes (<114 genomes)
acc <- acc[rowSums(acc) < ncol(acc), ]

### 4. Transpose to samples × genes
acc_t <- t(acc)

### 5. Load metadata
meta <- read.table("/Users/sjablonska/Desktop/Sanotherone.txt",
                   header=TRUE, sep="\t")

### Rename columns properly
colnames(meta)[1] <- "Sample"
colnames(meta)[2] <- "Site"

### 6. Ensure exact matching
meta$sample <- trimws(meta$sample)
rownames(acc_t) <- trimws(rownames(acc_t))

### Align metadata order to acc_t order
meta <- meta[match(rownames(acc_t), meta$sample), ]

### 7. Compute Jaccard distance
library(vegan)
dist_acc <- vegdist(acc_t, method="jaccard")

### 8. Run PCoA
pcoa_acc <- cmdscale(dist_acc, eig=TRUE, k=2)

### Variance explained
pc1_var <- round(100 * pcoa_acc$eig[1] / sum(pcoa_acc$eig), 2)
pc2_var <- round(100 * pcoa_acc$eig[2] / sum(pcoa_acc$eig), 2)

### 9. Plot dataframe
pcoa_df <- data.frame(
  sample = rownames(pcoa_acc$points),
  PC1 = pcoa_acc$points[,1],
  PC2 = pcoa_acc$points[,2],
  site = meta$site
)

### 10. Plot
library(ggplot2)
ggplot(pcoa_df, aes(PC1, PC2, color=site)) +
  geom_point(size=4) +
  theme_bw() +
  labs(
    title="PCoA of Accessory Genome (Presence/Absence)",
    x=paste0("PC1 (", pc1_var, "%)"),
    y=paste0("PC2 (", pc2_var, "%)")
  )

### 11. PERMANOVA
adonis2(dist_acc ~ site, data=pcoa_df)


library(vegan)

# your accessory matrix
acc <- acc

# your metadata
site <- pcoa_df$site

# run enrichment gene-by-gene
pvals <- apply(acc, 1, function(gene) {
  df <- data.frame(gene = gene, site = site)
  fit <- chisq.test(table(df$gene, df$site))
  fit$p.value
})

# adjust for multiple testing
padj <- p.adjust(pvals, method = "fdr")

sig_genes <- names(padj[padj < 0.05])
sig_genes
