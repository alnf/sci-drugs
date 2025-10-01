library(Rsubread)
library(tidyverse)
library(yaml)
library(DESeq2)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(matrixStats)
library(dplyr)
library(tibble)
source("utils.R")

config.drugs <- unlist(read_yaml("../configs/config.yaml"))
metadata.drugs <- read_tsv("../metadata/drugs_metadata.tsv")
rownames(metadata.drugs) <- metadata.drugs$sample

bam.files <- file.path(config.drugs["data.star"],
                       paste0(metadata.drugs$sample, "/", metadata.drugs$sample, ".Aligned.sortedByCoord.out.bam"))

fc <- featureCounts(
  files = bam.files,
  annot.ext = config.drugs["star.sjdbGTFfile"],
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  isPairedEnd = TRUE,
  strandSpecific = 2,
  nthreads = 28
)

counts <- as.data.frame(fc$counts)
colnames(counts) <- meta$sample   # keep names clean & in metadata order
stopifnot(all(colnames(counts) == meta$sample))


### Differential expression between genotypes for each time point separately
dds <- DESeqDataSetFromMatrix(counts, column_to_rownames(metadata.drugs, "sample"),
                              design = ~ 0 + group)

# Filter very low-expressed genes
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds)
resultsNames(dds)
rownames(dds) <- sub("\\..*$", "", rownames(dds))

# Annotate genes
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

bm <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters    = "ensembl_gene_id",
  values     = unique(rownames(dds)),
  mart       = mart
) |>
  distinct(ensembl_gene_id, .keep_all = TRUE)


# Normalised counts
lcounts <- counts(dds, normalized = TRUE)
write.table(m, "../results/drugs_lcounts.tsv", sep="\t", row.names = T)

# Heatmap of all counts
mat <- t(scale(t(as.matrix(lcounts))))
mat[!is.finite(mat)] <- 0

top_n <- 10000
vars  <- matrixStats::rowVars(mat)
keep  <- order(vars, decreasing = TRUE)[seq_len(min(length(vars), top_n))]
mat <- mat[keep, , drop = FALSE]


pal_cat <- function(n, pal = "Set2") {
  x <- brewer.pal(max(3, min(n, 8)), pal)
  if (n > length(x)) x <- colorRampPalette(x)(n)
  if (n < length(x)) x <- colorRampPalette(x)(n)
  x
}
col_group    <- setNames(pal_cat(n_distinct(metadata.drugs$group), "Set2"),
                         sort(unique(metadata.drugs$group)))
col_time     <- setNames(pal_cat(n_distinct(metadata.drugs$timepoint), "Set3"),
                         sort(unique(metadata.drugs$timepoint)))
col_genotype <- setNames(pal_cat(n_distinct(metadata.drugs$genotype), "Dark2"),
                         sort(unique(metadata.drugs$genotype)))

top_ha <- HeatmapAnnotation(
  df  = metadata.drugs[, c("group","timepoint","genotype")],
  col = list(group = col_group, timepoint = col_time, genotype = col_genotype),
  annotation_legend_param = list(
    group    = list(title = "Group"),
    timepoint= list(title = "Timepoint"),
    genotype = list(title = "Genotype")
  )
)

col_fun <- colorRamp2(c(min(mat), 0, max(mat)), c("#2166AC", "white", "#B2182B"))

ht <- Heatmap(
  mat,
  name = "Expression",
  col = col_fun,
  top_annotation = top_ha,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  column_title = "Samples",
  row_title = "Genes",
  # handy extras:
  # column_split = meta2$timepoint,    # split columns by timepoint
  # row_km = 4                          # k-means rows into 4 clusters
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")


# Differential expression
res <- results(dds, contrast=c("genotype","499","499ctrl"))
res <- as.data.frame(res)

res$ens_gene <- sub("\\..*$", "", rownames(res))

anno <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = unique(rownames(dds)),
  keytype  = "ENSEMBL",
  columns  = c("SYMBOL", "GENENAME")
)

anno <- anno |> distinct(ENSEMBL, .keep_all = TRUE)
colnames(anno) <- c("ens_gene", "symbol", "description")

res_anno <- res |>
  mutate(ens_gene = factor(ens_gene, levels = ens_gene)) |>
  left_join(anno, by = "ens_gene") |>
  arrange(ens_gene) |>
  relocate(ens_gene, symbol, description)

res_anno <- res_anno %>% filter(baseMean != 0)


write.table(res_anno, "../results/res_annot_D0.tsv", sep="\t", row.names = F)

m <- counts(dds, normalized = TRUE)
m <- m[rownames(res), ]
m <- as.data.frame(m) %>%  dplyr::select(contains("D0"))

m <- merge(m, res_D0_sig[,c("ens_gene"),drop=F], by="row.names")
rownames(m) <- m$ens_gene
m <- m[,c(2:7)]


m <- m[res_annot$ens_gene,]
write.table(m, "../results/res_annot_D0_counts.tsv", sep="\t", row.names = T)