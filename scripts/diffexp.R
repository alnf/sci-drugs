library(Rsubread)
library(tidyverse)
library(yaml)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db) 

config.drugs <- unlist(read_yaml("../configs/config.yaml"))
metadata.drugs <- read_tsv("../metadata/drugs_metadata.tsv")

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

dds <- DESeqDataSetFromMatrix(counts, column_to_rownames(metadata.drugs, "sample"),
                              design = ~ timepoint + genotype + timepoint:genotype)

dds$timepoint <- relevel(dds$timepoint, "D0")
dds$genotype  <- relevel(dds$genotype, "499ctrl")

dds <- DESeq(dds)

res_D0 <- results(dds, contrast=c("genotype","499","499ctrl"))
res_D0 <- as.data.frame(res_D0)

res_D0_sig <- res_D0 %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) > 0)


res_D0_sig$ens_gene <- sub("\\..*$", "", rownames(res_D0_sig))

ann <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = unique(res_D0_sig$ens_gene),
  keytype  = "ENSEMBL",
  columns  = c("SYMBOL", "GENENAME")
)

ann <- ann |> distinct(ENSEMBL, .keep_all = TRUE)
colnames(ann) <- c("ens_gene", "symbol", "description")

res_annot <- res_D0_sig |>
  mutate(ens_gene = factor(ens_gene, levels = ens_gene)) |>
  left_join(ann, by = "ens_gene") |>
  arrange(ens_gene) |>
  relocate(ens_gene, symbol, description)

res_annot <- res_annot %>% filter(!is.na(symbol))

write.table(res_annot, "../results/res_annot_D0.tsv", sep="\t", row.names = F)

m <- counts(dds, normalized = TRUE)
m <- m[rownames(res_D0_sig), ]
m <- as.data.frame(m) %>%  dplyr::select(contains("D0"))

m <- merge(m, res_D0_sig[,c("ens_gene"),drop=F], by="row.names")
rownames(m) <- m$ens_gene
m <- m[,c(2:7)]


m <- m[res_annot$ens_gene,]
write.table(m, "../results/res_annot_D0_counts.tsv", sep="\t", row.names = T)