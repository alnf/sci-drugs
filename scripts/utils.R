getDEGs <- function(dds, contrast, anno, lFC = 1, sval.filter = TRUE, filename, filename_short) {
  res <- results(dds, contrast=contrast)
  res <- lfcShrink(dds, contrast=contrast, type="ashr", lfcThreshold = log2(1.2), alpha=0.05)
  resTable <- data.frame(res)
  resTable$baseMean <- round(resTable$baseMean, 4)
  
  baseMeanPerLvl <- sapply(levels(dds$Group),
                           function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$Group == lvl, drop=F]))   
  baseMeanPerLvl <- baseMeanPerLvl[, contrast[2:3]]
  baseMeanPerLvl <- round(baseMeanPerLvl, 4)
  
  colnames(baseMeanPerLvl) <- paste("mean", colnames(baseMeanPerLvl), sep="_")
  resTable <- cbind(baseMeanPerLvl, resTable)
  resTable$ens_gene <- rownames(resTable)
  
  resTable <- merge(t2g, resTable, by="ens_gene")
  resTable <- resTable[!duplicated(resTable$ens_gene),]
  resTable$svalue[which(resTable$svalue<0)] <- 0
  resTable <- resTable %>% arrange(desc(abs(log2FoldChange)), padj)
  
  write.table(resTable, filename, sep="\t", row.names = F)
  
  resTable <- resTable %>% filter(padj <= 0.05)
  if (sval.filter) {
    resTable <- resTable %>% filter(svalue <= 0.005)
  }
  resTable <- resTable %>% filter(!str_detect(description, 'ribosomal RNA'))
  resTable <- resTable %>% filter(symbol!="" | msymbol!="") 
  resTable$symbol[which(resTable$symbol=="")] <- resTable$msymbol[which(resTable$symbol=="")]
  
  mean <- rowMeans(resTable[,6:7])
  if (length(which(mean<20))>0) {
    resTable <- resTable[-which(mean<20),]
  }
  resTable <- resTable %>% filter(abs(log2FoldChange) > lFC)
  print(nrow(resTable))
  write.table(resTable, filename_short, sep="\t", row.names = F)
}