# immune loci
library(GenomicRanges)

immune <- read.csv("/Volumes/nfs_home/data/PanCan_genome_features_lustre/raw/immune_loci_mart_export.txt")
immune <- immune[!grepl(pattern = 'PATCH', immune$Chromosome.Name),]
immune <- GRanges(paste0("chr", immune$Chromosome.Name), IRanges(immune$Gene.Start..bp., immune$Gene.End..bp.), strand=immune$Strand, name=immune$Associated.Gene.Name)
names(immune) <- immune$name

# save to data/ dir
devtools::use_data(immune, overwrite=TRUE)


# genes
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exs <- exonsBy(txdb, "gene")
txs <- transcriptsBy(txdb, "gene")
gene_ids <- names(txs)
gene_name_map <- select(org.Hs.eg.db, gene_ids, c("SYMBOL"))
names(exs) <- gene_name_map$SYMBOL
names(txs) <- gene_name_map$SYMBOL
remove(txdb, gene_ids, gene_name_map)

# save to data/ dir
devtools::use_data(exs, overwrite=TRUE)
devtools::use_data(txs, overwrite=TRUE)
