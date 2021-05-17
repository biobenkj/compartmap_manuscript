## Convert the bedpe(ish) output from the sci-HiC ML3 library for GM12878 and K562
## so that cooler can take them and convert to .cool files

## Distribution of data in terms of cell number
## GM12878    K562
## 326     593

## Read in annotation files to parse each cell type
sci_hic_metadata <- read.delim("valid_cells.list.labeled.human", header = FALSE,
                               stringsAsFactors = FALSE)

colnames(sci_hic_metadata) <- c("file_name", "cell_type")

# Split into respective cell types
sci_hic_metadata.gm12878 <- sci_hic_metadata[sci_hic_metadata$cell_type == "GM12878",]
sci_hic_metadata.k562 <- sci_hic_metadata[sci_hic_metadata$cell_type == "K562",]

# Save the metadata
write.table(sci_hic_metadata.gm12878, file = "sci_hic_gm12878_meta.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')
write.table(sci_hic_metadata.k562, file = "sci_hic_k562_meta.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

# Read in and process matrices
# We need to tile the genome at the same resolution as we did with bin_schic.py
# In this case we did 100kb
library(SummarizedExperiment)
library(Homo.sapiens)

hg19.100kb.tile <- tileGenome(seqlengths(Homo.sapiens)[standardChromosomes(Homo.sapiens)],
                              tilewidth = 100000, cut.last.tile.in.chrom = TRUE)

# Cast this to a vector for speedy bin reconstruction
hg19.100kb.tile.vec <- as.character(granges(hg19.100kb.tile))

# GM12878
sci_hic_gm12878 <- mclapply(sci_hic_metadata.gm12878$file_name, function(m) {
  #read in with no header
  #format is bin1, bin2, count, norm_count, chr1, chr2
  m.mat <- as.matrix(read.delim(m, header = FALSE, stringsAsFactors = FALSE))
  #reset the chrom names
  m.mat[,5:6] <- gsub("human_", "", m.mat[,5:6])
  #now we need to iterate through each line and rebuild bins
  m.mat.recon <- as.matrix(apply(m.mat, 1, function(b) {
    #chr1 from the tile vector
    b <- trimws(b)
    chr1 <- hg19.100kb.tile.vec
    #chr2 from the tile vector
    chr2 <- hg19.100kb.tile.vec
    #get the bin index as granges coordinates
    coord1 <- chr1[as.numeric(b[1])]
    coord2 <- chr2[as.numeric(b[2])]
    #get the count
    count <- as.numeric(b[3])
    #parse out the bedpe format that cooler expects
    #string that is split should look like: "chr21:9900001-10000000"
    coord1.spl <- strsplit(coord1, ":")
    bedpe.chr1 <- coord1.spl[[1]][1]
    coord1.spl <- strsplit(coord1.spl[[1]][2], "-")
    bedpe.start1 <- coord1.spl[[1]][1]
    bedpe.end1 <- coord1.spl[[1]][2]
    #second coordinate
    coord2.spl <- strsplit(coord2, ":")
    bedpe.chr2 <- coord2.spl[[1]][1]
    coord2.spl <- strsplit(coord2.spl[[1]][2], "-")
    bedpe.start2 <- coord2.spl[[1]][1]
    bedpe.end2 <- coord2.spl[[1]][2]
    bedpe.out <- c(bedpe.chr1, bedpe.start1, bedpe.end1,
                   bedpe.chr2, bedpe.start2, bedpe.end2, count)
    return(bedpe.out)
  }))
  return(t(m.mat.recon))
}, mc.cores = 16)

names(sci_hic_gm12878) <- gsub(".matrix", "", sci_hic_metadata.gm12878$file_name)

# Write out the new bedpe files with proper coordinates...
lapply(1:length(sci_hic_gm12878), function(bedpe) {
  write.table(sci_hic_gm12878[[bedpe]], file = paste0(names(sci_hic_gm12878)[bedpe], ".bedpe.bed"),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
})

# K562
# remove a bad cell(s)
bad_k562_cells <- c("human_2158_CTGAGGCA-ACCTCTTG_100000.matrix")
sci_hic_metadata.k562 <- sci_hic_metadata.k562[!sci_hic_metadata.k562$file_name %in% bad_k562_cells,]
sci_hic_k562 <- lapply(sci_hic_metadata.k562$file_name, function(m) {
  #read in with no header
  #format is bin1, bin2, count, norm_count, chr1, chr2
  message("Working on file ", m)
  m.mat <- as.matrix(read.delim(m, header = FALSE, stringsAsFactors = FALSE))
  #reset the chrom names
  m.mat[,5:6] <- gsub("human_", "", m.mat[,5:6])
  #now we need to iterate through each line and rebuild bins
  m.mat.recon <- as.matrix(apply(m.mat, 1, function(b) {
    #chr1 from the tile vector
    b <- trimws(b)
    chr1 <- hg19.100kb.tile.vec
    #chr2 from the tile vector
    chr2 <- hg19.100kb.tile.vec
    #get the bin index as granges coordinates
    coord1 <- chr1[as.numeric(b[1])]
    coord2 <- chr2[as.numeric(b[2])]
    #get the count
    count <- as.numeric(b[3])
    #parse out the bedpe format that cooler expects
    #string that is split should look like: "chr21:9900001-10000000"
    coord1.spl <- strsplit(coord1, ":")
    bedpe.chr1 <- coord1.spl[[1]][1]
    coord1.spl <- strsplit(coord1.spl[[1]][2], "-")
    bedpe.start1 <- coord1.spl[[1]][1]
    bedpe.end1 <- coord1.spl[[1]][2]
    #second coordinate
    coord2.spl <- strsplit(coord2, ":")
    bedpe.chr2 <- coord2.spl[[1]][1]
    coord2.spl <- strsplit(coord2.spl[[1]][2], "-")
    bedpe.start2 <- coord2.spl[[1]][1]
    bedpe.end2 <- coord2.spl[[1]][2]
    bedpe.out <- c(bedpe.chr1, bedpe.start1, bedpe.end1,
                   bedpe.chr2, bedpe.start2, bedpe.end2, count)
    return(bedpe.out)
  }))
  return(t(m.mat.recon))
})

names(sci_hic_k562) <- gsub(".matrix", "", sci_hic_metadata.k562$file_name)

# Write out the new bedpe files with proper coordinates...
lapply(1:length(sci_hic_k562), function(bedpe) {
  write.table(sci_hic_k562[[bedpe]], file = paste0(names(sci_hic_k562)[bedpe], ".bedpe.bed"),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
})
