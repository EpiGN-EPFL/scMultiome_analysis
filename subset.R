library(Seurat)
library(Signac)
library(Matrix)

set.seed(42)
n_cells_keep <- 8000

#### load the count matrix ####
count_aggr <- Read10X("data//filtered_feature_bc_matrix/")

# subset both matrices
cells_keep <- sample(colnames(count_aggr$`Gene Expression`), n_cells_keep)
count_genes_small <- count_aggr$`Gene Expression`[,cells_keep]
count_peaks_small <- count_aggr$`Peak`[, cells_keep]

# filter the peaks
min_cells <- round(0.05 * ncol(count_peaks_small))  # Keep peaks present in at least X% of cells
count_peaks_small <- count_peaks_small[rowSums(count_peaks_small) >= min_cells, ]
min_reads <- 100  # Keep peaks with at least Y reads
count_peaks_small <- count_peaks_small[rowSums(count_peaks_small) >= min_reads, ]
dim(count_peaks_small)  # Check new dimensions

#### save the subsetted matrices ####

# 1) Stack features (rows) but keep same cells (columns)
stopifnot(identical(colnames(count_genes_small), colnames(count_peaks_small)))

mat_all <- rbind(count_genes_small, count_peaks_small)

# 2) Make 10x v3 features.tsv.gz (3 columns)
features <- rbind(
  data.frame(
    feature_id   = rownames(count_genes_small),
    feature_name = rownames(count_genes_small),
    feature_type = "Gene Expression",
    stringsAsFactors = FALSE
  ),
  data.frame(
    feature_id   = rownames(count_peaks_small),
    feature_name = rownames(count_peaks_small),  # for peaks, name often same as id (e.g. chr1:...-...)
    feature_type = "Peaks",
    stringsAsFactors = FALSE
  )
)

barcodes <- data.frame(barcode = colnames(mat_all), stringsAsFactors = FALSE)

# 3) Write out like 10x "filtered_feature_bc_matrix/"
outdir <- "data_small/filtered_feature_bc_matrix-small"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

writeMM(mat_all, file.path(outdir, "matrix.mtx"))

write.table(features, file.path(outdir, "features.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(barcodes, file.path(outdir, "barcodes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 4) gzip them to mimic 10x exactly
system(paste("gzip -f", file.path(outdir, "matrix.mtx")))
system(paste("gzip -f", file.path(outdir, "features.tsv")))
system(paste("gzip -f", file.path(outdir, "barcodes.tsv")))

#### fragments ####

writeLines(cells_keep, "cells_keep.txt")

# peaks as BED (chr start end)
peaks_gr <- as(rownames(count_peaks_small), "GRanges")
peaks_bed <- data.frame(
  chr = as.character(GenomeInfoDb::seqnames(peaks_gr)),
  start = GenomicRanges::start(peaks_gr) - 1,  # BED is 0-based start
  end = GenomicRanges::end(peaks_gr)
)
write.table(peaks_bed, "peaks_keep.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
