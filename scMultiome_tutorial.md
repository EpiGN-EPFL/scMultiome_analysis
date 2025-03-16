# Introduction

With the advancement of single-cell technologies, we can now profile
multiple modalities in the **same** cell. There are many examples of
such approaches:

-   [Multiplexed spatial mapping of chromatin features, transcriptome
    and proteins in
    tissues](https://www.nature.com/articles/s41592-024-02576-0)
-   [Single-cell multiplex chromatin and RNA interactions in ageing
    human brain](https://www.nature.com/articles/s41586-024-07239-w)
-   [Linking genome structures to functions by simultaneous single-cell
    Hi-C and
    RNA-seq](https://www.science.org/doi/10.1126/science.adg3797)
-   …

One of the most popular methods is the [Epi Multiome ATAC + Gene
Expression](https://www.10xgenomics.com/products/epi-multiome) from 10x
Genomics. You can find the workflow
[here](https://www.10xgenomics.com/support/epi-multiome).

In this tutorial, we will use single-cell RNA-ATAC multiomic data from
18-day-old brain organoids, as described in the study [Inferring and
perturbing cell fate regulomes in human brain
organoids](https://www.nature.com/articles/s41586-022-05279-8#Sec5). We
will analyze the scMultiome data in `R` using
[`Seurat`](https://satijalab.org/seurat/) and
[`Signac`](https://stuartlab.org/signac/)

**Optional**: You can download the raw sequencing data (FASTQ files)
from
[E-MTAB-12002](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12002).
For this tutorial, we re-generated transcript count and peak
accessibility matrices using
[`cellranger-arc count`](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/single-library-analysis)
(v2.0.2) with refdata-cellranger-arc-GRCh38-2020-A-2.0.0.

**Note**: we typically run cellranger pipelines on HPC clusters.

# Part 1 - Setup

## Step 0 - Load Required Libraries

Begin by loading the necessary libraries:

    library(EnsDb.Hsapiens.v86) # Annotation database 
    library(biovizBase)
    library(Seurat)
    library(Signac)
    library(Matrix)
    library(dplyr)
    library(tidyr)
    library(knitr) # This is to make rendered Markdown files look pretty

**Note**: Both `Seurat` and `Signac` are under active development. To
ensure compatibility, verify your package versions:

    packageVersion("Seurat")

    ## [1] '5.2.1'

## Step 1 - Understand the Experiment Setup

Before we dive into the actual analysis, it is always good to examine
the dataset. You can check the detailed sample information
[here](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12002/sdrf),
or from the downloaded files:

    sample_info <- read.csv("data/E-MTAB-12002/E-MTAB-12002.sdrf.txt", sep = '\t')

To summarize the sample modalities and genetic modifications:

    # Summarize sample information by genetic modification and modality
    sample_info %>%
      mutate(
        Modality = case_when(
          grepl("ATAC", Assay.Name) ~ "ATAC",
          grepl("RNA", Assay.Name) ~ "RNA",
          TRUE ~ "Other"
        )
      ) %>%
      group_by(Characteristics.genetic.modification.) %>%
      summarize(
        ATAC_samples = paste(unique(Assay.Name[Modality == "ATAC"]), collapse = ", "),
        RNA_samples = paste(unique(Assay.Name[Modality == "RNA"]), collapse = ", ")
      ) %>%
      kable()

<table>
<colgroup>
<col style="width: 37%" />
<col style="width: 32%" />
<col style="width: 30%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Characteristics.genetic.modification.</th>
<th style="text-align: left;">ATAC_samples</th>
<th style="text-align: left;">RNA_samples</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">CRISPR Cas9 mediated GLI3 knockout (29 and
41 bp deletion)</td>
<td style="text-align: left;">MULTIOME_ATAC_GLI3_B4_S1,
MULTIOME_ATAC_GLI3_B4_S2</td>
<td style="text-align: left;">MULTIOME_RNA_GLI3_B4_S1,
MULTIOME_RNA_GLI3_B4_S2</td>
</tr>
<tr class="even">
<td style="text-align: left;">CRISPR Cas9 mediated GLI3 knockout (4 and
8 bp insertion)</td>
<td style="text-align: left;">MULTIOME_ATAC_GLI3_D3_S1,
MULTIOME_ATAC_GLI3_D3_S2</td>
<td style="text-align: left;">MULTIOME_RNA_GLI3_D3_S1,
MULTIOME_RNA_GLI3_D3_S2</td>
</tr>
<tr class="odd">
<td style="text-align: left;">none</td>
<td style="text-align: left;">MULTIOME_ATAC_GLI3_A4_S1,
MULTIOME_ATAC_GLI3_A4_S2</td>
<td style="text-align: left;">MULTIOME_RNA_GLI3_A4_S1,
MULTIOME_RNA_GLI3_A4_S2</td>
</tr>
</tbody>
</table>

In this dataset, `S1` and `S2` denote the same biological sample
sequenced across two different lanes. There are three paired RNA+ATAC
datasets: two CRISPR knockouts and one control. Each pair corresponds to
a specific sample ID:

-   A4: Control (WT)
-   B4: CRISPR KO with a 29 and 41 bp deletion in GLI3
-   D3: CRISPR KO with a 4 and 8 bp insertion in GLI3

**PS**: [GLI3](https://en.wikipedia.org/wiki/GLI3) is zinc finger
protein which is involved in Sonic hedgehog (Shh) signaling.

## Step 2 - Understand the output

To streamline the analysis, multiple runs (i.e., results from
`cellranger-arc count`) were aggregated using
[`cellranger-arc aggr`](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/running-pipelines/aggregating-multiple-gem-wells-aggr).
The aggregation was guided by a `libraries_aggr.csv` file:

    libraries <- read.csv('data/libraries_aggr.csv')
    libraries %>%
      kable()

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 31%" />
<col style="width: 32%" />
<col style="width: 31%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">library_id</th>
<th style="text-align: left;">atac_fragments</th>
<th style="text-align: left;">per_barcode_metrics</th>
<th style="text-align: left;">gex_molecule_info</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">count_A4</td>
<td
style="text-align: left;">/scratch/hhu/scMultiome/arc-runs_A4/count_A4/outs/atac_fragments.tsv.gz</td>
<td
style="text-align: left;">/scratch/hhu/scMultiome/arc-runs_A4/count_A4/outs/per_barcode_metrics.csv</td>
<td
style="text-align: left;">/scratch/hhu/scMultiome/arc-runs_A4/count_A4/outs/gex_molecule_info.h5</td>
</tr>
<tr class="even">
<td style="text-align: left;">count_B4</td>
<td
style="text-align: left;">/scratch/hhu/scMultiome/arc-runs_B4/count_B4/outs/atac_fragments.tsv.gz</td>
<td
style="text-align: left;">/scratch/hhu/scMultiome/arc-runs_B4/count_B4/outs/per_barcode_metrics.csv</td>
<td
style="text-align: left;">/scratch/hhu/scMultiome/arc-runs_B4/count_B4/outs/gex_molecule_info.h5</td>
</tr>
<tr class="odd">
<td style="text-align: left;">count_D3</td>
<td
style="text-align: left;">/scratch/hhu/scMultiome/arc-runs_D3/count_D3/outs/atac_fragments.tsv.gz</td>
<td
style="text-align: left;">/scratch/hhu/scMultiome/arc-runs_D3/count_D3/outs/per_barcode_metrics.csv</td>
<td
style="text-align: left;">/scratch/hhu/scMultiome/arc-runs_D3/count_D3/outs/gex_molecule_info.h5</td>
</tr>
</tbody>
</table>

**Key points**:

-   For RNA data, `aggr` combines count matrices (not an integration
    step).
-   For ATAC data, `aggr` merges all fragments and performs new peak
    calling, ensuring a unified peak set (check
    [here](https://kb.10xgenomics.com/hc/en-us/articles/6057890578829-Does-cellranger-atac-aggr-redo-peak-calling)).

**Note**: Alternatively, peak calling for each sample can be performed
using tools like
[`MACS3`](https://macs3-project.github.io/MACS/docs/callpeak.html#) or
the
[`Signac::CallPeaks`](https://stuartlab.org/signac/articles/peak_calling).
Subsequently, the union/intersection of peaks across samples can be
selected for downstream analysis.

# Part 2 - Load the Data

## Step 1 - Load the Count Matrix

Similar to the standard pipeline for single-cell RNA-seq (scRNA-seq)
data, we can use the output folder
[`filtered_feature_bc_matrix`](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/single-library-analysis).

    count_aggr <- Read10X("data/filtered_feature_bc_matrix/")

    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.

Since this dataset includes two modalities (RNA & ATAC), the output will
contain two matrices. We can check their dimensions as follows:

    dim(count_aggr$`Gene Expression`)

    ## [1] 36601 25327

    dim(count_aggr$Peaks)

    ## [1] 301573  25327

Before proceeding, we confirm whether the barcodes (cell identities) in
both matrices match:

    all(colnames(count_aggr$`Gene Expression`) == colnames(count_aggr$Peaks))

    ## [1] TRUE

If TRUE, the barcodes are perfectly aligned across both modalities. If
FALSE, further investigation is needed to identify discrepancies.

### Exploring the RNA Assay

The RNA assay generates a **GENE × CELL** matrix, which already provides
useful information:

    rownames(count_aggr$`Gene Expression`) %>% head()

    ## [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3" 
    ## [6] "AL627309.2"

    seurat <- CreateSeuratObject(counts = count_aggr$`Gene Expression`,
                                     assay = "RNA")
    seurat

    ## An object of class Seurat 
    ## 36601 features across 25327 samples within 1 assay 
    ## Active assay: RNA (36601 features, 0 variable features)
    ##  1 layer present: counts

### Exploring the ATAC Assay

Unlike RNA data, chromatin accessibility data is stored as a **PEAK ×
CELL** matrix:

    rownames(count_aggr$Peaks) %>% head()

    ## [1] "chr1:9795-10680"    "chr1:180701-181102" "chr1:181191-181694"
    ## [4] "chr1:267570-268463" "chr1:585750-586644" "chr1:629483-630394"

Each feature (row) represents a genomic region (peak) in the format:
`chromosome:start-end`. And we will notice that the matrix is still very
big:

    dim(count_aggr$Peaks)

    ## [1] 301573  25327

To **not kill this R session because of memory requirements**, we can do
some filtering first:

    min_cells <- round(0.05 * ncol(count_aggr$Peaks))  # Keep peaks present in at least X% of cells
    filtered_peaks <- count_aggr$Peaks[rowSums(count_aggr$Peaks > 0) >= min_cells, ]
    min_reads <- 100  # Keep peaks with at least Y reads
    filtered_peaks <- filtered_peaks[rowSums(filtered_peaks) >= min_reads, ]
    dim(filtered_peaks)  # Check new dimensions

    ## [1] 55566 25327

To better interpret the peaks, we can annotate them with the **nearest**
genes. Many genomic databases and tools provide R interfaces, making R a
powerful environment for genomic data analysis. Here we follow the
[`Seurat` WNN
tutorial](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis#wnn-analysis-of-10x-multiome-rna-atac):

    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    seqlevelsStyle(annotations) <- 'UCSC'
    genome(annotations) <- "hg38"

    ################ NOTE #################
    ## StringToGRanges keep crushing (?) ###
    ## So we do it mannually ##
    ########################################

    # grange.counts <- StringToGRanges(rownames(filtered_peaks), sep = c(":", "-"))

    # rownames(filtered_peaks) contains chromosome coordinates in 'chr:start-end' format
    peak_names <- rownames(filtered_peaks)  # Example: c("chr1:1-10", "chr2:12-3121")

    # Split based on ":" first (to separate 'chr' from 'start-end')
    split_peaks <- strsplit(peak_names, ":")

    # Extract chromosome names
    chr <- sapply(split_peaks, `[`, 1)

    # Further split the second part (start-end) by "-"
    ranges <- sapply(split_peaks, `[`, 2)
    ranges_split <- strsplit(ranges, "-")

    # Extract start and end positions as integers
    start <- as.integer(sapply(ranges_split, `[`, 1))
    end <- as.integer(sapply(ranges_split, `[`, 2))

    # Create GRanges object
    grange.counts <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))

    # we'll only use peaks in standard chromosomes
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    filtered_peaks <- filtered_peaks[as.vector(grange.use), ]

We can check the filtered matrix dimension:

    dim(filtered_peaks)

    ## [1] 55548 25327

**Important**: There are other ways of how we could reduce the number of
`peaks`. For example, as mentioned above, if you are using
`cellranger-arc count` outputs directly, you can manually curate a list
of peaks that are present in multiple samples (check
[here](https://github.com/stuart-lab/signac/discussions/461)). For large
datasets, consider leveraging HPC clusters with higher memory
capacities.

    frag.file <- "data/atac_fragments.tsv.gz"
    # Finally create the ATAC assay
    seurat[['ATAC']] <-  CreateChromatinAssay(
       counts = filtered_peaks,
       ranges = grange.counts[grange.use],
       genome = 'hg38',
       fragments = frag.file,
       annotation = annotations
     )

    ## Computing hash

Now we have the seurat object with both RNA and ATAC assays:

    seurat

    ## An object of class Seurat 
    ## 92149 features across 25327 samples within 2 assays 
    ## Active assay: RNA (36601 features, 0 variable features)
    ##  1 layer present: counts
    ##  1 other assay present: ATAC

**Optional**: Save the object.

    saveRDS(
      object = seurat,
      file = "out/raw_seurat.Rds"
    )
