# Introduction

With the advancement of single-cell technologies, we can now profile
multiple modalities in the **same** cell. There are many examples of
such approaches:

TODO: Add relevant papers.

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
(v2.0.2) and mapped to refdata-cellranger-arc-GRCh38-2020-A-2.0.0.
**Note**: we typically run the cellranger pipelines on HPC clusters.

# Part 1 - Setup

We start by loading the libraries:

    library(Seurat)
    library(Signac)
    library(Matrix)
    library(dplyr)
    library(tidyr)
    library(knitr) # This is to make rendered Markdown files look pretty

**Important**: Both `Seurat` and `Signac` are under active development.
You can check the package version like this:

    packageVersion("Seurat")

    ## [1] '5.2.1'

Before we dive into the actual analysis, it is always good to check what
are the samples we are dealing with. You can check the detailed sample
information
[here](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12002/sdrf),
or from the downloaded files:

    sample_info <- read.csv("data/E-MTAB-12002/E-MTAB-12002.sdrf.txt", sep = '\t')

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

S1 and S2 refers to the the same biological sample but were sequenced
across two different lanes. We can see that there are three paired
RNA+ATAC datasets: two CRISPR knockouts and one control. Each pair
corresponds to a specific sample ID:

-   B4: CRISPR KO with a 29 and 41 bp deletion in GLI3
-   D3: CRISPR KO with a 4 and 8 bp insertion in GLI3
-   A4: Control (WT)

**PS**: [GLI3](https://en.wikipedia.org/wiki/GLI3) is zinc finger
protein which is involved in Sonic hedgehog (Shh) signaling.
