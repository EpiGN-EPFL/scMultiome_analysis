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
    head(sample_info)

    ##             Source.Name Comment.ENA_SAMPLE. Comment.BioSD_SAMPLE.
    ## 1 MULTIOME_ATAC_GLI3_D3         ERS12569171        SAMEA110471145
    ## 2 MULTIOME_ATAC_GLI3_D3         ERS12569171        SAMEA110471145
    ## 3 MULTIOME_ATAC_GLI3_D3         ERS12569171        SAMEA110471145
    ## 4 MULTIOME_ATAC_GLI3_D3         ERS12569171        SAMEA110471145
    ## 5 MULTIOME_ATAC_GLI3_D3         ERS12569171        SAMEA110471145
    ## 6 MULTIOME_ATAC_GLI3_D3         ERS12569171        SAMEA110471145
    ##   Characteristics.organism. Characteristics.cell.line. Characteristics.sex.
    ## 1              Homo sapiens                      409B2               female
    ## 2              Homo sapiens                      409B2               female
    ## 3              Homo sapiens                      409B2               female
    ## 4              Homo sapiens                      409B2               female
    ## 5              Homo sapiens                      409B2               female
    ## 6              Homo sapiens                      409B2               female
    ##   Characteristics.age. Unit.time.unit. Characteristics.developmental.stage.
    ## 1                   36            year                                adult
    ## 2                   36            year                                adult
    ## 3                   36            year                                adult
    ## 4                   36            year                                adult
    ## 5                   36            year                                adult
    ## 6                   36            year                                adult
    ##   Characteristics.disease. Characteristics.organism.part.
    ## 1                   normal                         dermis
    ## 2                   normal                         dermis
    ## 3                   normal                         dermis
    ## 4                   normal                         dermis
    ## 5                   normal                         dermis
    ## 6                   normal                         dermis
    ##   Characteristics.progenitor.cell.type. Characteristics.genotype.
    ## 1         induced pluripotent stem cell    inducible Cas9 nickase
    ## 2         induced pluripotent stem cell    inducible Cas9 nickase
    ## 3         induced pluripotent stem cell    inducible Cas9 nickase
    ## 4         induced pluripotent stem cell    inducible Cas9 nickase
    ## 5         induced pluripotent stem cell    inducible Cas9 nickase
    ## 6         induced pluripotent stem cell    inducible Cas9 nickase
    ##                       Characteristics.genetic.modification.
    ## 1 CRISPR Cas9 mediated GLI3 knockout (4 and 8 bp insertion)
    ## 2 CRISPR Cas9 mediated GLI3 knockout (4 and 8 bp insertion)
    ## 3 CRISPR Cas9 mediated GLI3 knockout (4 and 8 bp insertion)
    ## 4 CRISPR Cas9 mediated GLI3 knockout (4 and 8 bp insertion)
    ## 5 CRISPR Cas9 mediated GLI3 knockout (4 and 8 bp insertion)
    ## 6 CRISPR Cas9 mediated GLI3 knockout (4 and 8 bp insertion)
    ##   Characteristics.growth.condition. Material.Type  Protocol.REF Protocol.REF.1
    ## 1         cerebral organoid culture       nucleus P-MTAB-124502  P-MTAB-124503
    ## 2         cerebral organoid culture       nucleus P-MTAB-124502  P-MTAB-124503
    ## 3         cerebral organoid culture       nucleus P-MTAB-124502  P-MTAB-124503
    ## 4         cerebral organoid culture       nucleus P-MTAB-124502  P-MTAB-124503
    ## 5         cerebral organoid culture       nucleus P-MTAB-124502  P-MTAB-124503
    ## 6         cerebral organoid culture       nucleus P-MTAB-124502  P-MTAB-124503
    ##   Protocol.REF.2 Protocol.REF.3 Protocol.REF.4          Extract.Name
    ## 1  P-MTAB-124504  P-MTAB-124505  P-MTAB-124506 MULTIOME_ATAC_GLI3_D3
    ## 2  P-MTAB-124504  P-MTAB-124505  P-MTAB-124506 MULTIOME_ATAC_GLI3_D3
    ## 3  P-MTAB-124504  P-MTAB-124505  P-MTAB-124506 MULTIOME_ATAC_GLI3_D3
    ## 4  P-MTAB-124504  P-MTAB-124505  P-MTAB-124506 MULTIOME_ATAC_GLI3_D3
    ## 5  P-MTAB-124504  P-MTAB-124505  P-MTAB-124506 MULTIOME_ATAC_GLI3_D3
    ## 6  P-MTAB-124504  P-MTAB-124505  P-MTAB-124506 MULTIOME_ATAC_GLI3_D3
    ##   Material.Type.1 Comment.LIBRARY_LAYOUT. Comment.LIBRARY_SELECTION.
    ## 1             DNA                  PAIRED                        PCR
    ## 2             DNA                  PAIRED                        PCR
    ## 3             DNA                  PAIRED                        PCR
    ## 4             DNA                  PAIRED                        PCR
    ## 5             DNA                  PAIRED                        PCR
    ## 6             DNA                  PAIRED                        PCR
    ##   Comment.LIBRARY_SOURCE. Comment.LIBRARY_STRATEGY.
    ## 1     GENOMIC SINGLE CELL                  ATAC-seq
    ## 2     GENOMIC SINGLE CELL                  ATAC-seq
    ## 3     GENOMIC SINGLE CELL                  ATAC-seq
    ## 4     GENOMIC SINGLE CELL                  ATAC-seq
    ## 5     GENOMIC SINGLE CELL                  ATAC-seq
    ## 6     GENOMIC SINGLE CELL                  ATAC-seq
    ##   Comment.single.cell.isolation. Comment.library.construction.
    ## 1                            10x                10x scATAC-seq
    ## 2                            10x                10x scATAC-seq
    ## 3                            10x                10x scATAC-seq
    ## 4                            10x                10x scATAC-seq
    ## 5                            10x                10x scATAC-seq
    ## 6                            10x                10x scATAC-seq
    ##   Comment.end.bias. Comment.input.molecule. Comment.primer. Comment.spike.in.
    ## 1              none             genomic DNA           other              none
    ## 2              none             genomic DNA           other              none
    ## 3              none             genomic DNA           other              none
    ## 4              none             genomic DNA           other              none
    ## 5              none             genomic DNA           other              none
    ## 6              none             genomic DNA           other              none
    ##   Comment.cdna.read. Comment.cdna.read.offset. Comment.cdna.read.size.
    ## 1       read1, read2                         0                      49
    ## 2       read1, read2                         0                      49
    ## 3       read1, read2                         0                      49
    ## 4       read1, read2                         0                      49
    ## 5       read1, read2                         0                      49
    ## 6       read1, read2                         0                      49
    ##   Comment.cell.barcode.read. Comment.cell.barcode.offset.
    ## 1                     index1                            0
    ## 2                     index1                            0
    ## 3                     index1                            0
    ## 4                     index1                            0
    ## 5                     index1                            0
    ## 6                     index1                            0
    ##   Comment.cell.barcode.size. Comment.umi.barcode.read.
    ## 1                         16                          
    ## 2                         16                          
    ## 3                         16                          
    ## 4                         16                          
    ## 5                         16                          
    ## 6                         16                          
    ##   Comment.umi.barcode.offset. Comment.umi.barcode.size.
    ## 1                          NA                        NA
    ## 2                          NA                        NA
    ## 3                          NA                        NA
    ## 4                          NA                        NA
    ## 5                          NA                        NA
    ## 6                          NA                        NA
    ##   Comment.sample.barcode.read. Comment.sample.barcode.offset.
    ## 1                       index2                              0
    ## 2                       index2                              0
    ## 3                       index2                              0
    ## 4                       index2                              0
    ## 5                       index2                              0
    ## 6                       index2                              0
    ##   Comment.sample.barcode.size. Protocol.REF.5               Performer
    ## 1                            8  P-MTAB-124507 DBSSE Genomics Facility
    ## 2                            8  P-MTAB-124507 DBSSE Genomics Facility
    ## 3                            8  P-MTAB-124507 DBSSE Genomics Facility
    ## 4                            8  P-MTAB-124507 DBSSE Genomics Facility
    ## 5                            8  P-MTAB-124507 DBSSE Genomics Facility
    ## 6                            8  P-MTAB-124507 DBSSE Genomics Facility
    ##                 Assay.Name  Technology.Type Comment.ENA_EXPERIMENT.
    ## 1 MULTIOME_ATAC_GLI3_D3_S1 sequencing assay              ERX9577217
    ## 2 MULTIOME_ATAC_GLI3_D3_S1 sequencing assay              ERX9577217
    ## 3 MULTIOME_ATAC_GLI3_D3_S1 sequencing assay              ERX9577217
    ## 4 MULTIOME_ATAC_GLI3_D3_S1 sequencing assay              ERX9577217
    ## 5 MULTIOME_ATAC_GLI3_D3_S2 sequencing assay              ERX9577217
    ## 6 MULTIOME_ATAC_GLI3_D3_S2 sequencing assay              ERX9577217
    ##                                       Scan.Name
    ## 1 MULTIOME_ATAC_GLI3_D3_S1_L001_I1_001.fastq.gz
    ## 2 MULTIOME_ATAC_GLI3_D3_S1_L001_R1_001.fastq.gz
    ## 3 MULTIOME_ATAC_GLI3_D3_S1_L001_R2_001.fastq.gz
    ## 4 MULTIOME_ATAC_GLI3_D3_S1_L001_R3_001.fastq.gz
    ## 5 MULTIOME_ATAC_GLI3_D3_S2_L002_I1_001.fastq.gz
    ## 6 MULTIOME_ATAC_GLI3_D3_S2_L002_R1_001.fastq.gz
    ##                    Comment.SUBMITTED_FILE_NAME. Comment.ENA_RUN.
    ## 1 MULTIOME_ATAC_GLI3_D3_S1_L001_I1_001.fastq.gz      ERR10036853
    ## 2 MULTIOME_ATAC_GLI3_D3_S1_L001_R1_001.fastq.gz      ERR10036853
    ## 3 MULTIOME_ATAC_GLI3_D3_S1_L001_R2_001.fastq.gz      ERR10036853
    ## 4 MULTIOME_ATAC_GLI3_D3_S1_L001_R3_001.fastq.gz      ERR10036853
    ## 5 MULTIOME_ATAC_GLI3_D3_S2_L002_I1_001.fastq.gz      ERR10036854
    ## 6 MULTIOME_ATAC_GLI3_D3_S2_L002_R1_001.fastq.gz      ERR10036854
    ##                                                                                  Comment.FASTQ_URI.
    ## 1 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10036853/MULTIOME_ATAC_GLI3_D3_S1_L001_I1_001.fastq.gz
    ## 2 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10036853/MULTIOME_ATAC_GLI3_D3_S1_L001_R1_001.fastq.gz
    ## 3 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10036853/MULTIOME_ATAC_GLI3_D3_S1_L001_R2_001.fastq.gz
    ## 4 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10036853/MULTIOME_ATAC_GLI3_D3_S1_L001_R3_001.fastq.gz
    ## 5 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10036854/MULTIOME_ATAC_GLI3_D3_S2_L002_I1_001.fastq.gz
    ## 6 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10036854/MULTIOME_ATAC_GLI3_D3_S2_L002_R1_001.fastq.gz
    ##   Comment.READ_TYPE. Comment.READ_INDEX.
    ## 1     sample_barcode              index2
    ## 2             paired               read1
    ## 3       cell_barcode              index1
    ## 4             paired               read2
    ## 5     sample_barcode              index2
    ## 6             paired               read1
    ##                          Factor.Value.genetic.modification.
    ## 1 CRISPR Cas9 mediated GLI3 knockout (4 and 8 bp insertion)
    ## 2 CRISPR Cas9 mediated GLI3 knockout (4 and 8 bp insertion)
    ## 3 CRISPR Cas9 mediated GLI3 knockout (4 and 8 bp insertion)
    ## 4 CRISPR Cas9 mediated GLI3 knockout (4 and 8 bp insertion)
    ## 5 CRISPR Cas9 mediated GLI3 knockout (4 and 8 bp insertion)
    ## 6 CRISPR Cas9 mediated GLI3 knockout (4 and 8 bp insertion)
    ##   Factor.Value.library.construction.
    ## 1                     10x scATAC-seq
    ## 2                     10x scATAC-seq
    ## 3                     10x scATAC-seq
    ## 4                     10x scATAC-seq
    ## 5                     10x scATAC-seq
    ## 6                     10x scATAC-seq

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
        ATAC_samples = paste(Assay.Name[Modality == "ATAC"], collapse = ", "),
        RNA_samples = paste(Assay.Name[Modality == "RNA"], collapse = ", ")
      )

    ## # A tibble: 3 × 3
    ##   Characteristics.genetic.modification.                 ATAC_samples RNA_samples
    ##   <chr>                                                 <chr>        <chr>      
    ## 1 CRISPR Cas9 mediated GLI3 knockout (29 and 41 bp del… MULTIOME_AT… MULTIOME_R…
    ## 2 CRISPR Cas9 mediated GLI3 knockout (4 and 8 bp inser… MULTIOME_AT… MULTIOME_R…
    ## 3 none                                                  MULTIOME_AT… MULTIOME_R…

From this, we can see that there are three paired RNA+ATAC datasets: two
CRISPR knockouts and one control. Each pair corresponds to a specific
sample ID: - B4: CRISPR KO with a 29 and 41 bp deletion in GLI3 - D3:
CRISPR KO with a 4 and 8 bp insertion in GLI3 - A4: Control (WT)

**PS**: [GLI3](https://en.wikipedia.org/wiki/GLI3) is zinc finger
protein which is involved in Sonic hedgehog (Shh) signaling.
