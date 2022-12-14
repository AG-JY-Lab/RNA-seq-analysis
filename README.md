# RNA-seq end to end analysis in R

![](./images/rnaseq.png)

### This repository maintains an end to end downstream workflow (post FASTQ alignment) of bulk RNA sequencing analysis from gene read counts to functional enrichment analysis.

<br>

> ### This workflow has been developed and tested on an Ubuntu 22.04 and 20.04 LTS machine so the individual R packages might not fully work on a Windows PC or on MacOS but users are more than welcome to try running it and take note of specific dependencies which might not be cross platform compatible (for possible alternatives) and make improvements to support other operating systems.

<br>

## Before you begin:

> Ensure that you have [renv](https://rstudio.github.io/renv/index.html) installed as this will the main environment management tool used in this workflow. The packages used in the analysis can be easily obtained with `renv::restore()` and this will download them by referencing the `renv.lock` file.

> **For Windows:** Some of the tools / packages in this analysis workflow might need to be built from source and the respective compilers for C/C++ & Fortran are nicely bundled within [`Rtools`](https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html) for Windows. [`Rtools`](https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html) can be installed by running the Rtools42 installer.

> Although this is not a standalone package for bulk RNA-seq analysis but more of a pipeline that integrates multiple tools together, there are some helper functions that were written as R source code within the `/src` directory. In order to load them for use in the analysis, run `source(here::here("src", "RNA_seq_helper_functions.R"))` at the start of each Rmd notebook to ensure the necessary helper functions are loaded too.

> The read count tab data are expected as outputs from [STAR](https://github.com/alexdobin/STAR) with its `--quantMode GeneCounts` flag turned on. As such, the format for the data would be a few `_outputs_ReadsPerGene.out.tab` files according to the number of RNA-seq samples. These should be placed in the `/data` folder along with a tab delimited metadata text file which the user has to provide for DESeq2 to build it's deseq_object.
> <br>

    ??????RNA-seq-analysis
        ?????????data
            ?????????Wildtype_1_outputs_ReadsPerGene.out.tab
            ?????????Wildtype_2_outputs_ReadsPerGene.out.tab
            ?????????Wildtype_3_outputs_ReadsPerGene.out.tab
            ?????????Wildtype_4_outputs_ReadsPerGene.out.tab
            ?????????Gene_xyz_KO_1_outputs_ReadsPerGene.out.tab
            ?????????Gene_xyz_KO_2_outputs_ReadsPerGene.out.tab
            ?????????Gene_xyz_KO_3_outputs_ReadsPerGene.out.tab
            ?????????Gene_xyz_KO_4_outputs_ReadsPerGene.out.tab
            ?????????treatment_vs_control.meta

> For the metadata file `treatment_vs_control.meta`, ensure that the row-names match the sample names of the individual RNA-seq samples.

          -------------   -------------   -------------
        |               | sampletype    | phenotype     |
        | ------------- | ------------- | ------------- |
        | Wildtype_1    | Control       | quiescent     |
        | ------------- | ------------- | ------------- |
        | Wildtype_2    | Control       | quiescent     |
        | ------------- | ------------- | ------------- |
        | Wildtype_3    | Control       | quiescent     |
        | ------------- | ------------- | ------------- |
        | Wildtype_4    | Control       | quiescent     |
        | ------------- | ------------- | ------------- |
        | Gene_xyz_KO_1 | Treatment     | proliferative |
        | ------------- | ------------- | ------------- |
        | Gene_xyz_KO_2 | Treatment     | proliferative |
        | ------------- | ------------- | ------------- |
        | Gene_xyz_KO_3 | Treatment     | proliferative |
        | ------------- | ------------- | ------------- |
        | Gene_xyz_KO_4 | Treatment     | proliferative |
          -------------   -------------   -------------

> Within the `config.yml` file, ensure that the order of sample naming matches the metadata file. This would be the order that the read count data is read into R and merged into the count `data.frame`. For the example metadata file above, the config file should look like:

```yaml
Paths:
  Treatment_vs_control_meta: "./data/treatment_vs_control.meta"
  Treatment_vs_control_data: "./data"

Samples:
  Treatment_vs_control:
    Wildtype_1: "Wildtype_1"
    Wildtype_2: "Wildtype_2"
    Wildtype_3: "Wildtype_3"
    Wildtype_4: "Wildtype_4"
    Gene_xyz_KO_1: "Gene_xyz_KO_1"
    Gene_xyz_KO_2: "Gene_xyz_KO_2"
    Gene_xyz_KO_3: "Gene_xyz_KO_3"
    Gene_xyz_KO_4: "Gene_xyz_KO_4"
```

<br>

## The main stages in the analysis of RNA sequencing data for this study are:

<br>

### 0. Merging STAR quantmode read counts

- The RNA seq count data was obtained by running STAR aligner with `--quantMode GeneCounts` enabled, which outputs the gene counts as a tab delimited file.
- Each file is then merged into a master data frame with its respective sample name as the column name.

  ### Dependencies:

  - [here](https://cran.r-project.org/web/packages/here/)

  - [yaml](https://cran.r-project.org/web/packages/yaml/)

<br>

### 1. Loading in count data, normalisation & QC

- The merged count data in a single data-frame is then compared against the `metadata.txt` file to ensure that its column names match the metedata row names. This is essential for DESeq2 to recognize the input dataframe as valid.
- Count distribution for each sample can also be plotted to QC checks on the read count distribution.
- The trend of variance vs mean is also plotted to show the overdispersion phenomenon in RNA-seq count data.
- An option for CPM filtering with edgeR is also included.
- The size factors for DESeq2 normalisation are calculated during the initializtion of the DESeq2 object and the DESeq2 normalised counts can be written out to disk.

  ### Dependencies:

  - [yaml](https://cran.r-project.org/web/packages/yaml/)

  - [tidyverse](https://www.tidyverse.org/)

  - [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

  - [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

  - [ggplot2](https://ggplot2.tidyverse.org/)

  - [pheatmap](https://cran.r-project.org/web/packages/pheatmap/)

  - [ggrepel](https://cran.r-project.org/web/packages/ggrepel/)

  - [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/)

<br>

### 2. DGE Hierarchical Clustering Heatmap & PCA

- This step performs log transformation of the count data for PCA and other correlative visualisations.
- Hierarchical clustering and a heatmap of pairwise correlation values of the gene wise counts per sample is performed to show how each sample relates to another. Ideally, samples within each condition should show the highest similarity to one another and most difference to those of the other condition.

  ### Dependencies:

  - [ggplot2](https://ggplot2.tidyverse.org/)

  - [pheatmap](https://cran.r-project.org/web/packages/pheatmap/)

  - [ggrepel](https://cran.r-project.org/web/packages/ggrepel/)

  - [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/)

<br>

### 3. Model fitting and hypothesis testing

- The differential expression testing with DESeq2's fitting of a general linear model takes place in this stage with `RNA_seq_dds <- DESeq(RNA_seq_dds)`. The dispersion estimates can subsequently be visualised.
- The contrast is defined to calculate the Log2FoldChange with the appropriate direction (i.e Treatment vs Control) and based on the dispersion estimates, the Log2FoldChange is shrunk appropriately.
- The MA plot can be visualised here as well as a quick summary of the differential expression results.
- After setting a Log2FoldChange threshold, the differentially expressed genes can be exported as a CSV file.

  ### Dependencies:

  - [here](https://cran.r-project.org/web/packages/here/)

  - [yaml](https://cran.r-project.org/web/packages/yaml/)

  - [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

  - [apeglm](https://bioconductor.org/packages/release/bioc/html/apeglm.html)

  - [dplyr](https://dplyr.tidyverse.org/)

  - [tidyverse](https://www.tidyverse.org/)

<br>

### 4. Visualization of DEGs

- The significant differentially expressed genes can be visualised as a bar-plot, volcano plot and a heatmap.

![](./images/heatmap.png)

- Above is an example of a heatmap of significant differentially expressed genes.

  ### Dependencies:

  - [ggplot2](https://ggplot2.tidyverse.org/)

  - [ggrepel](https://cran.r-project.org/web/packages/ggrepel/)

  - [ComplexHeatmap](https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)

  - [circlize](https://cran.r-project.org/web/packages/circlize/)

  - [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/)

  - [tidyverse](https://www.tidyverse.org/)

<br>

### 5. Post DESeq2 filtering and result annotation

This additional step post DESeq2 differential expression is mainly to accomplish these few objectives:
1. Filter out genes with less than a CPM of 0.2 and persist that data for downstream usage.
2. Annotate the filtered results with HGCN gene symbols, Entrez IDs & gene descriptions, and export as a csv.
3. Other filtering based on Log2FoldChange or Adjust p-value.
4. Export CPM edgeR normalized counts after CPM > 0.2 filter for [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp).

    ### Dependencies:

    - [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

    - [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

    - [dplyr](https://dplyr.tidyverse.org/)

    - [tidyverse](https://www.tidyverse.org/)

<br>

### 6. Functional enrichment analysis

The CPM filtered significant differentially expressed genes (sig_DEGs) will be piped into [ClusterProfiler](https://guangchuangyu.github.io/software/clusterProfiler/), a tool for performing functional enrichment analysis. ClusterProfiler houses a collection of tools and methods to perform a wholistic semantic analysis (drawing meaning from the data) of omics data by looking at global (and local) changes in the data and making relevant comparisons with publicly available curated omics and gene ontology databases.

A detailed walk-through of the major functionality of ClusterProfiler can be found [here](https://yulab-smu.top/biomedical-knowledge-mining-book/) as well as the code using the tools within ClusterProfiler as separate Rmd files [here](https://github.com/YuLab-SMU/biomedical-knowledge-mining-book).

For this workflow, we only use a number of the tools within ClusterProfiler that might be more commonly used in a typical RNA-seq analysis:

- `enrichGO()` function with:
  - Biological Processes (BP)
  - Molecular Function (MF)
  - Cellular Component (CC)

- Bar plots & dot plots of enrichment results for gene ontology over representation analysis (ORA).

- A number of disease-gene enrichment:
  - `enrichDGN()` Disease Gene Network
  - `enrichNCG()` Network of Cancer Genes
  - `enrichDO()` Disease Ontology

- Reactome pathway `enrichPathway()` and KEGG pathway `enrichKEGG()` over-representation analysis.

- Enrichment maps for gene sets from GO-BP, GO-MF & GO-CC.

- [More](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html) to come soon . . .

    ### Dependencies:

    - [tidyverse](https://www.tidyverse.org/)

    - [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)

    - [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

    - [DOSE](https://www.bioconductor.org/packages/release/bioc/html/DOSE.html)

    - [ReactomePA](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html)

    - [enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html)

    - [ggnewscale](https://cran.r-project.org/web/packages/ggnewscale/)

    - [ggplot2](https://ggplot2.tidyverse.org/)

    - [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

    - [msigdbr](https://cran.r-project.org/web/packages/msigdbr/)

<br>

## Some rules to follow for contributing to this code:

- All R scripts should be within `/src` and functions should be imported from there as to not clutter code.
- No data should be committed to git and file paths for data needed to be imported into R should be included in the `config.yaml` file. If file paths contain sensitive information, they can be added to gitignore.
- All analysis code should be in R markdown format as a `.rmd` file. Comments regarding analysis can be written as plain text in the R markdown file.
- renv.lock files can be committed to git for package version tracking and restoring a new renv on another system.

<br>

## Certain environment and user profile variables can be tweaked in the .Rprofile file but that file is gitignored for now.

Common issues can be tracked here along with their respective fixes:

- For the issue of not being able to fetch packages from Bioconductor, it might be attributed to a connection timeout. The fix for that is to add the `options(timeout=600)` for a longer timeout duration of 10 minutes. [Example](https://stackoverflow.com/questions/35282928/how-do-i-set-a-timeout-for-utilsdownload-file-in-r)

<br>

## TODOs:

- [x] Proper README and first PR.
- [x] In depth explanation of each step of the analysis.
- [x] Implement full functional enrichment analysis.
- [ ] Implement some version of Pathway Topology analysis (either NetGSA or SPIA or others if found).
- [ ] Comparative analysis with other pan cancer omics databases (e.g. TCGA).
