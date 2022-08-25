# RNA-seq end to end analysis in R

### This repository maintains and end to end downstream workflow of bulk RNA sequencing analysis from gene read counts to functional enrichment analysis.

<br>

## Before you begin:
> ### Ensure that you have [renv](https://rstudio.github.io/renv/index.html) installed as this will the main environment management tool used in this workflow. The packages used in the analysis can be easily obtained with `renv::restore()` and this will download them by referencing the `renv.lock` file.

> ### Although this is not a standalone package for bulk RNA-seq analysis but more of a pipeline that integrates multiple tools together, there are some helper functions that were written as R source code within the `/src` directory. In order to load them for use in the analysis, run `source(here::here("src", "RNA_seq_helper_functions.R"))` at the start of each Rmd notebook to ensure the necessary helper functions are loaded too.

> ### The read count data and metadata are expected to come from the outputs of [STAR](https://github.com/alexdobin/STAR). As such the format for the data would be a a number of `_outputs_ReadsPerGene.out.tab` files which should be place in the `/data` folder and a tab delimited metadata text while which the user has to write for DESeq2 to build it's deseq_object.
<br>

    ──RNA-seq-analysis
        └──data
            ├──Control_1_outputs_ReadsPerGene.out.tab
            ├──Control_2_outputs_ReadsPerGene.out.tab
            ├──Control_3_outputs_ReadsPerGene.out.tab
            ├──Control_4_outputs_ReadsPerGene.out.tab
            ├──Treatment_1_outputs_ReadsPerGene.out.tab
            ├──Treatment_2_outputs_ReadsPerGene.out.tab
            ├──Treatment_3_outputs_ReadsPerGene.out.tab
            ├──Treatment_4_outputs_ReadsPerGene.out.tab
            └──treatment_vs_control.meta
> ### For the metadata file `treatment_vs_control.meta`, ensure that the row-names match the sample names of the individual RNA-seq samples.
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
        | Gene_xyz_KO_  | Treatment     | proliferative |
        | ------------- | ------------- | ------------- |
        | Treatment_3   | Treatment     | proliferative |
        | ------------- | ------------- | ------------- |
        | Treatment_4   | Treatment     | proliferative |
          -------------   -------------   ------------- 



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

- Some explanation

<br>

### 5. Post DESeq2 filtering and thresholding

- Some explanation

<br>

### 6. Functional enrichment analysis

- Some explanation

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
- [ ] In depth explanation of each step of the analysis.
- [ ] Implement full functional enrichment analysis.
- [ ] Implement some version of Pathway Topology analysis (either NetGSA or SPIA or others if found).
- [ ] Comparative analysis with other pan cancer omics databases (e.g. TCGA).
