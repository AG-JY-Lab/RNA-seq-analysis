enrichedGenes <- function(enrichment_results, enrichment_term,
                          database = 'org.Hs.eg.db',
                          genetype = 'ENSEMBL',
                          delimiter = "/",
                          set_readable = FALSE){
  # If the set_readable argument is TRUE, the enrichment results have their IDs
  # converted to HGNC symbols with DOSE::setReadable() function.
  if (set_readable){
    enrichment_results_readable <- DOSE::setReadable(x = enrichment_results,
                                                     OrgDb = database,
                                                     keyType = genetype)
  } else {
    enrichment_results_readable <- enrichment_results
  }
  # Transform the enrichment results into a dataframe.
  enrichment_df <- data.frame(enrichment_results_readable)
  # Isolate the row which the enrichment dataframe cell in the Description column
  # matches the specified enrichment term argument and obtain all the genes in that row.
  enrichment_genes <- enrichment_df[enrichment_df$Description == enrichment_term, ]$geneID
  # Convert that single string of gene IDs into a vector of individual gene IDs by
  # splitting them with the delimiter, which defaults to "/"
  enrichment_genes_vector <- as.vector(unlist(strsplit(enrichment_genes, delimiter)))
  return(enrichment_genes_vector)
}


gene_set_heamap <- function(gene_vector, count_df,
                            gene_set_name = "Gene set XYZ",
                            experiment_name = "Treatment condition ABC vs Control",
                            gene_naming = "Ensembl ID",
                            map_to_HGNC_symbol = FALSE){

  # Convert the necessary column of the gene count dataframe to the row names.
  # This defaults to the column named "Ensembl ID".
  count_df_rn <- column_to_rownames(count_df, var = gene_naming)
  
  # Filter the rows within the gene count dataframe whose row names match the gene IDs
  # in the gene_vector argument, which is a character vector of gene IDs.
  heatmap_counts <- count_df_rn %>% dplyr::filter(row.names(count_df_rn) %in% gene_vector) %>%
    data.frame()
  
  # Quick-fix
  # TODO improve this Ensembl ID to HGNC mapping
  if (map_to_HGNC_symbol){
    # Get a mapping of the Ensembl IDs to HGNC symbols returned as a non unique
    # character vector of HGNC symbols which is the same length as the dataframe.
    symbols <- mapIds(org.Hs.eg.db, keys = row.names(heatmap_counts), keytype = "ENSEMBL", column = "SYMBOL")
    # Assign that to a new column in the count dataframe named Gene_symbol.
    heatmap_counts$Gene_symbol <- symbols
    # Consolidate the gene counts by grouping counts by HGNC symbols and taking
    # the sum of counts per HGNC symbol.
    heatmap_counts_consolidated <- {heatmap_counts %>% group_by(Gene_symbol) %>% summarise_all(sum)}
    # Remove the original row names of the count dataframe and add the HGNC
    # symbol column as the new row names.
    heatmap_counts <- {heatmap_counts_consolidated %>% remove_rownames %>% column_to_rownames(var = "Gene_symbol")}
  }

  # Apply the scaling function to the rows of the heat map counts dataframe to
  # normalize them to be suitable for the heat map.
  heatmap_counts_matrix_scaled = t(apply(X = heatmap_counts, MARGIN = 1, FUN = scale))
  # Calculate the row means, which would be the average gene count.
  base_mean <- rowMeans(heatmap_counts)
  
  # The heat map annotation would be the column names of the count dataframe,
  # i.e. the sample names.
  heatmap_annotation <- ComplexHeatmap::HeatmapAnnotation(Sample = colnames(count_df_rn))
  
  # The heat map title a combination of the experiment name and the gene set.
  heatmap_title <- paste(experiment_name, gene_set_name, sep = "\n")
  
  # Create the heat map list object with the scaled count matrix.
  # The heat map_title variable would be the title for the heat map.
  heatmap_list <- ComplexHeatmap::Heatmap(as.matrix(heatmap_counts_matrix_scaled),
                                          column_title = heatmap_title,
                                          name = "Relative \nexpression",
                                          # A green to red to black to green colour palette is chosen.
                                          col = circlize::colorRamp2(c(-2, 0, 2), c("green", "black", "red")),
                                          top_annotation = heatmap_annotation,
                                          # Row names are hidden to not conflict with the base mean heat map.
                                          show_row_names = FALSE) +
    # Join the heat map of the gene base mean to the main heat map.
    ComplexHeatmap::Heatmap(base_mean,
                            name = "Gene counts \nbase mean",
                            show_row_names = TRUE,
                            show_column_names = FALSE)
  
  return(heatmap_list)
}


gene_frequency <- function(enrichment_result,
                           delimiter = "/"){
  # Convert enrichment results into dataframe
  enirchment_result_df <- data.frame(enrichment_result)
  # Extract the gene ID column from the enrichment result dataframe
  enirchment_geneID_column <- enirchment_result_df$geneID
  # Initialize and empty vector
  all_genes <- vector()
  # Iterate through the column of merged Ensembl IDS delimited by "/"
  # ENSG00000037280/ENSG00000066056/ENSG00000106178/ENSG00000112715
  for (ensembl_string in enirchment_geneID_column) {
    # Split the strings by "/" or whichever relevant delimiter
    # Obtain output of vector of Ensembl (or other) gene IDS
    enrichment_genes_vector <- as.vector(unlist(strsplit(ensembl_string, delimiter)))
    # Append vector of Ensembl (or other) gene IDS to main gene vector
    all_genes <- c(all_genes, enrichment_genes_vector)
  }
  # Count the occurrences of genes in the main gene vector
  all_genes_table <- table(all_genes)
  # Return the gene frequencies
  return(all_genes_table)
}
