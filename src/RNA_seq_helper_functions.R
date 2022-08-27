STARo2df<-function(sample_names, samples_path,
                   gene_col=1, count_col=2, gene_naming="Ensembl",
                   file_type="\\.tab$", start=5){
  
  # Obtain individual data file names in full data paths
  star_counts_full_paths <- (list.files(path=samples_path, pattern=file_type, full.names=TRUE))
  
  # Check if number of samples provided is equal or less than the data source.
  if (length(star_counts_full_paths) == length(sample_names)){
    message("Same number of samples as data in sample data path. \n")
  } else if (length(star_counts_full_paths) > length(sample_names)){
    message("Fewer samples specified than data in sample data path. \n")
  }
  
  samples_to_process <- list()
  
  # Matches the sample names to the file paths for the sample data.
  for (sample_path in star_counts_full_paths){
    for (sample_name in sample_names){
      if (grepl(sample_name, sample_path)){
        names(sample_path) <- sample_name
        samples_to_process <- c(samples_to_process, sample_path)
      }
    }
  }
  
  count_subset_list <- list()
  genes_from_data <- list()
  
  idx <- 1
  
  for (sample in samples_to_process){
    # Iterate through samples and read TSV into R, subset rows of interest from default arguments.
    count_data <- read.table(sample)
    n_rows <- nrow(count_data)
    # Subset the columns of interest with default arguments.
    counts_subset <- data.frame(count_data[start: n_rows, gene_col:count_col])
    message("Loaded data from <", names(samples_to_process)[idx], "> and subsetting columns of interest.\n")
    # Name the column of gene names with the default argument.
    colnames(counts_subset)[gene_col] <- gene_naming
    # Name the gene count column of the sub-setted data-frame with the sample name.
    colnames(counts_subset)[count_col] <- names(samples_to_process)[idx]
    count_subset_list <- append(count_subset_list, list(counts_subset))
    genes_from_data <- append(genes_from_data, counts_subset[gene_col])
    idx <- idx + 1
  }
  
  # Check if all input data have the same set of genes by ensuring that the gene list are identical,
  # hence length(unique(x)) == 1.
  if (length(unique(genes_from_data)) == 1){
    message("Genes from all input data are a match. \n")
    STAR_samples_merged <- plyr::join_all(count_subset_list,
                                          by = gene_naming,
                                          type = "full",
                                          match = "all")
    
  } else {
    cat(crayon::red("Warning: Genes from certain input data are missing. \n"))
    cat(crayon::red("Warning: Performing reduce opreation to obtain list of unique genes. \n"))
    message()
    
    unique_genes <- Reduce(union, genes_from_data)
    
    STAR_samples_merged <- plyr::join_all(count_subset_list,
                                          by = gene_naming,
                                          type = "full",
                                          match = "all")
    # Check if joined dataframe has the same number of genes as
    # the union of all the unique genes.
    if (length(unique_genes) == nrow(STAR_samples_merged)){
      message("Unqiue genes successfully matches output dataframe size. \n")
    } else{
      # if not break and exit the function.
      cat(crayon::red("Warning: Missing genes cannot be successfully mapped! \n"))
      break
    }
  }
  return(STAR_samples_merged)
}
