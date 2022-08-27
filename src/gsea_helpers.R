format2gsea <- function(normalised_counts_df,
                        gene_name_col = "NAME",
                        additional_col = "Description",
                        additional_col_contents = "NA"){
  # Get the number of rows from the input dataframe to specify as the length of the vector of "NA"s
  rows <- nrow(normalised_counts_df)
  # Create the vector of "NA" repeats to be n rows long.
  add_col <- rep(additional_col_contents, each = rows)
  # cbind that vector in front of the input dataframe so that it occurs on the left.
  normalised_counts_with_NA <- cbind(add_col, normalised_counts_df)
  
  # Name that column of "NA"s with "Description" or any other specified input from the
  # additional_col argument..
  colnames(normalised_counts_with_NA)[1] <- additional_col
  
  # Convert results to dataframe if not already in that format.
  normalised_counts_with_NA_df <- data.frame(normalised_counts_with_NA)
  
  # Convert rownames to the first column with tht column name as  "NAME" or
  # anything else specified along with the gene_name_col argument.
  normalised_counts_with_NA_no_rownames_df <- tibble::rownames_to_column(normalised_counts_with_NA_df,
                                                                         gene_name_col)
  
  return(normalised_counts_with_NA_no_rownames_df)
}