library(tidyverse)

##### DOWNLOAD DATA FROM ONCOKB
#     https://oncodb.org/data_download.html
#     It will download all the cancer download data.

# it will extract Smoking data from clinical and then merge
merge_smoking_data <- function(fnames, output_file) {
  # Check if exactly three files are provided
  if (length(fnames) != 3) stop("Exactly three file names must be provided.")
  
  # Check if files exist
  if (!all(file.exists(fnames))) stop("One or more files do not exist: ", paste(fnames[!file.exists(fnames)], collapse = ", "))
  
  # Define prefixes and filtering conditions
  prefixes <- c("ClinExp_", "ClinMeth_", "DifExp_")
  filter_smoking <- c(TRUE, TRUE, FALSE) # Filter 'Smoking' for first two files only
  
  # Read and process each file
  dfs <- map2(fnames, seq_along(fnames), function(file, idx) {
    tryCatch({
      df <- read_tsv(file, show_col_types = FALSE)
      
      # Filter for 'Smoking' if required
      if (filter_smoking[idx]) {
        df <- df %>% filter(`Clinical parameters` == "Smoking")
      }
      
      # Sanitize and rename columns
      df <- df %>% rename_with(~paste0(prefixes[idx], make.names(.x)), 
                               -c(`NCBI gene id`, `Gene symbol`, `Cancer type`))
      
      # Ensure required columns exist
      required_cols <- c("NCBI gene id", "Gene symbol", "Cancer type")
      missing_cols <- setdiff(required_cols, colnames(df))
      if (length(missing_cols) > 0) {
        stop(sprintf("Missing columns %s in %s", paste(missing_cols, collapse = ", "), file))
      }
      
      df
    }, error = function(e) {
      stop("Error processing ", file, ": ", e$message)
    })
  })
  
  # Merge data frames on specified columns
  merged_df <- tryCatch({
    dfs[[1]] %>%
      inner_join(dfs[[2]], by = c("NCBI gene id", "Gene symbol", "Cancer type")) %>%
      inner_join(dfs[[3]], by = c("NCBI gene id", "Gene symbol", "Cancer type"))
  }, error = function(e) {
    stop("Error merging data frames: ", e$message)
  })
  
  # Save merged data to output file
  write_tsv(merged_df, output_file)
  
  # Return merged data frame
  merged_df
}

# LUSC
fnames <- c(
  "clinical/Differential_Clinical_Expression_LUSC.txt",
  "clinical/Differential_Clinical_Methylation_LUSC.txt",
  "expression/LUSC_Differential_Gene_Expression_Table.txt"
)
merged_data <- merge_smoking_data(fnames, output_file = 'lusc.txt')
# LUAD
fnames <- c(
  "clinical/Differential_Clinical_Expression_LUAD.txt",
  "clinical/Differential_Clinical_Methylation_LUAD.txt",
  "expression/LUAD_Differential_Gene_Expression_Table.txt"
)
merged_data <- merge_smoking_data(fnames, output_file = 'luad.txt')
# CESC
fnames <- c(
  "clinical/Differential_Clinical_Expression_CESC.txt",
  "clinical/Differential_Clinical_Methylation_CESC.txt",
  "expression/CESC"
)
merged_data <- merge_smoking_data(fnames, output_file = 'cesc.txt')
