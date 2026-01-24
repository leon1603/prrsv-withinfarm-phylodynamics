#' @title Merge Infection Classification and Identity Data with Structured Metadata
#' @description This function reads a main metadata file and merges it with
#' infection classification files and a sequence identity summary file.
#' @param metadata_path The full path to the main structured metadata.csv file.
#' @param classification_files_list A character vector of full paths to the
#' infection classification CSV files.
#' @param identity_summary_path The full path to the identity_summary_by_pigID.csv file.
#' @param output_path The full path where the final merged CSV will be saved.
#' @author Leon Balthaus, Gemini

merge_data <- function(metadata_path, classification_files_list, identity_summary_path, output_path) {
  
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  
  # --- 1. Input Validation ---
  if (!file.exists(metadata_path)) {
    stop("Error: Metadata file not found at: ", metadata_path)
  }
  if (!file.exists(identity_summary_path)) {
    stop("Error: Identity summary file not found at: ", identity_summary_path)
  }
  
  files_exist <- sapply(classification_files_list, file.exists)
  if (!all(files_exist)) {
    missing_files <- classification_files_list[!files_exist]
    stop("Error: The following classification files were not found:\n",
         paste(missing_files, collapse = "\n"))
  }
  
  cat("--- Starting metadata merge process ---\n")
  
  # --- 2. Load and Prepare Main Metadata ---
  metadata_df <- read_csv(
    metadata_path,
    col_types = cols(pigID_weekOfSampling = col_character()),
    show_col_types = FALSE
  ) %>%
    # Ensure new_infection is a logical column for reliable checks
    mutate(new_infection = as.logical(new_infection))
  
  cat(paste(" -> Successfully loaded main metadata with", nrow(metadata_df), "rows.\n"))
  
  # --- 3. Process and Merge Classification Files ---
  for (file_path in classification_files_list) {
    file_name <- basename(file_path)
    cat(paste(" -> Processing file:", file_name, "\n"))
    
    analysis_type <- case_when(
      str_detect(file_name, "WGS") ~ "WGS",
      str_detect(file_name, "ORF5") ~ "ORF5",
      TRUE ~ "unknown"
    )
    
    if (analysis_type == "unknown") {
      warning(paste("Could not determine analysis type for", file_name, "(or type has been excluded). Skipping."))
      next
    }
    
    classification_data <- read_csv(file_path, col_types = cols(pigID_weekOfSampling = col_character()), show_col_types = FALSE)
    
    classification_to_merge <- classification_data %>%
      select(pigID_weekOfSampling, Clade, Classification)
    
    names(classification_to_merge)[names(classification_to_merge) == "Clade"] <- paste0(analysis_type, "_clade")
    names(classification_to_merge)[names(classification_to_merge) == "Classification"] <- paste0(analysis_type, "_classification")
    
    # Use full_join to keep rows that are in classification_data but not yet in metadata_df
    metadata_df <- metadata_df %>%
      full_join(classification_to_merge, by = "pigID_weekOfSampling")
  }
  
  cat(" -> Classification merge complete. Expanding metadata to include all found samples.\n")
  
  # --- 4. Load and Merge Identity Summary Data ---
  identity_summary_df <- read_csv(identity_summary_path, show_col_types = FALSE)
  
  # Exclude the summary row and explicitly remove nsp2/nsp9 columns
  identity_data_to_merge <- identity_summary_df %>%
    filter(pigID_weekOfSampling != "summary") %>%
    select(-matches("nsp2"), -matches("nsp9"))
  
  metadata_df <- metadata_df %>%
    left_join(identity_data_to_merge, by = "pigID_weekOfSampling")
  
  cat(paste(" -> Merged sequence identity data.\n"))
  
  # --- 5. Final Data Processing ---
  classification_cols <- names(metadata_df)[str_detect(names(metadata_df), "_classification$")]
  clade_cols_to_convert <- names(metadata_df)[str_detect(names(metadata_df), "_clade$")]
  
  if (length(clade_cols_to_convert) > 0) {
    metadata_df <- metadata_df %>%
      mutate(across(all_of(clade_cols_to_convert), as.character))
  }
  
  all_new_cols <- c(clade_cols_to_convert, classification_cols)
  
  metadata_processed <- metadata_df %>%
    mutate(across(all_of(all_new_cols), ~tidyr::replace_na(., "none")))
  
  # --- Move Notes to the last column ---
  metadata_processed <- metadata_processed %>%
    relocate(matches("(?i)^notes$"), .after = last_col())
  
  # --- 7. Save Final Output ---
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  
  write_csv(metadata_processed, output_path)
  cat(paste(" -> Merge complete. Final metadata saved to:", output_path, "\n"))
  cat("--- Process finished successfully! ---\n")
  
  return(invisible(TRUE))
}