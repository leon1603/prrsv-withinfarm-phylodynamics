#' @title Perform Multiple Sequence Alignment on Specific Files
#' @description This function accepts a list of specific FASTA files, performs a
#' separate multiple sequence alignment  on each file using the DECIPHER
#' package, and saves the resulting aligned FASTA files to a specified output directory.
#' @param input_files A character vector containing the full paths to the FASTA
#' files to be aligned.
#' @param output_dir The directory where all aligned FASTA files will be saved.
#' If it doesn't exist, it will be created.
#' @author Leon Balthaus, Gemini

# --- 1. FUNCTION DEFINITION ---

perform_msa <- function(input_files, output_dir) {
  
  # --- 1a. SETUP: Load Packages and Validate Inputs ---
  if (!requireNamespace("DECIPHER", quietly = TRUE) || !requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Packages 'DECIPHER' and 'Biostrings' are required. Please install them.")
  }
  library(DECIPHER)
  library(Biostrings)
  
  # Check if input is valid
  if (missing(input_files) || length(input_files) == 0) {
    stop("Error: 'input_files' must be a non-empty vector of file paths.")
  }
  
  # Check if files actually exist
  missing_files <- input_files[!file.exists(input_files)]
  if (length(missing_files) > 0) {
    stop("Error: The following input files were not found:\n  ", paste(missing_files, collapse = "\n  "))
  }
  
  # --- 1b. Create Output Directories ---
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # --- 2. MAIN LOGIC ---
  tryCatch({
    cat("--- Starting Multiple Sequence Alignment on Selected Files ---\n")
    cat(sprintf("-> Found %d files to process.\n", length(input_files)))
    cat(sprintf("-> Saving results to: %s\n\n", output_dir))
    
    # --- Step 2b: Loop Through Each File, Align, and Save ---
    for (file_path in input_files) {
      file_basename <- basename(file_path)
      cat(sprintf("--- Processing file: %s ---\n", file_basename))
      
      # Load sequences from the current file
      sequences_to_align <- readDNAStringSet(file_path)
      cat(sprintf(" -> Found %d sequences to align.\n", length(sequences_to_align)))
      
      # Perform the alignment
      cat(" -> Aligning sequences (this may take a moment)...\n")
      aligned_sequences <- AlignSeqs(sequences_to_align, verbose = TRUE)
      cat(" -> Alignment complete.\n")
      
      # --- Step 2c: Save Outputs for the Current File ---
      
      # Define new output filenames
      base_name_sans_ext <- tools::file_path_sans_ext(file_basename)
      output_fasta_name <- paste0(base_name_sans_ext, "_Aligned.fasta")
      
      # Set full path for saving
      output_fasta_path <- file.path(output_dir, output_fasta_name)
      
      # Save files
      writeXStringSet(aligned_sequences, file = output_fasta_path)
      cat(sprintf(" -> Saved aligned FASTA to: %s\n", output_fasta_path))
      
      cat("\n")
    }
    
    cat("--- All MSA tasks completed successfully. ---\n")
    return(invisible(TRUE))
    
  }, error = function(e) {
    # --- Error Handling ---
    cat("\n--------------------------------------------------\n")
    cat("An error occurred during the MSA process.\n")
    message(e)
    cat("--------------------------------------------------\n\n")
    return(invisible(FALSE))
  })
}