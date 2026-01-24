#' @title Validate Sequences from FASTA Files
#' @description This script validates a given list of FASTA files against a
#' series of checks (e.g., duplicates, invalid characters, ambiguities) and
#' generates summary and detailed reports.
#' @param fasta_files_to_validate A character vector of full paths to the FASTA
#' files that need to be validated.
#' @param output_dir The directory where the results (CSV files) will be stored.
#' @author Leon Balthaus, Gemini

# --- 1. FUNCTION DEFINITION ---

validate_sequences <- function(fasta_files_to_validate, output_dir) {
  
  # --- 1a. SETUP: LOAD PACKAGES ---
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required but not installed.")
  }
  library(Biostrings)
  
  # --- 1b. PRE-CHECKS: VALIDATE INPUTS ---
  
  # Check that the file list is valid and all files exist
  if (!is.character(fasta_files_to_validate) || length(fasta_files_to_validate) == 0) {
    stop("Error: 'fasta_files_to_validate' must be a non-empty character vector of file paths.")
  }
  files_not_found <- fasta_files_to_validate[!sapply(fasta_files_to_validate, file.exists)]
  if (length(files_not_found) > 0) {
    stop("Error: The following input FASTA files were not found:\n  ", paste(files_not_found, collapse = "\n  "))
  }
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    cat("Output directory '", output_dir, "' not found. Creating it.\n", sep = "")
    dir.create(output_dir, recursive = TRUE)
  }
  
  # --- 2. NESTED HELPER FUNCTION ---
  
  #' @description Performs validation checks on a DNAStringSet and returns a structured list.
  validate_sequence_set <- function(sequences, description) {
    cat(paste("--- Starting validation for:", description, "---\n"))
    
    if (length(sequences) == 0) {
      warning(paste("Validation Warning for", description, ": No sequences found. Skipping."))
      return(NULL)
    }
    cat(paste("    - Found", length(sequences), "sequence(s) to validate.\n"))
    
    # Check for valid names and duplicates
    if (any(is.null(names(sequences))) || any(nchar(names(sequences)) == 0)) warning(paste("Warning for", description, ": Missing names found."))
    if (any(duplicated(names(sequences)))) warning(paste("Warning for", description, ": Duplicated names found."))
    
    # Check for empty sequences
    seq_widths <- width(sequences)
    if (any(seq_widths == 0)) stop(paste("Error for", description, ": Zero-length sequences found."))
    cat(paste("    - Sequence lengths range from", min(seq_widths), "to", max(seq_widths), "bp.\n"))
    
    # Check for non-standard characters
    freq <- alphabetFrequency(sequences, baseOnly = FALSE, collapse = TRUE)
    valid_alphabet <- c(names(IUPAC_CODE_MAP), "-", "+", ".")
    present_chars <- names(freq[freq > 0])
    invalid_chars <- setdiff(present_chars, valid_alphabet)
    if (length(invalid_chars) > 0) {
      warning(paste("Warning for", description, ": Non-standard characters found:", paste(invalid_chars, collapse=", ")))
    }
    
    # Check for and report ambiguity codes
    cat("    - Checking for ambiguity codes...\n")
    ambiguity_codes <- c("M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "N")
    search_chars <- union(ambiguity_codes, invalid_chars)
    ambiguity_report <- data.frame(Sequence_Name=character(), Character=character(), Position=integer(), stringsAsFactors=FALSE)
    
    for (i in seq_along(sequences)) {
      for (char in search_chars) {
        locations <- unlist(gregexpr(char, as.character(sequences[[i]]), ignore.case = TRUE))
        if (locations[1] != -1) {
          for (pos in locations) {
            ambiguity_report <- rbind(ambiguity_report, data.frame(Sequence_Name=names(sequences)[i], Character=char, Position=pos))
          }
        }
      }
    }
    
    # Create summary data frame
    summary_df <- data.frame(
      Description = description,
      Num_Sequences = length(sequences),
      Min_Length = min(seq_widths),
      Max_Length = max(seq_widths),
      Ambiguity_Codes_Found = nrow(ambiguity_report),
      Ambiguity_Codes = if(nrow(ambiguity_report) > 0) paste(sort(unique(ambiguity_report$Character)), collapse = ", ") else "None"
    )
    
    cat(paste("    - Validation complete for:", description, "\n\n"))
    return(list(summary = summary_df, report = ambiguity_report))
  }
  
  # --- 3. MAIN LOGIC ---
  
  tryCatch({
    cat("--- Found", length(fasta_files_to_validate), "FASTA files to validate. ---\n\n")
    
    # Initialize lists to store results
    validation_summaries <- list()
    ambiguity_reports <- list()
    
    # Loop through each FASTA file provided in the list and validate it
    for (file_path in fasta_files_to_validate) {
      description <- tools::file_path_sans_ext(basename(file_path))
      sequences <- readDNAStringSet(file_path)
      
      result <- validate_sequence_set(sequences, description)
      
      if (!is.null(result)) {
        validation_summaries[[description]] <- result$summary
        if (nrow(result$report) > 0) {
          result$report$Description <- description
          ambiguity_reports[[description]] <- result$report
        }
      }
    }
    
    cat("--- All sequence sets loaded and validated successfully! ---\n")
    
    # Save the reports to CSV files
    cat("\n--- Saving validation results to CSV files... ---\n")
    
    final_summary_df <- do.call(rbind, validation_summaries)
    summary_csv_path <- file.path(output_dir, "validation_summary.csv")
    write.csv(final_summary_df, file = summary_csv_path, row.names = FALSE)
    cat(paste("    - Summary saved to '", summary_csv_path, "'.\n", sep = ""))
    
    if (length(ambiguity_reports) > 0) {
      final_ambiguity_df <- do.call(rbind, ambiguity_reports)
      ambiguity_csv_path <- file.path(output_dir, "ambiguity_details.csv")
      write.csv(final_ambiguity_df, file = ambiguity_csv_path, row.names = FALSE)
      cat(paste("    - Ambiguity details saved to '", ambiguity_csv_path, "'.\n", sep = ""))
    }
    
    return(invisible(TRUE))
    
  }, error = function(e) {
    cat("\n--------------------------------------------------\n")
    cat("An error occurred during the data validation process.\n")
    message(e)
    cat("--------------------------------------------------\n\n")
    return(invisible(FALSE))
  })
}