#' @title Reformat, Split, and Correct Dates in FASTA Sequences
#' @description Processes FASTA data to extract genomic regions, reformat headers, 
#' split by batch/region, and optionally correct sampling dates for specific batches/weeks/regions.
#' @param input_data A character vector of file paths or a named list of sequence data.
#' @param output_dir The directory where all return files will be saved.
#' @param target_batch (Optional) String specifying the batch ID to target for date correction (e.g., "L1").
#' @param target_week (Optional) String specifying the week to target for date correction (e.g., "0").
#' @param target_region (Optional) String specifying the genomic region to target ("WGS" or "ORF5").
#' @param new_date (Optional) String for the new date in 'YYYY-MM-DD' format.
#' @author Leon Balthaus, Gemini

reformat_and_split_sequences <- function(input_data, 
                                         output_dir, 
                                         target_batch = NULL, 
                                         target_week = NULL, 
                                         target_region = NULL,
                                         new_date = NULL) {
  
  # --- Auto-detect input type ---
  if (is.character(input_data) && !is.list(input_data)) {
    cat("Detected vector of file paths. Converting to named list...\n")
    sequence_list <- as.list(input_data)
    names(sequence_list) <- basename(input_data)
  } else {
    sequence_list <- input_data
  }
  
  # Validate input structure
  if (!is.list(sequence_list) || is.null(names(sequence_list))) {
    stop("Error: Input must be either a character vector of file paths OR a named list.")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Validate Date Correction Arguments
  perform_correction <- FALSE
  if (!is.null(target_batch) || !is.null(target_week) || !is.null(new_date) || !is.null(target_region)) {
    # Check completeness (Region is technically optional here, but usually required for specific targeting)
    if (is.null(target_batch) || is.null(target_week) || is.null(new_date)) {
      warning("Date correction arguments (batch, week, or date) incomplete. Skipping date correction.")
    } else {
      if (!grepl("^\\d{4}-\\d{2}-\\d{2}$", new_date)) {
        stop("Error: 'new_date' must be in 'YYYY-MM-DD' format.")
      }
      
      # Validate target_region if provided
      if (!is.null(target_region) && !target_region %in% c("WGS", "ORF5")) {
        stop("Error: 'target_region' must be either 'WGS' or 'ORF5'.")
      }
      
      perform_correction <- TRUE
      target_week <- as.character(target_week)
      
      msg_region <- if(is.null(target_region)) "ALL Regions" else target_region
      cat(sprintf("--- Date Correction Active: Changing Batch %s Week %s (Region: %s) to %s ---\n", 
                  target_batch, target_week, msg_region, new_date))
    }
  }
  
  cat("Found", length(sequence_list), "items to process...\n\n")
  
  # Initialize lists
  sequences_by_batch <- list(Combined = list(), b1 = list(), b3 = list())
  correction_counter <- 0
  
  for (item_name in names(sequence_list)) {
    
    # 1. Region Extraction
    clean_name_for_region <- sub("^Combined_", "", item_name)
    genomic_region <- strsplit(clean_name_for_region, "_")[[1]][1]
    
    cat(" -> Processing '", item_name, "' (Region: '", genomic_region, "')... ", sep = "")
    
    # 2. Content Reading
    content <- sequence_list[[item_name]]
    if (length(content) == 1 && file.exists(content[1])) {
      fasta_lines <- readLines(content[1], warn = FALSE)
      cat("(Read from file)\n")
    } else {
      fasta_lines <- content
      cat("(Processed raw lines)\n")
    }
    
    # 3. Parsing Logic
    header_indices <- grep(">", fasta_lines)
    if (length(header_indices) == 0) {
      warning(paste("\n   No headers found in", item_name, "- skipping."))
      next
    }
    
    sequence_boundaries <- c(header_indices, length(fasta_lines) + 1)
    
    for (i in 1:(length(sequence_boundaries) - 1)) {
      start_index <- sequence_boundaries[i]
      end_index <- sequence_boundaries[i + 1] - 1
      
      original_header <- fasta_lines[start_index]
      if (start_index >= end_index) { sequence_block <- "" } 
      else { sequence_block <- fasta_lines[(start_index + 1):end_index] }
      
      clean_header <- sub(">", "", original_header)
      
      # Regex extraction
      accession_number <- strsplit(clean_header, " |\\|")[[1]][1]
      isolate_match <- regexec("BA-(L\\d+)-(\\d+)-(\\d+)w", clean_header)
      isolate_parts <- regmatches(clean_header, isolate_match)
      all_fields <- strsplit(clean_header, "\\|")[[1]]
      
      date_of_sampling <- NA
      if (length(all_fields) >= 2) {
        date_of_sampling <- trimws(all_fields[length(all_fields) - 1])
      }
      
      # --- CHECK PARSING SUCCESS ---
      if (length(isolate_parts[[1]]) == 4 && 
          !is.na(date_of_sampling) && nzchar(date_of_sampling) && 
          !is.na(accession_number) && nzchar(accession_number)) {
        
        batch <- isolate_parts[[1]][2]
        pigID <- isolate_parts[[1]][3]
        week  <- isolate_parts[[1]][4]
        
        # --- DATE CORRECTION ---
        if (perform_correction) {
          # Check if the region matches (only if target_region is provided)
          region_matches <- is.null(target_region) || (genomic_region == target_region)
          
          if (region_matches && batch == target_batch && week == target_week) {
            date_of_sampling <- new_date
            correction_counter <- correction_counter + 1
          }
        }
        
        new_header <- paste(
          accession_number, batch, pigID, week, genomic_region, date_of_sampling, sep = "_"
        )
        new_header <- paste0(">", new_header)
        
      } else {
        new_header <- original_header
      }
      
      formatted_chunk <- c(new_header, sequence_block)
      
      # Add to Combined
      sequences_by_batch$Combined[[genomic_region]] <- c(sequences_by_batch$Combined[[genomic_region]], formatted_chunk)
      
      # Add to Specific Batches (based on original header presence of L1/L3)
      if (grepl("L1", original_header)) {
        sequences_by_batch$b1[[genomic_region]] <- c(sequences_by_batch$b1[[genomic_region]], formatted_chunk)
      }
      if (grepl("L3", original_header)) {
        sequences_by_batch$b3[[genomic_region]] <- c(sequences_by_batch$b3[[genomic_region]], formatted_chunk)
      }
    }
  }
  
  if (perform_correction) {
    cat(sprintf("\n--- Date Correction Summary: Updated %d sequences ---\n", correction_counter))
  }
  
  # --- Write consolidated files ---
  cat("\n--- Writing consolidated FASTA files ---\n")
  for (batch_name in names(sequences_by_batch)) {
    for (region_name in names(sequences_by_batch[[batch_name]])) {
      sequences_to_write <- sequences_by_batch[[batch_name]][[region_name]]
      
      if (length(sequences_to_write) > 0) {
        output_filename <- sprintf("%s_%s_Cohort_Formatted_Aligned.fasta", batch_name, region_name)
        output_filepath <- file.path(output_dir, output_filename)
        writeLines(sequences_to_write, output_filepath)
        cat(sprintf(" -> Saved file to: %s\n", output_filepath))
      }
    }
  }
  
  return(invisible(TRUE))
}