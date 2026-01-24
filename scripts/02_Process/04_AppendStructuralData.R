#' @title Insert Pen and Room into FASTA Headers with Batch Filtering
#' @description This function reads a FASTA file and a metadata CSV, matches
#' sequences to pen and room numbers based on the pig ID in the header, and
#' writes a new FASTA file.
#'
#' It allows filtering sequences by "Batch" (e.g., L1 vs L3) found in the header
#' and supports custom output filenames.
#'
#' @param metadata_csv_path The full path to the metadata CSV file.
#' @param input_fasta_path The full path to the input FASTA file.
#' @param output_dir The base path for the output directory.
#' @param batch_to_remove (Optional) A character string indicating a batch to
#' remove (e.g., "L1" or "L3"). The function looks for this string surrounded
#' by underscores in the header (e.g., "_L1_") and removes those sequences.
#' @param custom_output_filename (Optional) The specific name for the output file
#' (e.g., "MyFilteredData.fasta"). If NULL, a default name is generated.
#'
#' @return Creates a new FASTA file in the output_dir.
#' @author Gemini
#' @date 2025-10-23 (Updated 2025-11-24)

# --- 1. FUNCTION DEFINITION ---

append_pen_and_room <- function(metadata_csv_path,
                                input_fasta_path,
                                output_dir,
                                batch_to_remove = NULL,
                                custom_output_filename = NULL) {
  
  # --- 2. Load Libraries ---
  required_packages <- c("dplyr", "readr", "Biostrings", "stringr", "tools")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(paste("Package", pkg, "is needed."))
    library(pkg, character.only = TRUE)
  }
  
  # --- 3. Main Logic ---
  tryCatch({
    
    cat("--- Starting: Insert Pen and Room into FASTA Headers ---\n")
    
    # --- 3.1. Validate inputs ---
    if (!file.exists(metadata_csv_path)) stop(paste("Metadata file not found:", metadata_csv_path))
    if (!file.exists(input_fasta_path)) stop(paste("Input FASTA file not found:", input_fasta_path))
    
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    cat(paste(" -> Output will be saved to:", output_dir, "\n"))
    
    # --- 3.2. Load Metadata (Pen and Room) ---
    cat(paste(" -> Loading metadata from:", metadata_csv_path, "\n"))
    
    meta_map <- read_csv(metadata_csv_path, show_col_types = FALSE) %>%
      select(pigID_weekOfSampling, pen, room) %>%
      distinct() %>%
      mutate(
        pen = as.character(pen),
        room = as.character(room)
      )
    
    cat(paste(" -> Metadata loaded.", nrow(meta_map), "unique pig/pen/room combinations found.\n"))
    
    # --- 3.3. Load FASTA ---
    cat(paste(" -> Loading FASTA sequences from:", input_fasta_path, "\n"))
    sequences <- readDNAStringSet(input_fasta_path)
    all_headers <- names(sequences)
    cat(paste(" -> Loaded", length(sequences), "sequences.\n"))
    
    # --- 3.4. Extract Pig IDs ---
    header_df <- data.frame(OriginalHeader = all_headers, stringsAsFactors = FALSE)
    # Extracts ID between underscores (e.g., in '..._518_0_...')
    header_df$pigID_weekOfSampling <- str_match(all_headers, "_(\\d+_\\d+)_")[, 2]
    
    num_extracted <- sum(!is.na(header_df$pigID_weekOfSampling))
    cat(paste(" -> Extracted", num_extracted, "pig IDs from headers.\n"))
    
    # --- 3.5. Join with Metadata ---
    header_df <- header_df %>%
      left_join(meta_map, by = "pigID_weekOfSampling") %>%
      mutate(
        pen_str = ifelse(is.na(pen), "Unknown", pen),
        room_str = ifelse(is.na(room), "Unknown", room)
      )
    
    num_matched <- sum(header_df$pen_str != "Unknown")
    cat(paste(" -> Matched", num_matched, "sequences to metadata.\n"))
    
    # --- 3.6. Create New Headers ---
    parts <- str_split_fixed(header_df$OriginalHeader, "_", 2)
    part1 <- parts[, 1]
    part2 <- parts[, 2]
    has_underscore <- str_detect(header_df$OriginalHeader, "_")
    
    # Construct: [Part1]_P[Pen]_R[Room]_[Part2]
    header_df$NewHeader <- ifelse(
      has_underscore,
      paste0(part1, "_P", header_df$pen_str, "_R", header_df$room_str, "_", part2),
      paste0(part1, "_P", header_df$pen_str, "_R", header_df$room_str)
    )
    
    names(sequences) <- header_df$NewHeader
    cat(" -> New headers created (Pen and Room appended).\n")
    
    # --- 3.7. Filter by Batch (NEW) ---
    if (!is.null(batch_to_remove)) {
      cat(paste(" -> Filtering: Removing sequences belonging to batch", shQuote(batch_to_remove), "\n"))
      
      # We look for "_L1_" or "_L3_" specifically to avoid partial matches
      batch_pattern <- paste0("_", batch_to_remove, "_")
      
      # Determine which headers contain the batch pattern
      # Note: We check OriginalHeader to be safe, though NewHeader works too
      has_batch <- str_detect(header_df$OriginalHeader, batch_pattern)
      
      # Keep only those that DO NOT match the batch
      keep_index <- !has_batch
      
      num_removed <- sum(has_batch)
      
      if (num_removed > 0) {
        sequences <- sequences[keep_index]
        cat(paste(" -> Removed", num_removed, "sequences matching batch", batch_to_remove, "\n"))
        cat(paste(" ->", length(sequences), "sequences remaining.\n"))
      } else {
        cat(paste(" -> Warning: Batch", shQuote(batch_to_remove), "was specified but not found in any headers.\n"))
      }
    }
    
    # --- 3.8. Save Output (Updated Naming) ---
    
    if (!is.null(custom_output_filename)) {
      # Use the custom name provided by the user
      output_filename <- custom_output_filename
    } else {
      # Fallback to automatic naming
      fasta_basename <- tools::file_path_sans_ext(basename(input_fasta_path))
      output_filename <- paste0(fasta_basename, "_AppendStructuralData.fasta")
    }
    
    output_fasta_path <- file.path(output_dir, output_filename)
    writeXStringSet(sequences, output_fasta_path)
    
    cat(paste(" -> Saved to:", output_fasta_path, "\n"))
    cat("\n--- Completed Successfully ---\n")
    return(invisible(TRUE))
    
  }, error = function(e) {
    cat(paste("\n--- ERROR:", e$message, "\n"))
    return(invisible(FALSE))
  })
}