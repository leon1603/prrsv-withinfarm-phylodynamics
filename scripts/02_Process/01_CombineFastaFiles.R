#' @title Combine and Filter Multiple FASTA Files
#' @description This function reads content from multiple FASTA files, combines them into a single output FASTA file, and filters out duplicate sequences.
#' @param input_files A character vector containing the file paths of the FASTA files to be combined.
#' @param output_file A string specifying the file path for the combined output FASTA file.
#' @author Leon Balthaus, Gemini

# --- 1. FUNCTION DEFINITION ---

combine_fasta_files <- function(input_files, output_file) {
  # A vector to keep track of pigID_week combinations already encountered.
  seen_ids <- character(0)
  
  # Open the output file connection for writing.
  out_conn <- file(output_file, "w")
  
  # A counter for the total sequences written.
  sequences_written <- 0
  
  # Loop through each file path in the input_files vector.
  for (file in input_files) {
    # Check if the file exists before trying to read it.
    if (!file.exists(file)) {
      warning(paste("File not found:", file, "- Skipping."))
      next
    }
    
    # Read all lines from the current file.
    lines <- readLines(file)
    if (length(lines) == 0) {
      next
    }
    
    # Process the lines of the current file.
    i <- 1
    while (i <= length(lines)) {
      line <- lines[i]
      
      # Check if the line is a FASTA header.
      if (startsWith(line, ">")) {
        header <- line
        
        # Attempt to extract the pigID and week from the header.
        match <- regmatches(header, regexec("_(\\d+)_(\\d+)_", header))
        
        is_duplicate <- FALSE
        
        # Check if the pattern was found.
        if (length(match[[1]]) == 3) {
          # Construct the pigID_week string
          pig_id_week <- paste(match[[1]][2], match[[1]][3], sep = "_")
          
          # Check if this ID has been seen before.
          if (pig_id_week %in% seen_ids) {
            is_duplicate <- TRUE
          } else {
            # If it's a new ID, record it.
            seen_ids <- c(seen_ids, pig_id_week)
          }
        }
        
        # If it's not a duplicate, write the header and the sequence to the output file.
        if (!is_duplicate) {
          writeLines(header, out_conn)
          sequences_written <- sequences_written + 1
          
          # Move to the next line(s), which contain the sequence data.
          i <- i + 1
          while (i <= length(lines) && !startsWith(lines[i], ">")) {
            writeLines(lines[i], out_conn)
            i <- i + 1
          }
        } else {
          # If it is a duplicate, skip the header and all its sequence lines.
          i <- i + 1
          while (i <= length(lines) && !startsWith(lines[i], ">")) {
            i <- i + 1
          }
        }
      } else {
        # This handles cases where a file might not start with a header.
        i <- i + 1
      }
    }
  }
  
  close(out_conn)
  
  cat("Successfully combined", length(input_files), "files into", output_file, "\n")
  cat("A total of", sequences_written, "sequences were written, representing", length(seen_ids), "unique pigID_week combinations.\n")
}

