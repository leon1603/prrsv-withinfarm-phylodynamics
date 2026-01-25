#' @title Summarize Multiple BEAST Log Files with Individual Burn-in
#' @description Parses multiple BEAST .log files, removes individual burn-in, and calculates
#' summary statistics (Mean, 95% HPD, ESS) for each parameter.
#' @param log_file_paths A character vector containing full paths to the input BEAST .log files.
#' @param output_dir The directory where the resulting .tsv files should be saved.
#' @param burn_in_fraction A numeric value or vector (0 to 1). If a vector, it must 
#' match the length of log_file_paths. Defaults to 0.1.
#' @author Leon Balthaus, Gemini

summarize_beast_logs <- function(log_file_paths, output_dir, burn_in_fraction = 0.1) {
  
  library(coda)
  
  # 1. Validate Inputs
  n_files <- length(log_file_paths)
  if (n_files == 0) {
    stop("Error: No input log files provided.")
  }
  
  # Handle burn_in_fraction logic
  # If user provides one value, repeat it for all files
  if (length(burn_in_fraction) == 1) {
    burn_in_vector <- rep(burn_in_fraction, n_files)
  } else if (length(burn_in_fraction) == n_files) {
    burn_in_vector <- burn_in_fraction
  } else {
    stop("Error: 'burn_in_fraction' must be a single value or a vector of the same length as 'log_file_paths'.")
  }
  
  if (any(burn_in_vector < 0 | burn_in_vector >= 1)) {
    stop("Error: All 'burn_in_fraction' values must be between 0 and 1.")
  }
  
  # Validate and Create Output Directory
  if (!dir.exists(output_dir)) {
    cat("Warning: Output directory not found. Trying to create it: ", output_dir, "\n")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) {
      stop("Error: Could not create output directory: ", output_dir)
    }
  }
  
  overall_success <- TRUE
  
  # --- 2. MAIN LOOP ---
  cat("Starting batch processing of", n_files, "files...\n")
  cat("--------------------------------------------------\n")
  
  # Using a loop with index 'i' to track the corresponding burn-in
  for (i in seq_along(log_file_paths)) {
    
    log_file_path <- log_file_paths[i]
    current_burn_in <- burn_in_vector[i]
    
    tryCatch({
      # Generate Output Filename
      file_name <- basename(log_file_path)
      
      if (grepl("\\.log$", file_name, ignore.case = TRUE)) {
        new_file_name <- sub("\\.log$", ".tsv", file_name, ignore.case = TRUE)
      } else {
        new_file_name <- paste0(file_name, ".tsv")
      }
      
      output_file_path <- file.path(output_dir, new_file_name)
      
      if (!file.exists(log_file_path)) {
        stop("File not found: ", log_file_path)
      }
      
      cat("\nProcessing:", file_name, "(Burn-in:", current_burn_in * 100, "%)\n")
      
      # --- Read the log file ---
      log_data <- read.delim(log_file_path, comment.char = "#", stringsAsFactors = FALSE)
      
      # --- Apply burn-in ---
      n_states <- nrow(log_data)
      if (n_states == 0) stop("Log file contains no data rows.")
      
      burn_in_states <- floor(current_burn_in * n_states)
      log_data_burned <- log_data[(burn_in_states + 1):n_states, ]
      
      if (nrow(log_data_burned) == 0) {
        stop("No states left after removing burn-in.")
      }
      
      # --- Remove the 'state' or 'Sample' column ---
      if ("state" %in% colnames(log_data_burned)) {
        log_data_burned$state <- NULL
      } else if ("Sample" %in% colnames(log_data_burned)) {
        log_data_burned$Sample <- NULL
      }
      
      # --- Calculate Stats ---
      param_names <- colnames(log_data_burned)
      
      stats_list <- lapply(param_names, function(param_name) {
        suppressWarnings({
          column_data <- log_data_burned[[param_name]]
          mean_val <- mean(column_data, na.rm = TRUE)
          
          hpd_low <- mean_val 
          hpd_up <- mean_val
          ess_val <- NA
          
          data_variance <- var(column_data, na.rm = TRUE)
          
          if (!is.na(data_variance) && data_variance > 0) {
            tryCatch({
              mcmc_obj <- as.mcmc(na.omit(column_data))
              
              if(length(mcmc_obj) > 0) {
                hpd <- HPDinterval(mcmc_obj, prob = 0.95)
                hpd_low <- hpd[1, "lower"] 
                hpd_up <- hpd[1, "upper"]
                ess_val <- effectiveSize(mcmc_obj)
              }
            }, error = function(e_calc) {})
          }
          
          return(data.frame(
            Parameter = param_name,
            Mean = mean_val,
            `95%_HPD_Lower` = hpd_low,
            `95%_HPD_Upper` = hpd_up,
            ESS = ess_val,
            check.names = FALSE,
            stringsAsFactors = FALSE
          ))
        })
      }) 
      
      # --- Assemble and Save ---
      final_table <- do.call(rbind, stats_list)
      row.names(final_table) <- NULL
      write.table(final_table, output_file_path, row.names = FALSE, sep = "\t", quote = FALSE)
      
      cat("-> Saved to:", new_file_name, "\n")
      
    }, error = function(e) {
      cat("!!! Error processing '", basename(log_file_path), "': ", e$message, "\n", sep = "")
      overall_success <<- FALSE
    })
  }
  
  cat("\n--------------------------------------------------\n")
  cat("Batch processing complete.\n")
  return(invisible(overall_success))
}