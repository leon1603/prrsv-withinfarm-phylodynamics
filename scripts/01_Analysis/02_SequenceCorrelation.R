# --- 1. HELPER FUNCTION: CREATE SINGLE MATRIX ---

#' @title Fasta Sequence Percentage Identity Matrix Generator
#' @description This script provides a function to read multiple FASTA
#'              files, calculate pairwise percentage identity matrices, create
#'              a detailed summary table and generate combined heatmaps.
#' @param file_path The path to the input FASTA file.
#' @param num_threads The number of CPU threads to use for calculation.
#' @return A square matrix containing the pairwise percentage identities.
#' @author Leon Balthaus, Gemini

create_identity_matrix_optimized <- function(file_path, num_threads) {
  # Read sequences
  tryCatch({
    sequences <- Biostrings::readDNAStringSet(file_path)
  }, error = function(e) {
    stop("Error reading FASTA file: ", e$message)
  })
  
  if (length(sequences) < 2) {
    warning(paste("Skipping", basename(file_path), "- must contain at least two sequences."))
    return(NULL)
  }
  
  seq_names <- names(sequences)
  
  # Use stringdistmatrix for a fast, parallelized calculation of all pairwise distances
  cat(sprintf(" -> Calculating distance matrix using %d threads...\n", num_threads))
  distance_matrix <- stringdist::stringdistmatrix(sequences, method = "lv", nthread = num_threads)
  
  # Convert the distance matrix to a percentage identity matrix
  seq_lengths <- Biostrings::width(sequences)
  max_len_matrix <- outer(seq_lengths, seq_lengths, pmax)
  distance_matrix <- as.matrix(distance_matrix)
  
  identity_matrix <- (1 - distance_matrix / max_len_matrix) * 100
  diag(identity_matrix) <- 100
  
  dimnames(identity_matrix) <- list(seq_names, seq_names)
  
  return(identity_matrix)
}


# --- 2. MAIN FUNCTION DEFINITION ---

#' @title Analyze Sequence Identity, Generate Summary & Create Heatmaps
#' @description Reads FASTA files, calculates identity matrices, generates a detailed
#' summary table and saves a combined heatmap image.
#' @param fasta_files A character vector of full paths to the input FASTA files.
#' @param figure_dir The directory where the heatmap image will be saved.
#' @param table_dir The directory where the output CSV files will be saved.
#' @param num_threads The number of CPU threads for distance calculation.
#' @return A data frame with the detailed identity summary.
#' @author Gemini, Leon Balthaus

analyze_sequence_identity <- function(fasta_files, figure_dir, table_dir, num_threads = parallel::detectCores() - 1) {
  
  # --- 2a. Load Packages, Validate Inputs ---
  
  for (pkg in c("Biostrings", "stringdist", "parallel", "dplyr", "tidyr", "stringr", "ggplot2", "reshape2", "patchwork", "scales", "igraph")) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package '", pkg, "' is required. Please install it.", sep = ""))
    }
  }
  
  # Create both directories
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
  
  extract_pig_id <- function(header) {
    parts <- unlist(stringr::str_split(header, "_"))
    if (length(parts) >= 4) return(paste(parts[3], parts[4], sep = "_"))
    return(NA)
  }
  
  all_region_summaries <- list()
  all_heatmaps <- list()
  file_stats_summary <- list() 
  
  # --- 2b. MAIN LOGIC ---
  tryCatch({
    cat("--- Starting Sequence Identity Analysis ---\n")
    
    for (file_path in fasta_files) {
      if (!file.exists(file_path)) {
        cat(sprintf(" -> Warning: File not found, skipping: '%s'\n", file_path))
        next
      }
      
      cat(sprintf("\n--- Processing: '%s' ---\n", basename(file_path)))
      
      # 1. Calculate Matrix
      correlation_matrix <- create_identity_matrix_optimized(file_path, num_threads)
      
      if (is.null(correlation_matrix)) next
      
      total_sequences <- nrow(correlation_matrix)
      
      # Dynamically calculate axis text size
      dynamic_size <- 300 / total_sequences 
      axis_label_size <- pmin(6, pmax(4, dynamic_size)) 
      
      min_val <- min(correlation_matrix)
      
      matrix_filename <- paste0(tools::file_path_sans_ext(basename(file_path)), "_identity_matrix.csv")
      
      # Save Matrix CSV to TABLE directory
      write.csv(round(correlation_matrix, 2), file = file.path(table_dir, matrix_filename))
      cat(sprintf(" -> Saved identity matrix to '%s'\n", file.path(table_dir, matrix_filename)))
      region_name <- sub(".*(ORF5|WGS).*", "\\1", basename(file_path))
      
      # --- Heatmap Generation ---
      melted_matrix <- reshape2::melt(correlation_matrix)
      short_labels <- sapply(rownames(correlation_matrix), extract_pig_id)
      
      p <- ggplot2::ggplot(melted_matrix, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradientn(
          name = "Identity (%)",
          limits = c(min_val, 100),
          colours = c("blue", "lightblue", "red"), 
          values = scales::rescale(c(min_val, 99.99, 100)) 
        ) +
        ggplot2::labs(title = region_name, x = "", y = "") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = axis_label_size),
                       axis.text.y = ggplot2::element_text(size = axis_label_size),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::scale_x_discrete(labels = short_labels) +
        ggplot2::scale_y_discrete(labels = short_labels)
      
      all_heatmaps[[region_name]] <- p
      
      # --- Detailed Summary Logic ---
      seq_headers <- rownames(correlation_matrix)
      diag(correlation_matrix) <- 0
      match_indices <- which(correlation_matrix == 100, arr.ind = TRUE)
      
      # --- CLUSTERING & STATS LOGIC START ---
      
      if (nrow(match_indices) == 0) {
        cat(" -> No 100% identical pairs found.\n")
        
        file_stats_summary[[region_name]] <- data.frame(
          Region = region_name,
          Total_Sequences = total_sequences,
          Unique_100_Percent_Pairs = 0,
          Num_Identical_Groups = 0,
          Group_Size_Summary = "NA",
          Num_Unique_Sequences = total_sequences,
          Redundancy_Percentage = 0
        )
        
        next
      }
      
      # Get unique pairs (upper triangle only)
      unique_pairs <- match_indices[match_indices[, "row"] < match_indices[, "col"], , drop = FALSE]
      num_unique_pairs <- nrow(unique_pairs)
      
      # --- Clustering with igraph ---
      
      edges_df <- data.frame(
        from = seq_headers[unique_pairs[, "row"]],
        to = seq_headers[unique_pairs[, "col"]]
      )
      
      g <- igraph::graph_from_data_frame(edges_df, directed = FALSE)
      
      all_seqs_in_graph <- igraph::V(g)$name
      singletons <- setdiff(seq_headers, all_seqs_in_graph)
      g <- igraph::add_vertices(g, nv = length(singletons), attr = list(name = singletons))
      
      components <- igraph::components(g)
      group_sizes <- components$csize
      
      identical_groups <- group_sizes[group_sizes > 1]
      num_identical_groups <- length(identical_groups)
      
      if (num_identical_groups > 0) {
        group_size_table <- table(identical_groups)
        group_size_summary <- paste(paste0(group_size_table, "x", names(group_size_table)), collapse = ", ")
      } else {
        group_size_summary <- "NA"
      }
      
      num_unique_sequences <- components$no
      redundancy_pct <- (total_sequences - num_unique_sequences) / total_sequences * 100
      
      file_stats_summary[[region_name]] <- data.frame(
        Region = region_name,
        Total_Sequences = total_sequences,
        Unique_100_Percent_Pairs = num_unique_pairs,
        Num_Identical_Groups = num_identical_groups,
        Group_Size_Summary = group_size_summary,
        Num_Unique_Sequences = num_unique_sequences,
        Redundancy_Percentage = round(redundancy_pct, 2)
      )
      
      # --- CLUSTERING & STATS LOGIC END ---
      
      if (nrow(unique_pairs) == 0) {
        cat(" -> No unique 100% identical pairs found.\n")
        next
      }
      
      match_df <- data.frame(Header1 = seq_headers[unique_pairs[, "row"]],
                             Header2 = seq_headers[unique_pairs[, "col"]])
      match_df$PigID1 <- sapply(match_df$Header1, extract_pig_id)
      match_df$PigID2 <- sapply(match_df$Header2, extract_pig_id)
      
      matches_long <- dplyr::bind_rows(match_df %>% dplyr::select(ID = PigID1, Matched_To = PigID2),
                                       match_df %>% dplyr::select(ID = PigID2, Matched_To = PigID1))
      
      region_summary <- matches_long %>%
        dplyr::group_by(ID) %>%
        dplyr::summarise(match_count = dplyr::n(),
                         matching_samples = paste(sort(unique(Matched_To)), collapse = ", "),
                         .groups = 'drop') %>%
        dplyr::rename(pigID_weekOfSampling = ID)
      
      cat(sprintf(" -> Found %d samples with 100%% identity matches.\n", nrow(region_summary)))
      names(region_summary) <- c("pigID_weekOfSampling", paste0(region_name, "_match_count"), paste0(region_name, "_matching_samples"))
      all_region_summaries[[region_name]] <- region_summary
    }
    
    # --- 2c. Combine, Finalize, and Save ---
    
    if (length(all_region_summaries) > 0) {
      final_summary <- Reduce(function(x, y) dplyr::full_join(x, y, by = "pigID_weekOfSampling"), all_region_summaries)
      
      final_summary <- final_summary %>%
        dplyr::mutate(dplyr::across(dplyr::ends_with("_match_count"), ~tidyr::replace_na(., 0))) %>%
        dplyr::mutate(dplyr::across(dplyr::ends_with("_matching_samples"), ~tidyr::replace_na(., "")))
      
      # Save Summary CSV to TABLE directory
      summary_output_path <- file.path(table_dir, "identity_summary_by_pigID.csv")
      write.csv(final_summary, file = summary_output_path, row.names = FALSE)
      cat(sprintf("\n--- Saved final summary table to '%s' ---\n", summary_output_path))
    } else {
      cat("\n--- No 100% matches found. No detailed pigID summary table generated. ---\n")
    }
    
    if (length(file_stats_summary) > 0) {
      stats_summary_df <- dplyr::bind_rows(file_stats_summary)
      
      # Save Stats CSV to TABLE directory
      stats_output_path <- file.path(table_dir, "identity_stats_summary.csv")
      write.csv(stats_summary_df, file = stats_output_path, row.names = FALSE)
      cat(sprintf("\n--- Saved identity stats summary to '%s' ---\n", stats_output_path))
    } else {
      cat("\n--- No files processed. No identity stats summary generated. ---\n")
    }
    
    if (length(all_heatmaps) > 0) {
      combined_plot <- patchwork::wrap_plots(all_heatmaps, ncol = 2)
      
      # Save Heatmap Image to figure directory
      heatmap_path <- file.path(figure_dir, "identity_heatmaps.png")
      ggplot2::ggsave(heatmap_path, combined_plot, width = 12, height = 5 * ceiling(length(all_heatmaps)/2), limitsize = FALSE)
      cat(sprintf("--- Saved combined heatmaps to '%s' ---\n", heatmap_path))
    }
    
    cat("--- All identity analysis tasks completed successfully! ---\n")
    return(final_summary)
    
  }, error = function(e) {
    cat("\n--------------------------------------------------\n")
    cat("An error occurred during the sequence identity analysis.\n")
    message(e)
    cat("--------------------------------------------------\n\n")
    return(invisible(NULL))
  })
}