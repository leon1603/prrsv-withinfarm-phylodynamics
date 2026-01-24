#' @title Classify Infections by Clade
#' @description This function classifies infections based on clade assignment and
#' saves a tree image, a data table, and a FASTA file in specified directories.
#' @param tree_file_path The full path to the input phylogenetic tree file.
#' @param table_output_dir The directory where the classification CSV file will be stored.
#' @param fasta_output_dir The directory where the filtered FASTA file will be stored.
#' @param image_output_dir The directory where the tree plot image will be stored.
#' @param num_clades Integer specifying the number of clades to split the sequences in.
#' @param node_size Size of the circular tips on the tree.
#' @param label_size Size of the sequence text labels on the tree.
#' @param legend_size Size of the legend text and box.
#' @author Leon Balthaus, Gemini

# --- 1. FUNCTION DEFINITION ---

classify_infections_by_clade <- function(tree_file_path, 
                                         table_output_dir, 
                                         fasta_output_dir, 
                                         image_output_dir,
                                         num_clades,
                                         node_size = 1.2,
                                         label_size = 0.7,
                                         legend_size = 1.2) {
  
  library(dplyr)
  library(ape)
  library(Biostrings)
  library(grDevices)
  
  # Validate num_clades argument
  if (missing(num_clades) || !is.numeric(num_clades) || num_clades < 2) {
    stop("Error: 'num_clades' must be a numeric value of 2 or greater.")
  }
  
  # Infer and check file paths
  msa_file_path <- sub("\\.treefile$", "", tree_file_path, ignore.case = TRUE)
  if (!file.exists(tree_file_path)) stop("Error: Tree file not found at: ", tree_file_path)
  if (!file.exists(msa_file_path)) stop("Error: Inferred MSA file not found at: ", msa_file_path)
  
  # --- 2. MAIN LOGIC ---
  tryCatch({
    
    # --- 2a. LOAD DATA AND DEFINE CLADES ---
    cat("--- Loading input files and defining clades... ---\n")
    
    phylo_tree <- read.tree(tree_file_path) 
    all_sequences <- readDNAStringSet(msa_file_path)
    cat(paste(" -> Successfully loaded tree and", length(all_sequences), "sequences.\n"))
    
    patristic_distances <- cophenetic.phylo(phylo_tree)
    dist_matrix <- as.dist(patristic_distances)
    hierarchical_clustering <- hclust(dist_matrix, method = "complete")
    clade_assignments <- cutree(hierarchical_clustering, k = num_clades)
    
    plotting_tree <- as.phylo(hierarchical_clustering)
    
    clade_df <- data.frame(
      SequenceName = names(clade_assignments),
      Clade = clade_assignments
    )
    cat(paste(" -> Successfully assigned all sequences to", num_clades, "clades.\n"))
    
    # Create Metadata and New Display Name
    metadata_df <- data.frame(SequenceName = phylo_tree$tip.label) %>%
      mutate(
        AnimalID = sapply(strsplit(SequenceName, "_"), `[`, 3),
        Week = as.numeric(gsub("w$", "", sapply(strsplit(SequenceName, "_"), `[`, 4))),
        DisplayName = paste0(AnimalID, " (", Week, ")") 
      ) %>%
      left_join(clade_df, by = "SequenceName")
    
    cat(" -> Metadata parsed and merged with clades successfully.\n\n")
    
    # --- 2b. INFECTION CLASSIFICATION LOGIC ---
    cat("--- Classifying infections for each animal based on clades... ---\n")
    unique_pigs <- unique(metadata_df$AnimalID)
    results_list <- list()
    
    for (current_pig in unique_pigs) {
      animal_sequences <- metadata_df %>% filter(AnimalID == current_pig) %>% arrange(Week)
      pig_results <- data.frame()
      if (nrow(animal_sequences) > 0) {
        first_infection <- animal_sequences[1,]
        pig_results <- rbind(pig_results, data.frame(
          pigID_weekOfSampling = paste(first_infection$AnimalID, first_infection$Week, sep = "_"),
          SequenceName = first_infection$SequenceName,
          Clade = first_infection$Clade, Classification = "Initial Infection"
        ))
        if (nrow(animal_sequences) > 1) {
          for (i in 2:nrow(animal_sequences)) {
            current_seq <- animal_sequences[i,]
            previous_seq <- animal_sequences[i - 1,]
            classification <- ifelse(current_seq$Clade == previous_seq$Clade, "Persistent Infection", "Reinfection")
            pig_results <- rbind(pig_results, data.frame(
              pigID_weekOfSampling = paste(current_seq$AnimalID, current_seq$Week, sep = "_"),
              SequenceName = current_seq$SequenceName,
              Clade = current_seq$Clade, Classification = classification
            ))
          }
        }
      }
      results_list[[current_pig]] <- pig_results
    }
    
    final_classification_df <- bind_rows(results_list)
    cat(" -> Classification complete.\n\n")
    
    # --- 2c. PREPARE AND SAVE OUTPUTS ---
    msa_basename <- tools::file_path_sans_ext(basename(msa_file_path))
    
    # Create the requested directories if they don't exist
    dir.create(table_output_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(fasta_output_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(image_output_dir, showWarnings = FALSE, recursive = TRUE)
    
    cat("--- Saving results... ---\n")
    
    # 1. Save CSV Table
    csv_filename <- paste0(msa_basename, "_infection_classification_details.csv")
    csv_path <- file.path(table_output_dir, csv_filename)
    write.csv(final_classification_df, csv_path, row.names = FALSE)
    cat(paste(" -> Classification table saved to:", csv_path, "\n"))
    
    # 2. Save the annotated dendogram
    plot_path <- file.path(image_output_dir, paste0(msa_basename, "_tree.png"))
    
    # A. Prepare colors for branches
    clade_palette <- rainbow(num_clades, s = 1, v = 0.6)
    
    tip_clades_ordered <- clade_df$Clade[match(plotting_tree$tip.label, clade_df$SequenceName)]
    internal_node_descendants <- prop.part(plotting_tree)
    edge_colors <- rep("black", nrow(plotting_tree$edge))
    
    for(i in 1:nrow(plotting_tree$edge)) {
      child_node <- plotting_tree$edge[i, 2]
      if (child_node <= Ntip(plotting_tree)) {
        clade_id <- tip_clades_ordered[child_node]
        edge_colors[i] <- clade_palette[clade_id]
      } else {
        list_index <- child_node - Ntip(plotting_tree)
        descendant_tips_indices <- internal_node_descendants[[list_index]]
        descendant_clades <- tip_clades_ordered[descendant_tips_indices]
        if (length(unique(descendant_clades)) == 1) {
          common_clade_id <- unique(descendant_clades)
          edge_colors[i] <- clade_palette[common_clade_id]
        }
      }
    }
    
    # B. Prepare symbol colors for nodes
    status_levels <- c("Initial Infection", "Persistent Infection", "Reinfection")
    status_palette <- c("blue", "lightgreen", "firebrick")
    names(status_palette) <- status_levels
    status_tip_data <- final_classification_df[match(plotting_tree$tip.label, final_classification_df$SequenceName), ]
    node_symbol_colors <- status_palette[status_tip_data$Classification]
    
    # Rename Tree Tips for Plotting
    current_tip_labels <- plotting_tree$tip.label
    new_tip_labels <- metadata_df$DisplayName[match(current_tip_labels, metadata_df$SequenceName)]
    plotting_tree$tip.label <- new_tip_labels
    
    # C. Create and save the plot
    png(plot_path, width = 11, height = 16, units = "in", res = 150)
    par(mar = c(5, 1, 3, 1)) 
    
    max_x_depth <- max(node.depth.edgelength(plotting_tree))
    plot_xlim <- c(0, max_x_depth * 1.7)
    
    plot(plotting_tree,
         type = "phylogram",
         show.tip.label = FALSE,
         edge.color = edge_colors,
         edge.width = 3,
         x.lim = plot_xlim)
    tiplabels(pch = 21, bg = node_symbol_colors, cex = node_size)
    tiplabels(text = plotting_tree$tip.label,
              frame = "none",
              adj = 0,
              cex = label_size,
              offset = max_x_depth * 0.02)
    
    title(main = "Phylogenetic Tree")
    
    status_counts <- table(final_classification_df$Classification)
    ordered_status_counts <- status_counts[names(status_palette)]
    ordered_status_counts[is.na(ordered_status_counts)] <- 0
    status_legend_labels <- paste0(names(status_palette), " (n=", ordered_status_counts, ")")
    legend_pt_size <- node_size * 1.25
    legend(
      "topleft",
      legend = c(
        as.expression(bquote(bold("Branch Clade"))),
        paste("Clade", 1:num_clades),
        "",
        as.expression(bquote(bold("Infection Status"))),
        status_legend_labels
      ),
      lwd = c(NA, rep(3, num_clades), NA, NA, rep(NA, length(status_palette))),
      pch = c(NA, rep(NA, num_clades), NA, NA, rep(21, length(status_palette))),
      pt.bg = c(NA, rep(NA, num_clades), NA, NA, status_palette),
      col = c(NA, clade_palette, NA, NA, rep("black", length(status_palette))),
      pt.cex = c(NA, rep(NA, num_clades), NA, NA, rep(legend_pt_size, length(status_palette))),
      bty = "n",
      cex = legend_size
    )
    
    dev.off()
    cat(paste(" -> Dual-coded tree image saved to:", plot_path, "\n"))
    
    # 3. Save unique infection sequences as fasta file
    unique_infection_names <- final_classification_df %>%
      filter(Classification %in% c("Initial Infection", "Reinfection")) %>%
      pull(SequenceName)
    unique_infection_sequences <- all_sequences[names(all_sequences) %in% unique_infection_names]
    
    # Suffix appended to fasta file becomes _UniqueInfection.fasta
    fasta_filename <- paste0(msa_basename, "_UniqueInfection.fasta")
    fasta_path <- file.path(fasta_output_dir, fasta_filename)
    writeXStringSet(unique_infection_sequences, fasta_path)
    cat(paste(" -> Unique infection sequences saved to:", fasta_path, "\n"))
    
    cat("\n--- All tasks completed successfully! ---\n")
    return(invisible(TRUE))
    
  }, error = function(e) {
    cat("\n--------------------------------------------------\n")
    cat("An error occurred during script execution:\n")
    message(e)
    cat("--------------------------------------------------\n\n")
    return(invisible(FALSE))
  })
}