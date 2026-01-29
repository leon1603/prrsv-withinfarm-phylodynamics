#' @title Plot Annotated Rooted Phylogram
#' @description This function loads a phylogenetic tree, automatically midpoint-roots it,
#' and saves a phylogram. It plots bootstrap values, colors branches by Clade, 
#' and annotates tips with Pen (Triangle) and Room (Square) metadata.
#' @param tree_file_path The full path to the input phylogenetic tree file.
#' @param csv_file_path The full path to the input CSV file containing classification metadata.
#' @param output_dir The base path for the output directory.
#' @param output_filename The name of the output PNG file.
#' @param legend_cex Font size scaling factor for the legends.
#' @param tip_cex Font size scaling factor for the sequence names/tips.
#' @param bootstrap_cex Font size scaling factor for the bootstrap values.
#' @param branch_width Numeric value controlling the thickness of the tree branches.
#' @param visualize_clades If TRUE, branches are colored by Clade.
#' @param tip_label_offset Manually sets the space between the tree tip and the sequence name. 
#' @param tree_x_offset Manually shifts the start of the tree to the right.
#' @param node_symbol_cex Scaling factor for the size of the symbols (triangles, squares).
#' @param star_cex Scaling factor for the size of the star (*) symbols for infections.
#' @author Leon Balthaus, Gemini

phylogram <- function(tree_file_path, csv_file_path, output_dir, output_filename,
                      legend_cex = 0.8, tip_cex = 1.0, bootstrap_cex = 0.8,
                      branch_width = 1.5, visualize_clades = TRUE, 
                      tip_label_offset = NULL, tree_x_offset = 0,
                      node_symbol_cex = 1.2, star_cex = 1.5) { 
  
  library(ape)
  library(grDevices)
  library(phytools)
  
  # Check file paths
  if (!file.exists(tree_file_path)) stop("Error: Tree file not found at: ", tree_file_path)
  if (!file.exists(csv_file_path)) stop("Error: CSV file not found at: ", csv_file_path)
  
  tryCatch({
    cat("--- Generating phylogram with Pen & Room structure... ---\n")
    
    # 1. Load data
    phylo_tree <- read.tree(tree_file_path)
    meta_data <- read.csv(csv_file_path, stringsAsFactors = FALSE)
    
    # 1b. Automatic midpoint rooting
    if (!is.rooted(phylo_tree)) {
      cat(" -> Tree is unrooted. Performing midpoint rooting...\n")
      phylo_tree <- midpoint.root(phylo_tree)
    } else {
      cat(" -> Tree is already rooted.\n")
    }
    
    plotting_tree <- ladderize(phylo_tree)
    
    # Store original labels for parsing metadata
    original_labels <- plotting_tree$tip.label
    num_tips <- length(original_labels)
    
    # Create lookup maps
    meta_map <- setNames(meta_data$Classification, meta_data$pigID_weekOfSampling)
    
    # Check if Clade column exists
    has_clade <- "Clade" %in% colnames(meta_data)
    if (has_clade) {
      clade_map <- setNames(meta_data$Clade, meta_data$pigID_weekOfSampling)
    } else {
      if(visualize_clades) warning("Column 'Clade' not found in CSV. Branches will not be colored by clade.")
      clade_map <- NULL
    }
    
    star_chars <- rep(NA, num_tips)
    star_colors <- rep(NA, num_tips)
    tip_clade_ids <- rep(NA, num_tips) 
    
    extract_raw_key <- function(x) {
      parts <- strsplit(x, "_")[[1]]
      if (length(parts) >= 6) {
        return(paste0(parts[5], "_", parts[6]))
      }
      return(NA)
    }
    
    for (i in 1:num_tips) {
      key <- extract_raw_key(original_labels[i])
      if (!is.na(key)) {
        if (key %in% names(meta_map)) {
          classification <- meta_map[[key]]
          
          # Infection classification logic
          if (classification == "Initial Infection") {
            star_chars[i] <- "*"
            star_colors[i] <- "black"
          } else if (classification == "Reinfection") {
            star_chars[i] <- "*"
            star_colors[i] <- "red"
          } 
        }
        
        if (has_clade && key %in% names(clade_map)) {
          tip_clade_ids[i] <- clade_map[[key]]
        }
      }
    }
    
    # Parsing functions
    extract_pen <- function(x) {
      parts <- strsplit(x, "_")[[1]]
      if (length(parts) >= 2 && grepl("^P", parts[2])) return(parts[2]) else return("Other")
    }
    extract_room <- function(x) {
      parts <- strsplit(x, "_")[[1]]
      if (length(parts) >= 3 && grepl("^R", parts[3])) return(parts[3]) else return("Other")
    }
    
    extract_week <- function(x) {
      parts <- strsplit(x, "_")[[1]]
      if (length(parts) >= 6) {
        val <- parts[6]
        if (val == "3") return("2")
        return(val)
      } else return("Unknown")
    }
    
    extract_pig_id <- function(x) {
      parts <- strsplit(x, "_")[[1]]
      if (length(parts) >= 5) return(parts[5]) else return(x)
    }
    
    # Data processing
    numeric_sort <- function(ids) {
      if (length(ids) == 0) return(ids)
      nums <- suppressWarnings(as.numeric(gsub("[^0-9]", "", ids)))
      return(ids[order(nums, ids)])
    }
    
    pen_ids <- sapply(original_labels, extract_pen)
    valid_pen_indices <- which(pen_ids != "PUnknown" & pen_ids != "Other")
    unique_pens <- numeric_sort(unique(pen_ids[valid_pen_indices]))
    
    room_ids <- sapply(original_labels, extract_room)
    valid_room_indices <- which(room_ids != "RUnknown" & room_ids != "Other")
    unique_rooms <- numeric_sort(unique(room_ids[valid_room_indices]))
    
    farrowing_idx <- which(grepl("Farrowing", unique_rooms, ignore.case = TRUE))
    if (length(farrowing_idx) > 0) {
      unique_rooms <- c(unique_rooms[farrowing_idx], unique_rooms[-farrowing_idx])
    }
    
    week_ids <- sapply(original_labels, extract_week)
    unique_weeks <- unique(week_ids)
    num_weeks <- suppressWarnings(as.numeric(unique_weeks))
    if (!any(is.na(num_weeks))) unique_weeks <- unique_weeks[order(num_weeks)] else unique_weeks <- sort(unique_weeks)
    
    # Color generation
    distinct_color_generator <- colorRampPalette(c("red", "blue", "green", "darkorange", "purple", "cyan", "magenta", "brown", "darkgreen", "navy", "pink"))
    if(length(unique_pens) > 0) {
      pen_palette <- distinct_color_generator(length(unique_pens))
      names(pen_palette) <- unique_pens
      tip_pen_colors <- pen_palette[pen_ids] 
    }
    
    fallback_colors <- c("blue", "green", "gray", "purple", "cyan", "orange", "magenta")
    if(length(unique_rooms) > 0) {
      room_palette <- character(length(unique_rooms))
      names(room_palette) <- unique_rooms
      fallback_counter <- 1
      for(i in 1:length(unique_rooms)) {
        r <- unique_rooms[i]
        clean_name <- sub("^R", "", r) 
        if (grepl("Farrowing", clean_name, ignore.case = TRUE)) room_palette[i] <- "yellow"
        else if (clean_name == "Nursery1") room_palette[i] <- "brown"
        else if (clean_name == "Nursery2") room_palette[i] <- "red"
        else if (clean_name == "Nursery3") room_palette[i] <- "#FF7F7F" 
        else if (clean_name == "Nursery4") room_palette[i] <- "pink"
        else {
          room_palette[i] <- fallback_colors[((fallback_counter - 1) %% length(fallback_colors)) + 1]
          fallback_counter <- fallback_counter + 1
        }
      }
      tip_room_colors <- room_palette[room_ids]
    }
    
    # Clade colors
    edge_colors <- rep("black", nrow(plotting_tree$edge)) 
    clade_palette <- NULL
    
    if (has_clade && visualize_clades) {
      unique_clades <- sort(unique(na.omit(tip_clade_ids)))
      n_clades <- length(unique_clades)
      if (n_clades > 0) {
        raw_clade_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628")
        if (n_clades <= length(raw_clade_colors)) clade_palette <- raw_clade_colors[1:n_clades]
        else clade_palette <- colorRampPalette(raw_clade_colors)(n_clades)
        names(clade_palette) <- unique_clades
        
        get_tips_for_node <- function(node, tree) {
          if (node <= length(tree$tip.label)) return(node)
          children <- tree$edge[tree$edge[,1] == node, 2]
          return(unlist(lapply(children, get_tips_for_node, tree = tree)))
        }
        
        for (i in 1:nrow(plotting_tree$edge)) {
          child_node <- plotting_tree$edge[i, 2]
          descendant_tips <- get_tips_for_node(child_node, plotting_tree)
          descendant_clades <- tip_clade_ids[descendant_tips]
          descendant_clades <- descendant_clades[!is.na(descendant_clades)]
          if (length(descendant_clades) > 0) {
            first_clade <- descendant_clades[1]
            if (all(descendant_clades == first_clade)) edge_colors[i] <- clade_palette[as.character(first_clade)]
          }
        }
      }
    }
    
    # Update tip labels
    final_pig_ids <- sapply(original_labels, extract_pig_id)
    final_week_nums <- sapply(original_labels, extract_week)
    final_labels <- paste0(final_pig_ids, " (", final_week_nums, ")")
    plotting_tree$tip.label <- final_labels
    
    # 2. Prepare output
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    plot_path <- file.path(output_dir, output_filename)
    png(plot_path, width = 16, height = 20, units = "in", res = 150)
    
    max_x_depth <- max(node.depth.edgelength(plotting_tree))
    
    left_sidebar_width <- (max_x_depth * 0.3) + tree_x_offset 
    
    if (is.null(tip_label_offset)) {
      text_offset <- max_x_depth * 0.01 
    } else {
      text_offset <- tip_label_offset
    }
    
    plot_width <- (max_x_depth + text_offset) * 1.5
    right_limit <- max(max_x_depth * 1.8, plot_width)
    plot_xlim <- c(-left_sidebar_width, right_limit)
    
    plot(plotting_tree,
         type = "phylogram",    
         show.tip.label = TRUE,
         align.tip.label = TRUE,
         label.offset = text_offset, 
         edge.color = edge_colors, 
         edge.width = branch_width, 
         cex = tip_cex,                
         x.lim = plot_xlim, 
         tip.color = "black", 
         main = "") 
    
    # Add stars for infection classification
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    label_widths <- strwidth(final_labels, cex = tip_cex)
    
    text_start_x <- max(lastPP$xx) + text_offset
    star_x_coords <- text_start_x + label_widths + (max_x_depth * 0.02)
    
    for (i in 1:num_tips) {
      if (!is.na(star_chars[i])) {
        text(x = star_x_coords[i], y = lastPP$yy[i], 
             labels = star_chars[i], col = star_colors[i], 
             cex = star_cex, font = 2) 
      }
    }
    
    # --- ADD SYMBOLS ---
    
    symbol_width_user <- strwidth("M", cex = node_symbol_cex)
    symbol_gap <- symbol_width_user * 0.7
    
    # 1. Triangle (Pen)
    offset_pen <- symbol_gap
    
    if (length(valid_pen_indices) > 0) {
      tiplabels(tip = valid_pen_indices, pch = 24, bg = tip_pen_colors[valid_pen_indices], 
                col = "black", cex = node_symbol_cex, offset = offset_pen) 
    }
    
    # 2. Square (Room)
    offset_room <- symbol_gap * 2.2
    
    if (length(valid_room_indices) > 0) {
      tiplabels(tip = valid_room_indices, pch = 22, bg = tip_room_colors[valid_room_indices], 
                col = "black", cex = node_symbol_cex, offset = offset_room) 
    }
    
    # Legends
    legend_x <- -left_sidebar_width * 0.95 
    current_y <- Ntip(plotting_tree)
    legend_gap <- Ntip(plotting_tree) * 0.03
    
    if (!is.null(clade_palette)) {
      l_clade <- legend(x = legend_x, y = current_y,
                        legend = paste("Clade", names(clade_palette)), col = clade_palette,
                        lty = 1, lwd = branch_width, cex = legend_cex, 
                        bty = "n", title = "Clade (Branches)", title.adj = 0)
      current_y <- current_y - l_clade$rect$h - legend_gap
    }
    
    if (length(unique_pens) > 0) {
      l2 <- legend(x = legend_x, y = current_y,
                   legend = sub("^P", "", unique_pens), pt.bg = pen_palette,
                   pch = 24, col = "black", 
                   pt.cex = node_symbol_cex, cex = legend_cex, 
                   bty = "n", title = "Pen Location", title.adj = 0)
      current_y <- current_y - l2$rect$h - legend_gap
    }
    
    if (length(unique_rooms) > 0) {
      room_lbls <- sub("Nursery([0-9])", "Nursery \\1", sub("^R", "", unique_rooms))
      l3 <- legend(x = legend_x, y = current_y,
                   legend = room_lbls, pt.bg = room_palette,
                   pch = 22, col = "black", 
                   pt.cex = node_symbol_cex, cex = legend_cex, 
                   bty = "n", title = "Room Location", title.adj = 0)
      current_y <- current_y - l3$rect$h - legend_gap
    }
    
    # Bootstrap visualization
    bs_values <- plotting_tree$node.label
    if (!is.null(bs_values)) {
      # 1. Identify valid nodes and labels
      valid_indices <- which((suppressWarnings(as.numeric(bs_values)) >= 71) | 
                               (bs_values != "" & bs_values != "NA" & is.na(suppressWarnings(as.numeric(bs_values)))))
      
      if (length(valid_indices) > 0) {
        node_ids <- Ntip(plotting_tree) + valid_indices
        labels   <- bs_values[valid_indices]
        
        # 2. Get coordinates
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        x_coords <- lastPP$xx[node_ids]
        y_coords <- lastPP$yy[node_ids]
        
        # 3. Calculate text dimensions
        char_width  <- strwidth(labels, cex = bootstrap_cex, font = 2)
        char_height <- strheight(labels, cex = bootstrap_cex, font = 2)
        
        padding_x <- strwidth("M", cex = bootstrap_cex, font = 2) * 0.2
        padding_y <- strheight("M", cex = bootstrap_cex, font = 2) * 0.2
        
        box_w <- char_width + 2 * padding_x
        box_h <- char_height + 2 * padding_y
        
        # 4. Initial positions
        gap_from_node <- strwidth("M", cex = bootstrap_cex) * 0.3
        
        final_x_right <- x_coords - gap_from_node
        final_x_left  <- final_x_right - box_w
        final_y_top   <- y_coords + (box_h / 2)
        final_y_bot   <- y_coords - (box_h / 2)
        
        # 5. Collision Detection & Shifting
        boxes <- data.frame(
          id = 1:length(node_ids),
          left = final_x_left,
          right = final_x_right,
          top = final_y_top,
          bottom = final_y_bot,
          label = labels,
          stringsAsFactors = FALSE
        )
        
        # Iterative pass to resolve overlaps
        for (iter in 1:5) { 
          hit <- FALSE
          for (i in 1:nrow(boxes)) {
            for (j in 1:nrow(boxes)) {
              if (i == j) next
              
              # Check Y overlap
              y_overlap <- (boxes$top[i] > boxes$bottom[j]) && (boxes$bottom[i] < boxes$top[j])
              
              if (y_overlap) {
                # Check X overlap
                x_overlap <- (boxes$right[i] > boxes$left[j]) && (boxes$left[i] < boxes$right[j])
                
                if (x_overlap) {
                  hit <- TRUE
                  # Resolve: Move the leftmost one further left
                  center_i <- (boxes$left[i] + boxes$right[i])/2
                  center_j <- (boxes$left[j] + boxes$right[j])/2
                  
                  shift <- (min(boxes$right[i], boxes$right[j]) - max(boxes$left[i], boxes$left[j])) + gap_from_node
                  
                  if (center_i <= center_j) {
                    boxes$left[i]  <- boxes$left[i] - shift
                    boxes$right[i] <- boxes$right[i] - shift
                  } else {
                    boxes$left[j]  <- boxes$left[j] - shift
                    boxes$right[j] <- boxes$right[j] - shift
                  }
                }
              }
            }
          }
          if (!hit) break 
        }
        
        # 6. Draw squares for bootsrap values
        rect(boxes$left, boxes$bottom, boxes$right, boxes$top, 
             col = "white", border = "black")
        
        text_x <- (boxes$left + boxes$right) / 2
        text_y <- (boxes$bottom + boxes$top) / 2
        
        text(text_x, text_y, labels = boxes$label, 
             cex = bootstrap_cex, font = 2, col = "black", adj = c(0.5, 0.5))
      }
    }
    
    add.scale.bar(x = legend_x, y = 0, cex = legend_cex)
    dev.off()
    
    cat(paste(" -> Tree image saved to:", plot_path, "\n"))
    cat("--- Task complete! ---\n")
    return(invisible(TRUE))
    
  }, error = function(e) {
    cat("\nError during plotting:\n")
    message(e)
    return(invisible(FALSE))
  })
}