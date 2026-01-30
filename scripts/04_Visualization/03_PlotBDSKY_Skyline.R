#' @title Plot BEAST BDSKY Skylines
#' @description This function reads multiple BEAST results and FASTA files 
#' and generates a single PDF containing all skyline plots arranged in a grid.
#' @param input_groups A list of lists. Each sub-list must contain:
#'    - `results`: Vector of file paths to BEAST .tsv files.
#'    - `fasta`: Single file path to the corresponding FASTA file.
#' @param output_dir The base path for the output directory.
#' @param output_filename The specific name for the output PDF file.
#' @param t_present_date The calendar date of the "present" (most recent sample).
#' @param density_bandwidth Controls smoothness of density curve.
#' @param show_batch_bars If TRUE, visualizes batch start/end dates as colored vertical bars.
#' @param outbreak_date (Optional) "YYYY-MM-DD" string for outbreak start.
#' @param batch1_start_date (Optional) Start of Batch 1.
#' @param batch1_end_date (Optional) End of Batch 1.
#' @param batch2_start_date (Optional) Start of Batch 2. 
#' @param batch2_end_date (Optional) End of Batch 2. 
#' @param batch3_start_date (Optional) Start of Batch 3.
#' @param batch3_end_date (Optional) End of Batch 3.
#' @param title_font_size Font size for the plot titles (A, B, C...)
#' @param axis_font_size Font size for axis tick labels (dates/numbers)
#' @param legend_font_size Font size for the legend text.
#' @author Leon Balthaus, Gemini

Re_plot <- function(input_groups, 
                    output_dir, 
                    output_filename, 
                    t_present_date,
                    density_bandwidth = 0.2, 
                    show_batch_bars = TRUE, 
                    outbreak_date = NULL,
                    batch1_start_date = NULL,
                    batch1_end_date = NULL,
                    batch2_start_date = NULL, 
                    batch2_end_date = NULL, 
                    batch3_start_date = NULL,
                    batch3_end_date = NULL,
                    title_font_size = 24,
                    axis_font_size = 10,
                    legend_font_size = 8) {
  
  # --- 1a. Load Required Libraries ---
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(lubridate) 
  library(stringr) 
  library(gridExtra)
  library(RColorBrewer) 
  
  # --- 1b. Validate Shared Arguments between ORF5 and WGS plots ---
  validate_batch_dates <- function(start_str, end_str, batch_name) {
    if (is.null(start_str) && is.null(end_str)) return(NULL) 
    if (is.null(start_str) || is.null(end_str)) stop(paste("Error:", batch_name, "requires both start and end."))
    start_obj <- ymd(start_str)
    end_obj <- ymd(end_str)
    if (is.na(start_obj) || is.na(end_obj)) stop(paste("Error:", batch_name, "dates invalid."))
    if (start_obj > end_obj) stop(paste("Error:", batch_name, "start must be before end."))
    return(list(start = start_obj, end = end_obj))
  }
  
  if (missing(t_present_date)) stop("Error: 't_present_date' is required.")
  present_date_obj <- ymd(t_present_date)
  outbreak_date_obj <- if (!is.null(outbreak_date)) ymd(outbreak_date) else NULL
  batch1_dates <- validate_batch_dates(batch1_start_date, batch1_end_date, "Batch 1")
  batch2_dates <- validate_batch_dates(batch2_start_date, batch2_end_date, "Batch 2")
  batch3_dates <- validate_batch_dates(batch3_start_date, batch3_end_date, "Batch 3")
  
  if (missing(input_groups) || length(input_groups) == 0) stop("Error: 'input_groups' must be provided.")
  
  # --- 1c. Process FASTA Dates ---
  get_fasta_dates <- function(f_path) {
    if (is.null(f_path)) return(NULL)
    cat(paste("--- Processing FASTA for dates:", basename(f_path), "---\n"))
    tryCatch({
      fasta_lines <- readLines(f_path)
      headers <- fasta_lines[grep("^>", fasta_lines)]
      date_strings <- str_extract(headers, "[0-9]{4}-[0-9]{2}-[0-9]{2}$")
      sampling_dates <- ymd(date_strings)
      sampling_dates <- sampling_dates[!is.na(sampling_dates)]
      
      if(length(sampling_dates) > 0) {
        cat(paste(" -> Found", length(sampling_dates), "dates.\n"))
        return(data.frame(date = sampling_dates))
      } else {
        warning(paste("No dates found in", basename(f_path)))
        return(NULL)
      }
    }, error = function(e) {
      warning(paste("Failed to read FASTA:", e$message))
      return(NULL)
    })
  }
  
  # --- 2. MAIN LOGIC ---
  all_plots <- list()
  
  for (group in input_groups) {
    current_results_files <- group$results
    current_fasta_path <- group$fasta
    sampling_dates_df <- get_fasta_dates(current_fasta_path)
    
    for (current_file_path in current_results_files) {
      tryCatch({
        cat(paste("--- Processing Log:", basename(current_file_path), "---\n"))
        summary_data <- read_tsv(current_file_path, show_col_types = FALSE)
        
        # --- Parse Parameters ---
        # 1. Re Values
        re_params <- summary_data %>% 
          filter(grepl("^ReSPEpi\\.i[0-9]+$", Parameter)) %>%
          mutate(index = as.integer(gsub("ReSPEpi\\.i", "", Parameter))) %>%
          arrange(index)
        
        # 2. Change Times (End Times)
        time_params <- summary_data %>%
          filter(grepl("^ReSPEpi\\.i[0-9]+_endtime$", Parameter)) %>%
          mutate(index = as.integer(gsub("ReSPEpi\\.i|_endtime", "", Parameter))) %>%
          arrange(index)
        
        # 3. Origin
        origin_param <- summary_data %>% filter(Parameter == "originBDMMPrime")
        
        t_present_numeric   <- origin_param$Mean
        t_present_lower_hpd <- origin_param$`95%_HPD_Lower`
        t_present_upper_hpd <- origin_param$`95%_HPD_Upper`
        
        # --- Calculate Dates for Skyline ---
        date_origin_mean <- as.Date(present_date_obj - duration(t_present_numeric, "years"))
        date_origin_earliest <- as.Date(present_date_obj - duration(t_present_upper_hpd, "years"))
        date_origin_latest <- as.Date(present_date_obj - duration(t_present_lower_hpd, "years"))
        
        dates_end <- as.Date(date_origin_mean + duration(time_params$Mean, "years"))
        plot_start_dates <- c(date_origin_mean, dates_end) 
        plot_end_dates <- c(dates_end, present_date_obj)
        
        ribbon_data <- data.frame(xmin = plot_start_dates, xmax = plot_end_dates, ymin = re_params$`95%_HPD_Lower`, ymax = re_params$`95%_HPD_Upper`)
        line_data <- data.frame(x_time = c(plot_start_dates, present_date_obj), y_mean = c(re_params$Mean, re_params$Mean[nrow(re_params)]))
        
        
        # 1. Prepare tMRCA (Origin) Data
        ct_hpd_data <- data.frame(
          x_min = date_origin_earliest,
          x_max = date_origin_latest,
          x_mean = date_origin_mean,
          Transition = "tMRCA"
        )
        
        # 2. Prepare Epoch Data and combine
        if(nrow(time_params) > 0) {
          ct_lower_dates <- as.Date(date_origin_mean + duration(time_params$`95%_HPD_Lower`, "years"))
          ct_upper_dates <- as.Date(date_origin_mean + duration(time_params$`95%_HPD_Upper`, "years"))
          ct_mean_dates  <- as.Date(date_origin_mean + duration(time_params$Mean, "years"))
          
          epoch_labels <- paste("End of Epoch", time_params$index + 1)
          
          epoch_data <- data.frame(
            x_min = ct_lower_dates,
            x_max = ct_upper_dates,
            x_mean = ct_mean_dates,
            Transition = epoch_labels
          )
          
          # Combine tMRCA and Epochs into one dataframe
          ct_hpd_data <- bind_rows(ct_hpd_data, epoch_data)
        }
        
        # 3. Sort by x_min to optimize packing order
        ct_hpd_data <- ct_hpd_data %>% arrange(x_min)
        
        # Checks collisions for ALL bars (tMRCA and Epochs alike)
        y_step <- 0.15 
        y_base <- 0.1 
        
        ct_hpd_data$y_pos <- 0
        levels_occupied <- list() 
        
        for (i in 1:nrow(ct_hpd_data)) {
          current_xmin <- ct_hpd_data$x_min[i]
          current_xmax <- ct_hpd_data$x_max[i]
          
          collision <- TRUE
          level_idx <- 0
          assigned_level <- 0
          
          while(collision) {
            collision <- FALSE
            if (length(levels_occupied) >= (level_idx + 1)) {
              existing_ranges <- levels_occupied[[level_idx + 1]]
              if (!is.null(existing_ranges)) {
                for (r in existing_ranges) {
                  if (current_xmin <= r$end && current_xmax >= r$start) {
                    collision <- TRUE
                    break
                  }
                }
              }
            }
            if (collision) {
              level_idx <- level_idx + 1 
            } else {
              assigned_level <- level_idx 
            }
          }
          
          ct_hpd_data$y_pos[i] <- y_base + (assigned_level * y_step)
          
          if (length(levels_occupied) < (assigned_level + 1)) {
            levels_occupied[[assigned_level + 1]] <- list()
          }
          levels_occupied[[assigned_level + 1]] <- append(levels_occupied[[assigned_level + 1]], list(list(start=current_xmin, end=current_xmax)))
        }
        
        # 5. Set Factor Levels for Legend (Ensure tMRCA is always first in legend)
        if(nrow(time_params) > 0) {
          epoch_levels <- paste("End of Epoch", sort(unique(time_params$index + 1)))
          ct_hpd_data$Transition <- factor(ct_hpd_data$Transition, levels = c("tMRCA", epoch_levels))
        } else {
          ct_hpd_data$Transition <- factor(ct_hpd_data$Transition, levels = "tMRCA")
        }
        
        # --- Plot Setup ---
        batch1_breaks <- if(show_batch_bars && !is.null(batch1_dates)) c(batch1_dates$start, batch1_dates$end) else NULL
        batch2_breaks <- if(show_batch_bars && !is.null(batch2_dates)) c(batch2_dates$start, batch2_dates$end) else NULL
        batch3_breaks <- if(show_batch_bars && !is.null(batch3_dates)) c(batch3_dates$start, batch3_dates$end) else NULL
        
        # --- CLEAN AXIS BREAKS ---
        axis_breaks <- sort(unique(c(plot_start_dates, present_date_obj, 
                                     if(!is.null(outbreak_date_obj)) outbreak_date_obj, 
                                     batch1_breaks, batch2_breaks, batch3_breaks)))
        
        plot_min_date <- min(c(axis_breaks, date_origin_earliest)) 
        plot_max_date <- max(axis_breaks) 
        plot_max_extended <- plot_max_date + ceiling(as.numeric(plot_max_date - plot_min_date) * 0.05)
        
        skyline_plot <- ggplot() +
          # 1. Background Shading (Batches)
          (if (show_batch_bars && !is.null(batch1_dates)) annotate("rect", xmin = batch1_dates$start, xmax = batch1_dates$end, ymin = -Inf, ymax = Inf, fill = "darkseagreen", alpha = 0.4) else NULL) +
          (if (show_batch_bars && !is.null(batch2_dates)) annotate("rect", xmin = batch2_dates$start, xmax = batch2_dates$end, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.4) else NULL) +
          (if (show_batch_bars && !is.null(batch3_dates)) annotate("rect", xmin = batch3_dates$start, xmax = batch3_dates$end, ymin = -Inf, ymax = Inf, fill = "lightcoral", alpha = 0.4) else NULL) +
          
          # 2. Sampling Density
          (if (!is.null(sampling_dates_df)) geom_density(data = sampling_dates_df, aes(x = date, y = after_stat(scaled) * (2)), adjust = density_bandwidth, fill = "purple", color = NA, alpha = 0.8) else NULL) +
          
          # 3. Main Skyline (Re)
          geom_rect(data = ribbon_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "steelblue", alpha = 0.3) +
          geom_step(data = line_data, aes(x = x_time, y = y_mean), color = "black", linewidth = 0.8, direction = "hv") +
          geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick", linewidth = 1) +
          
          # 4. Outbreak Date
          (if (!is.null(outbreak_date_obj)) geom_vline(xintercept = outbreak_date_obj, linetype = "dotted", color = "darkorange", linewidth = 1) else NULL) +
          
          # 5. All Transition HPD Intervals (tMRCA + Epochs)
          (if (!is.null(ct_hpd_data)) 
            geom_segment(data = ct_hpd_data, 
                         aes(x = x_min, xend = x_max, y = y_pos, yend = y_pos, color = Transition), 
                         linewidth = 2, lineend = "round", alpha = 0.8)
           else NULL) +
          
          (if (!is.null(ct_hpd_data))
            geom_point(data = ct_hpd_data, aes(x = x_mean, y = y_pos), color="white", size=1)
           else NULL) +
          
          # 6. Styling
          labs(title = NULL, y = expression(R[e]), x = NULL, color = NULL) +
          guides(color = guide_legend(ncol = 3)) + 
          
          scale_x_date(breaks = axis_breaks, date_labels = "%Y-%m-%d", guide = guide_axis(check.overlap = TRUE)) +
          scale_color_brewer(palette = "Dark2") + 
          theme_bw() +
          theme(
            axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5, size = axis_font_size), 
            axis.text.y = element_text(size = axis_font_size),
            axis.title.y = element_text(size = axis_font_size + 2),
            
            plot.title = element_text(size = title_font_size, face = "bold"), 
            
            panel.grid.minor = element_blank(),
            
            legend.position = "bottom",
            legend.box.margin = margin(t = -5, b = 5),
            legend.text = element_text(size = legend_font_size),
            legend.title = element_blank() 
          ) +
          coord_cartesian(ylim = c(0, 6), xlim = c(plot_min_date, plot_max_extended), expand = FALSE)
        
        all_plots[[length(all_plots) + 1]] <- skyline_plot
      }, error = function(e) { cat(paste("\nError processing:", basename(current_file_path), "\n", e$message, "\n")) })
    } 
  } 
  
  # --- 3. SAVE TO PDF ---
  if (length(all_plots) > 0) {
    
    # --- 3b. Overwrite Titles with A, B, C... ---
    cat("--- Assigning panel labels (A, B, C...) ---\n")
    for (i in seq_along(all_plots)) {
      new_title <- if(i <= 26) LETTERS[i] else paste("Panel", i)
      all_plots[[i]] <- all_plots[[i]] + labs(title = new_title)
    }
    
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    pdf_path <- file.path(output_dir, output_filename)
    n_plots <- length(all_plots)
    n_cols <- 2
    n_rows <- ceiling(n_plots / n_cols)
    dynamic_height <- n_rows * 5.0 
    dynamic_width <- 14 
    pdf(pdf_path, width = dynamic_width, height = dynamic_height)
    cat(paste("\n--- arranging", n_plots, "plots into a", n_rows, "x", n_cols, "grid ---\n"))
    grid.arrange(grobs = all_plots, ncol = n_cols)
    dev.off()
    cat(paste("--- All plots saved to:", pdf_path, "---\n"))
  } else {
    warning("No valid plots were generated. PDF not created.")
  }
}