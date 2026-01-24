#' @title Plot BDSKY Model Results
#' @description Generates a comparison plot for all inferred posterior estimates across WGS, ORF5 and Sample From Prior analyses
#' @param figure_height Total height of the output image in inches. If NULL, height is calculated automatically.
#' @param figure_width Total width of the output image in inches.
#' @param axis_text_size Size of the X-axis labels.
#' @param y_axis_text_size Size of the Y-axis numbers.
#' @param legend_text_size Size of the Legend text and title
#' @param strip_text_size Size of the Facet Titles.
#' @param plot_title_size Size of the Main Plot Title
plot_bdsky_results <- function(results_file_paths, 
                               color_groups, 
                               data_groups, 
                               output_dir, 
                               output_filename,
                               plot_title = "Comparison of BDSKY Model Results",
                               
                               # --- DIMENSIONS ---
                               figure_width = 14,       
                               figure_height = NULL,   
                               row_height_base = 3,    
                               
                               # --- FONT SIZES ---
                               axis_text_size = 12,    
                               y_axis_text_size = 12,
                               legend_text_size = 12,  
                               strip_text_size = 14,   
                               plot_title_size = 16,
                               
                               point_size = 4,         
                               y_axis_padding = 0.05,
                               show_re_values = TRUE,
                               show_change_times = TRUE,
                               show_become_uninfectious = FALSE,
                               show_clock_rate = FALSE,
                               show_sampling_batch1 = FALSE,
                               show_sampling_batch2 = FALSE,
                               show_origin = FALSE,
                               show_cov = FALSE,
                               exclude_high_re_epoch2 = FALSE) {
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  
  # --- 1. Validation ---
  if(length(results_file_paths) != length(color_groups) || length(results_file_paths) != length(data_groups)) {
    stop("Error: The number of file paths must match the length of 'color_groups' and 'data_groups'.")
  }
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # --- 2. Data Processing Loop ---
  plot_data <- data.frame()
  
  for (i in seq_along(results_file_paths)) {
    file_path <- results_file_paths[i]
    color_label <- color_groups[i] 
    data_label  <- data_groups[i] 
    
    if(!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }
    
    # Read TSV
    df <- read_tsv(file_path, show_col_types = FALSE)
    
    # --- Filtering ---
    keep_rows <- rep(FALSE, nrow(df))
    
    if (show_re_values) keep_rows <- keep_rows | grepl("^ReSPEpi\\.i\\d+$", df$Parameter)
    if (show_change_times) keep_rows <- keep_rows | grepl("^ReSPEpi\\.i\\d+_endtime$", df$Parameter)
    if (show_become_uninfectious) keep_rows <- keep_rows | df$Parameter == "becomeUninfectiousRateSPEpi"
    if (show_clock_rate) keep_rows <- keep_rows | df$Parameter == "rate.mean"
    if (show_cov) keep_rows <- keep_rows | df$Parameter == "rate.coefficientOfVariation"
    if (show_origin) keep_rows <- keep_rows | df$Parameter == "originBDMMPrime"
    if (show_sampling_batch1) keep_rows <- keep_rows | df$Parameter == "samplingProportionSPEpi.i1"
    if (show_sampling_batch2) keep_rows <- keep_rows | df$Parameter == "samplingProportionSPEpi.i3"
    
    # Subset Data
    df_subset <- df[keep_rows, ] %>%
      select(Parameter, Mean, `95%_HPD_Lower`, `95%_HPD_Upper`) %>%
      mutate(
        ColorLabel = color_label, 
        DataSource = data_label   
      )
    
    # --- Transformation: Rate to Period to days ---
    df_subset <- df_subset %>%
      mutate(
        is_uninfectious = Parameter == "becomeUninfectiousRateSPEpi",
        temp_lower = ifelse(is_uninfectious, 365 / `95%_HPD_Upper`, `95%_HPD_Lower`),
        temp_upper = ifelse(is_uninfectious, 365 / `95%_HPD_Lower`, `95%_HPD_Upper`),
        Mean = ifelse(is_uninfectious, 365 / Mean, Mean)
      ) %>%
      mutate(
        `95%_HPD_Lower` = ifelse(is_uninfectious, temp_lower, `95%_HPD_Lower`),
        `95%_HPD_Upper` = ifelse(is_uninfectious, temp_upper, `95%_HPD_Upper`)
      ) %>%
      select(-is_uninfectious, -temp_lower, -temp_upper)
    
    # --- Labelling ---
    df_subset <- df_subset %>%
      mutate(
        Category = case_when(
          grepl("endtime", Parameter) ~ "Change Times (Years)",
          grepl("^ReSPEpi", Parameter) & !grepl("sampling", Parameter) ~ "Re Estimates",
          Parameter == "becomeUninfectiousRateSPEpi" ~ "Infectious Period (Days)",
          Parameter == "rate.mean" ~ "Molecular Clock - Rate",
          Parameter == "rate.coefficientOfVariation" ~ "Molecular Clock - CoV",
          
          # Separate Batches
          Parameter == "samplingProportionSPEpi.i1" ~ "Sampling Proportion - Batch 1",
          Parameter == "samplingProportionSPEpi.i3" ~ "Sampling Proportion - Batch 2",
          Parameter == "originBDMMPrime" ~ "Origin (Years From Present)",
          TRUE ~ "Other"
        ),
        
        Raw_Index = as.numeric(str_extract(Parameter, "\\d+")),
        Display_Index = Raw_Index + 1,
        
        Label = case_when(
          Category == "Re Estimates" & Raw_Index == 0 ~ "Epoch 1 (Recent)",
          Category == "Re Estimates" ~ paste0("Epoch ", Display_Index),
          Category == "Change Times (Years)" ~ paste0("End of Epoch ", Display_Index),
          
          # Blank labels for single-parameter plots
          Parameter == "becomeUninfectiousRateSPEpi" ~ "",
          Parameter == "rate.mean" ~ "",
          Parameter == "rate.coefficientOfVariation" ~ "",
          Parameter == "originBDMMPrime" ~ "",
          Parameter == "samplingProportionSPEpi.i1" ~ "", 
          Parameter == "samplingProportionSPEpi.i3" ~ "", 
          TRUE ~ Parameter
        )
      )
    
    plot_data <- rbind(plot_data, df_subset)
  }
  
  if (nrow(plot_data) == 0) {
    stop("Error: No matching parameters found.")
  }
  
  # --- 3. FILTERING & CLEANING ---
  if (exclude_high_re_epoch2) {
    plot_data <- plot_data %>%
      filter(!(Label == "Epoch 2" & ColorLabel %in% c("4 Epochs", "5 Epochs")))
  }
  
  # Handle NA labels
  plot_data$Label[is.na(plot_data$Label) | plot_data$Label == "NA"] <- ""
  
  # --- 4. CREATE  X-AXIS ---
  
  # A. Define Parameter Order
  unique_base_labels <- unique(plot_data$Label)
  re_labels   <- sort(unique_base_labels[grepl("^Epoch", unique_base_labels)], decreasing = TRUE)
  time_labels <- sort(unique_base_labels[grepl("^End of Epoch", unique_base_labels)], decreasing = TRUE)
  fixed_end   <- c("Origin", "") 
  param_order <- unique(c(re_labels, time_labels, fixed_end, unique_base_labels))
  
  # B. Define Data Source Order
  available_ds <- unique(plot_data$DataSource)
  preferred_order <- c("ORF5", "WGS", "Prior", "SampleFromPrior") 
  
  sorted_ds <- intersect(preferred_order, available_ds)
  remaining_ds <- setdiff(available_ds, preferred_order)
  ds_order <- c(sorted_ds, remaining_ds)
  
  plot_data$DataSource <- factor(plot_data$DataSource, levels = ds_order)
  
  # C. Generate Interleaved Combined Levels
  combined_levels <- c()
  
  for (lbl in param_order) {
    for (ds in ds_order) {
      if (lbl == "" || lbl == " ") {
        val <- as.character(ds)
      } else {
        val <- paste0(lbl, "\n", as.character(ds))
      }
      combined_levels <- c(combined_levels, val)
    }
  }
  
  # D. Map Data to New Labels
  plot_data <- plot_data %>%
    mutate(
      CombinedLabel = ifelse(Label == "" | Label == " ", 
                             as.character(DataSource), 
                             paste0(Label, "\n", as.character(DataSource)))
    )
  
  # E. Final Ordering
  final_levels_ordered <- combined_levels
  
  # --- 5. Final Factors ---
  plot_data$CombinedLabel <- factor(plot_data$CombinedLabel, levels = final_levels_ordered)
  plot_data$ColorLabel    <- factor(plot_data$ColorLabel, levels = unique(color_groups)) 
  
  # Updated cat_levels to match the new string
  cat_levels <- c("Re Estimates", "Change Times (Years)", "Origin (Years From Present)", 
                  "Infectious Period (Days)", "Molecular Clock - Rate", "Molecular Clock - CoV", 
                  "Sampling Proportion - Batch 1", "Sampling Proportion - Batch 2")
  
  plot_data$Category <- factor(plot_data$Category, levels = intersect(cat_levels, unique(plot_data$Category)))
  
  # --- 6. Plotting ---
  p <- ggplot(plot_data, aes(x = CombinedLabel, y = Mean, color = ColorLabel)) + 
    
    {if(any(plot_data$Category == "Re Estimates")) 
      geom_hline(data = filter(plot_data, Category == "Re Estimates"), 
                 aes(yintercept = 1), linetype = "dotted", color = "red", linewidth = 1)} +
    
    geom_errorbar(aes(ymin = `95%_HPD_Lower`, ymax = `95%_HPD_Upper`), 
                  width = 0.2, linewidth = 0.8, 
                  position = position_dodge(width = 0.6)) +
    
    geom_point(shape = 16, size = point_size, position = position_dodge(width = 0.6)) +
    
    facet_wrap(~Category, scales = "free", ncol = 3) + 
    scale_y_continuous(expand = expansion(mult = y_axis_padding)) +
    scale_x_discrete(labels = function(x) ifelse(x == "NA" | is.na(x), "", x)) +
    
    theme_bw(base_size = 12) + 
    labs(
      title = plot_title, 
      y = NULL,
      x = NULL,
      color = "Model Configuration"
    ) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = axis_text_size),
      axis.text.y = element_text(size = y_axis_text_size),
      legend.text = element_text(size = legend_text_size),
      legend.title = element_text(size = legend_text_size, face = "bold"),
      strip.text = element_text(size = strip_text_size, face = "bold"), 
      plot.title = element_text(size = plot_title_size, face = "bold", hjust = 0.5),
      
      legend.position = "bottom",
      legend.box = "vertical",
      panel.spacing = unit(1, "lines")
    ) +
    guides(color = guide_legend(ncol = 3))
  
  # --- 7. Saving ---
  full_output_path <- file.path(output_dir, output_filename)
  
  n_cats <- length(unique(plot_data$Category))
  n_rows <- ceiling(n_cats / 3)
  
  final_height <- if(!is.null(figure_height)) {
    figure_height
  } else {
    max(3, n_rows * row_height_base)
  }
  
  ggsave(full_output_path, plot = p, width = figure_width, height = final_height, dpi = 300)
  
  message(paste("Plot saved to:", full_output_path))
}