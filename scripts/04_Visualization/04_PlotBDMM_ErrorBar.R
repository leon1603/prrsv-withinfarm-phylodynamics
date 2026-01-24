library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(patchwork) # Required for combining plots with custom layouts

#' @title Plot BDMM Results (Model Comparison - Phase Centered)
#' @description Reads multiple BEAST BDMM results files and generates a single comparison plot
#'              where the 'Re' parameter is given full width, and all panels maintain equal height.
#'              "Global" labels on the x-axis are hidden.
#'
#' @param results_file_paths A character vector of full paths to the BDMM results TSV files.
#' @param model_names A character vector of labels corresponding to each file.
#' @param output_dir The base path for the output directory.
#' @param output_filename Filename for the plot.
#' @param include_continuous_models Boolean. If FALSE, excludes any model/file with "Continuous" in its name.
#' @param show_re Boolean. If TRUE, plots Re estimates.
#' @param show_ghost Boolean. If TRUE, plots Ghost Deme Re.
#' @param show_migration Boolean. If TRUE, plots Migration Rates.
#' @param show_sampling Boolean. If TRUE, plots Sampling Proportions.
#' @param show_become_uninfectious Boolean. If TRUE, plots Become Uninfectious Rates.
#' @param show_origin Boolean. If TRUE, plots Origin time.
#' @param show_clock Boolean. If TRUE, plots Clock Rate and CoV.
#' @param show_reservoir Boolean. If TRUE, visualizes estimates associated with the Reservoir phase.
#' @param base_font_size Numeric. Base font size for the plot. Default 12.
#' @param line_width Numeric. Thickness of the error bars. If NULL, defaults to 0.8.
#' 
#' @param axis_text_size Numeric. Font size for X-axis text. Default 16.
#' @param y_axis_text_size Numeric. Font size for Y-axis text. Default 18.
#' @param legend_text_size Numeric. Font size for Legend text. Default 18.
#' @param strip_text_size Numeric. Font size for Facet Strip text. Default 18.
#' @param plot_title_size Numeric. Font size for Plot Title (if added). Default 18.
#' @param point_size Numeric. Size of the estimate points. Default 4.
#' @param figure_height Numeric. Fixed height of the output image in inches. Default 12.
#' @param figure_width Numeric. Width of the output image in inches. Default 14.
#' @param row_height_base Numeric. Base height per row for dynamic calculation (used if figure_height is NULL). Default 3.
#' @param y_axis_padding Numeric. Multiplicative expansion factor for Y-axis limits. Default 0.05.
#' 
#' @return Saves one .png image file.
plot_bdmm_results <- function(results_file_paths, 
                               model_names, 
                               output_dir, 
                               output_filename,
                               include_continuous_models = TRUE,
                               show_re = TRUE,
                               show_ghost = TRUE,
                               show_migration = TRUE, 
                               show_sampling = TRUE,
                               show_become_uninfectious = TRUE, 
                               show_origin = TRUE,
                               show_clock = TRUE,
                               show_reservoir = TRUE,
                               # --- VISUAL SCALING ARGUMENTS ---
                               base_font_size = 12,
                               line_width = NULL,
                               # New Layout Arguments
                               axis_text_size = 16,    
                               y_axis_text_size = 18,
                               legend_text_size = 18,  
                               strip_text_size = 18,   
                               plot_title_size = 18,
                               figure_height = 12,
                               point_size = 4,         
                               figure_width = 14,        
                               row_height_base = 3,    
                               y_axis_padding = 0.05) {
  
  # --- 1. Validation ---
  if(length(results_file_paths) != length(model_names)) {
    stop("Error: The number of file paths must match the number of model names.")
  }
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # --- 1b. Filter Continuous Models (If Requested) ---
  if (!include_continuous_models) {
    continuous_indices <- grep("Continuous", model_names, ignore.case = TRUE)
    
    if (length(continuous_indices) > 0) {
      results_file_paths <- results_file_paths[-continuous_indices]
      model_names <- model_names[-continuous_indices]
      message(paste("Excluded", length(continuous_indices), "Continuous model(s) from the plot."))
    }
    
    if (length(results_file_paths) == 0) {
      stop("Error: All models were filtered out based on the 'include_continuous_models = FALSE' setting.")
    }
  }
  
  # --- 2. Build Target Parameter List ---
  target_params <- c()
  
  if (show_re) {
    target_params <- c(target_params,
                       "ReSPEpi.i1_P1",                  # Farrowing Within-Pen
                       "ReSPEpi.i2_P1",                  # Nursery Within-Pen
                       "ReAmongDemesSPEpi.i2_P1_to_P2",  # Nursery Within-Room
                       "ReAmongDemesSPEpi.i2_P1_to_P5"   # Nursery Between-Room
    )
    if (show_ghost) {
      target_params <- c(target_params, "ReSPEpi.i0_Ghost")
    }
  }
  
  if (show_migration) {
    target_params <- c(target_params,
                       "migrationRateSPEpi.i1",  # Sow -> Farrowing
                       "migrationRateSPEpi.i3",  # Farrowing -> Nursery
                       "migrationRateSPEpi.i4",  # Nursery Continuous
                       "migrationRateSPEpi.i5",  # Nursery Week 6
                       "migrationRateSPEpi.i7"   # Nursery Week 9
    )
  }
  
  if (show_sampling) {
    target_params <- c(target_params,
                       "samplingProportionSPEpi.i1",       
                       "samplingProportionSPEpi.i1_P1",        
                       "samplingProportionSPEpi.i2"             
    )
  }
  
  if (show_become_uninfectious) {
    target_params <- c(target_params,
                       "becomeUninfectiousRateSPEpi.i0", 
                       "becomeUninfectiousRateSPEpi.i1",
                       "becomeUninfectiousRateSPEpi.i2"  
    )
  }
  
  if (show_origin) {
    target_params <- c(target_params, "originBDMMPrime")
  }
  
  if (show_clock) {
    target_params <- c(target_params, "rate.mean", "rate.coefficientOfVariation")
  }
  
  if (length(target_params) == 0) {
    stop("Error: No parameters selected. Set at least one show_* argument to TRUE.")
  }
  
  # --- 3. Data Loading ---
  all_data <- data.frame()
  
  for (i in seq_along(results_file_paths)) {
    file_path <- results_file_paths[i]
    model_name <- model_names[i]
    
    if(!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }
    
    df <- read_tsv(file_path, show_col_types = FALSE)
    available_params <- intersect(target_params, df$Parameter)
    
    if(length(available_params) > 0) {
      df_filtered <- df %>%
        filter(Parameter %in% available_params) %>%
        select(Parameter, Mean, `95%_HPD_Lower`, `95%_HPD_Upper`) %>%
        mutate(Model = model_name)
      
      all_data <- bind_rows(all_data, df_filtered)
    }
  }
  
  if (nrow(all_data) == 0) {
    stop("Error: None of the requested parameters were found in the provided files.")
  }
  
  # --- 3b. Transformation: Rate to Period (DAYS) ---
  all_data <- all_data %>%
    mutate(
      is_uninfectious = grepl("becomeUninfectious", Parameter),
      temp_lower = ifelse(is_uninfectious, 365 / `95%_HPD_Upper`, `95%_HPD_Lower`),
      temp_upper = ifelse(is_uninfectious, 365 / `95%_HPD_Lower`, `95%_HPD_Upper`),
      Mean = ifelse(is_uninfectious, 365 / Mean, Mean)
    ) %>%
    mutate(
      `95%_HPD_Lower` = ifelse(is_uninfectious, temp_lower, `95%_HPD_Lower`),
      `95%_HPD_Upper` = ifelse(is_uninfectious, temp_upper, `95%_HPD_Upper`)
    ) %>%
    select(-is_uninfectious, -temp_lower, -temp_upper)
  
  # --- 4. Categorization & Metadata ---
  all_data <- all_data %>%
    mutate(
      Category = case_when(
        Parameter == "migrationRateSPEpi.i1" ~ "Mig: Sow -> Farrowing",
        Parameter == "migrationRateSPEpi.i3" ~ "Mig: Farrowing -> Nursery",
        Parameter %in% c("migrationRateSPEpi.i4", 
                         "migrationRateSPEpi.i5", 
                         "migrationRateSPEpi.i7") ~ "Mig: Within Nursery",
        
        grepl("^Re", Parameter) ~ "Re",
        
        Parameter %in% c("samplingProportionSPEpi.i1", 
                         "samplingProportionSPEpi.i1_P1") ~ "Sampling Proportion",
        
        Parameter == "samplingProportionSPEpi.i2" ~ "Sampling\n(Nursery)",
        
        grepl("becomeUninfectious", Parameter) ~ "Infectious Period (Days)",
        Parameter == "originBDMMPrime" ~ "Origin Time",
        Parameter == "rate.mean" ~ "Clock Rate",
        Parameter == "rate.coefficientOfVariation" ~ "Clock CoV",
        
        TRUE ~ "Other"
      ),
      
      Phase = case_when(
        Parameter == "migrationRateSPEpi.i1" ~ "Reservoir",        
        Parameter == "migrationRateSPEpi.i3" ~ "Farrowing",        
        Parameter == "migrationRateSPEpi.i4" ~ "Nursery",          
        Parameter == "migrationRateSPEpi.i5" ~ "Nursery (Wk6)", 
        Parameter == "migrationRateSPEpi.i7" ~ "Nursery (Wk9)", 
        
        Parameter == "ReSPEpi.i1_P1" ~ "Farrowing",
        
        grepl("becomeUninfectious.*i0", Parameter) ~ "Reservoir",
        grepl("becomeUninfectious.*i1", Parameter) ~ "Farrowing",
        grepl("becomeUninfectious.*i2", Parameter) ~ "Nursery",
        
        Parameter == "ReSPEpi.i2_P1" ~ "Nursery (Pen)",
        Parameter == "ReAmongDemesSPEpi.i2_P1_to_P2" ~ "Nursery (WR)",
        Parameter == "ReAmongDemesSPEpi.i2_P1_to_P5" ~ "Nursery (BR)",
        
        Parameter == "ReSPEpi.i0_Ghost" ~ "Reservoir",
        
        TRUE ~ "Global" 
      ),
      
      # Assign Color/Shape attributes based on Model Name
      InferenceType = ifelse(grepl("Prior", Model, ignore.case = TRUE), 
                             "Sample from Prior", "Data Informed"),
      
      ModelStructure = ifelse(grepl("Continuous", Model, ignore.case = TRUE), 
                              "Continuous", "Discrete")
    ) %>%
    filter(!(grepl("migrationRate", Parameter) & Mean == 0))
  
  # --- Filter by Reservoir Flag ---
  if (!show_reservoir) {
    all_data <- all_data %>% filter(Phase != "Reservoir")
    
    if (nrow(all_data) == 0) {
      stop("Error: All data was filtered out (everything was in Reservoir phase).")
    }
  }
  
  # --- 5. Factors and Ordering ---
  facet_order <- c(
    "Re",
    "Mig: Sow -> Farrowing",
    "Mig: Farrowing -> Nursery",
    "Mig: Within Nursery",
    "Sampling Proportion",  
    "Sampling\n(Nursery)",
    "Infectious Period (Days)",
    "Origin Time",
    "Clock Rate",
    "Clock CoV"
  )
  
  present_levels <- intersect(facet_order, unique(all_data$Category))
  all_data$Category <- factor(all_data$Category, levels = present_levels)
  # Model factor order still preserved based on input order for dodge grouping
  all_data$Model <- factor(all_data$Model, levels = model_names)
  
  # Define phase levels to ensure logical x-axis order
  phase_levels <- c("Global", "Reservoir", "Farrowing", 
                    "Nursery", "Nursery (Pen)", 
                    "Nursery (WR)", "Nursery (BR)", 
                    "Nursery (Wk6)", "Nursery (Wk9)")
  
  all_data$Phase <- factor(all_data$Phase, levels = intersect(phase_levels, unique(all_data$Phase)))
  
  # --- 6. Set Visual Defaults ---
  
  if (is.null(line_width)) {
    final_line_width <- 0.8
  } else {
    final_line_width <- line_width
  }
  
  # --- 7. Plot Generation ---
  
  # Define common aesthetic mappings and layers
  pd <- position_dodge(width = 0.6)
  
  # Helper list for common layers
  common_layers <- list(
    geom_point(size = point_size, position = pd), # Uses fixed point_size
    geom_errorbar(aes(ymin = `95%_HPD_Lower`, ymax = `95%_HPD_Upper`), 
                  width = 0.2, linewidth = final_line_width, position = pd),
    theme_bw(base_size = base_font_size),
    labs(y = NULL, x = NULL, color = "Inference", shape = "Model Type"),
    
    # Hide "Global" label on X axis
    scale_x_discrete(labels = function(x) ifelse(x == "Global", "", x)),
    
    # Custom scale expansion (padding)
    scale_y_continuous(expand = expansion(mult = y_axis_padding)),
    
    # Updated Theme with custom sizes
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = axis_text_size), 
      axis.text.y = element_text(color = "black", size = y_axis_text_size),                         
      axis.ticks.x = element_line(color = "black"),
      axis.title.y = element_text(size = rel(1.2)),                       
      strip.text = element_text(face = "bold", size = strip_text_size),          
      legend.text = element_text(size = legend_text_size),
      legend.title = element_text(size = rel(1.1), face = "bold"),
      legend.position = "right",
      plot.title = element_text(size = plot_title_size)
    ),
    scale_color_manual(values = c("Data Informed" = "#1b9e77", "Sample from Prior" = "#d95f02")),
    scale_shape_manual(values = c("Discrete" = 16, "Continuous" = 15))
  )
  
  # Split Data into 'Re' and 'Others'
  df_re <- all_data %>% filter(Category == "Re")
  df_other <- all_data %>% filter(Category != "Re")
  
  # -- Generate Re Plot (Full Width) --
  p_re <- NULL
  if(nrow(df_re) > 0) {
    p_re <- ggplot(df_re, aes(x = Phase, y = Mean, color = InferenceType, shape = ModelStructure, group = Model)) +
      common_layers +
      # Add Re-specific line (Re=1)
      geom_hline(yintercept = 1, linetype = "dotted", color = "red", linewidth = 1) +
      facet_wrap(~Category)
  }
  
  # -- Generate Other Plots (Grid) --
  p_other <- NULL
  # Calculate rows needed for 'Other' plot here for layout sizing
  n_rows_other <- 0
  
  if(nrow(df_other) > 0) {
    n_facets_other <- length(unique(df_other$Category))
    n_rows_other <- ceiling(n_facets_other / 3)
    
    p_other <- ggplot(df_other, aes(x = Phase, y = Mean, color = InferenceType, shape = ModelStructure, group = Model)) +
      common_layers +
      # Use scales = "free" so units don't conflict, but layout control will fix heights
      facet_wrap(~Category, scales = "free", ncol = 3)
  }
  
  # -- Combine using Patchwork --
  if(!is.null(p_re) && !is.null(p_other)) {
    # Dynamically assign height ratios: 
    # If Re is 1 row and Others are 3 rows, heights should be c(1, 3) 
    # so that every individual panel has the exact same physical height.
    final_plot <- p_re / p_other + 
      plot_layout(guides = "collect", heights = c(2, n_rows_other))
    
  } else if (!is.null(p_re)) {
    final_plot <- p_re
  } else {
    final_plot <- p_other
  }
  
  # --- 8. Save Output ---
  
  # Height Calculation logic
  n_rows_re <- if(!is.null(p_re)) 1 else 0
  total_visual_rows <- n_rows_re + n_rows_other
  
  if (!is.null(figure_height)) {
    final_height <- figure_height
  } else {
    final_height <- max(4, (total_visual_rows * row_height_base) + 1)
  }
  
  if (!grepl("\\.png$", output_filename, ignore.case = TRUE)) {
    output_filename <- paste0(output_filename, ".png")
  }
  
  full_path <- file.path(output_dir, output_filename)
  
  ggsave(full_path, plot = final_plot, width = figure_width, height = final_height, dpi = 300)
  
  cat(paste("Saved plot:", full_path, "\n"))
}