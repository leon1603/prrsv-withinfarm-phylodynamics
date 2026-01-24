library(knitr)
library(kableExtra)
library(dplyr)

# ==========================================
# USER CONFIGURATION
# ==========================================
output_dir <- "Thesis/Tables" 
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
table_name <- "Reinfection Filtering Results.png" 
full_path <- file.path(output_dir, table_name)

# ==========================================
# RAW DATA DEFINITION
# ==========================================
# We define the raw numbers first so we can calculate totals mathematically
# before converting them to the text format "X / Y".

# --- General Counts ---
B1_All <- c(6, 2, 8, 34, 36)
B1_Inc <- c(6, 0, 8, 25, 14)

B3_All <- c(0, 5, 0, 22, 21)
B3_Inc <- c(0, 5, 0, 18, 5)

# --- BDMM Spatial Counts (All) ---
# Week 2 All=2 assumed to be Farrowing based on exclusion notes
S_Far_All <- c(6, 2, 0, 0, 0) 
S_R1_All  <- c(0, 0, 0, 9, 12)
S_R2_All  <- c(0, 0, 3, 8, 4)
S_R3_All  <- c(0, 0, 5, 17, 11)
S_R4_All  <- c(0, 0, 0, 0, 9)

# --- BDMM Spatial Counts (Included) ---
S_Far_Inc <- c(6, 0, 0, 0, 0)
S_R1_Inc  <- c(0, 0, 0, 9, 5)
S_R2_Inc  <- c(0, 0, 3, 5, 1)
S_R3_Inc  <- c(0, 0, 5, 11, 2)
S_R4_Inc  <- c(0, 0, 0, 0, 6)

# ==========================================
# HELPER FUNCTIONS
# ==========================================

# 1. Combine Function: Creates "All / Inc" string
combine_cols <- function(all, inc) {
  paste0(all, " / ", inc)
}

# 2. Total Function: Sums vectors then combines
calc_total_str <- function(all_vec, inc_vec) {
  t_all <- sum(all_vec)
  t_inc <- sum(inc_vec)
  paste0(t_all, " / ", t_inc)
}

# ==========================================
# DATAFRAME CONSTRUCTION
# ==========================================

# Create the main body text
df_body <- data.frame(
  Week = c("Week 0", "Week 2", "Week 4", "Week 6", "Week 9"),
  
  # General Columns
  Batch1 = combine_cols(B1_All, B1_Inc),
  Batch3 = combine_cols(B3_All, B3_Inc),
  
  # Spacer
  Sp = "",
  
  # Spatial Columns
  Far = combine_cols(S_Far_All, S_Far_Inc),
  R1  = combine_cols(S_R1_All, S_R1_Inc),
  R2  = combine_cols(S_R2_All, S_R2_Inc),
  R3  = combine_cols(S_R3_All, S_R3_Inc),
  R4  = combine_cols(S_R4_All, S_R4_Inc)
)

# Create the Total row
total_row <- c(
  "Total",
  calc_total_str(B1_All, B1_Inc),
  calc_total_str(B3_All, B3_Inc),
  "",
  calc_total_str(S_Far_All, S_Far_Inc),
  calc_total_str(S_R1_All, S_R1_Inc),
  calc_total_str(S_R2_All, S_R2_Inc),
  calc_total_str(S_R3_All, S_R3_Inc),
  calc_total_str(S_R4_All, S_R4_Inc)
)

# Bind them
final_data <- rbind(df_body, total_row)

# ==========================================
# TABLE RENDERING
# ==========================================

final_data %>%
  kbl(
    col.names = c("Week", 
                  "Batch 1", "Batch 3", 
                  " ", # Spacer
                  "Farrowing", "Room 1", "Room 2", "Room 3", "Room 4"),
    caption = "Table 3: Metadata summary of ORF5 sequences before and after filtering for persistent infections. Data are presented as Total Sequences / Included Sequences. 'Included' sequences represent independent transmission events retained for phylodynamic analysis after the exclusion of persistent infections (within-host evolution). The table stratifies counts by batch for the long-term temporal analysis (BDSKY) and by spatial compartment for the fine-scale structured analysis (BDMM).",
    booktabs = TRUE,
    align = c("l", "c", "c", "c", "c", "c", "c", "c", "c") 
  ) %>%
  
  # --- HEADERS ---
  # Grouping the Spatial columns together
  add_header_above(c(" " = 1, 
                     "General Cohorts" = 2, 
                     " " = 1, 
                     "BDMM Spatial Distribution (Batch 1)" = 5), 
                   bold = TRUE) %>%
  
  # --- STYLING ---
  kable_classic(full_width = FALSE, html_font = "Cambria") %>%
  
  # Bold Total Row and add top border
  row_spec(nrow(final_data), bold = TRUE, extra_css = "border-top: 1px solid black;") %>%
  
  # Spacer Column styling (narrower)
  column_spec(4, width = "1em", border_left = FALSE, border_right = FALSE) %>%
  
  # Vertical separation line after General Cohorts
  column_spec(3, border_right = TRUE) %>%
  
  # Light background for Spatial section to visually group it
  column_spec(5:9, background = "#f9f9f9") %>%
  
  # Save
  save_kable(file = full_path, zoom = 3)

message(paste("Combined Slash-separated table saved successfully at:", full_path))