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
table_name <- "validation_summary.png" 
full_path <- file.path(output_dir, table_name)

# ==========================================
# DATA PREPARATION
# ==========================================
data_val <- data.frame(
  Description = c(
    "Cohort WGS Consensus",       # Row 1 (Total)
    "Batch 1",                    # Row 2 (Subset)
    "Batch 3",                    # Row 3 (Subset) -> Line after this
    "Cohort ORF5 Consensus",      # Row 4 (Total)
    "Batch 1",                    # Row 5 (Subset)
    "Batch 3",                    # Row 6 (Subset) -> Line after this
    "Regional Sequences",         # Row 7 -> Line after this
    "Licenced Vaccines in Spain", # Row 8 -> Line after this
    "WGS NCBIVirus (taxID 1965066)", # Row 9 -> Line after this
    "ORF5 NCBIVirus (taxID 1965066)" # Row 10 (End of table)
  ),
  Num_Sequences = c(
    "35", "24", "11",             
    "134", "86", "48",            
    "10", "5", "37", "197"
  ),
  Min_Length = c(
    "15098", "", "",              
    "606", "", "", 
    "14443", "14758", "14932", "606"
  ),
  Max_Length = c(
    "15098", "", "", 
    "606", "", "", 
    "14910", "15120", "15428", "606"
  ),
  Ambiguities_Total = c(
    "2", "", "", 
    "0", "", "", 
    "132", "0", "225", "97"
  ),
  Ambiguity_Chars = c(
    "Y", "", "", 
    "None", "", "", 
    "K, M, N, R, S, W, Y", "None", "K, N, R, S, Y", "K, M, R, S, W, Y"
  )
)

# ==========================================
# TABLE GENERATION
# ==========================================
data_val %>%
  kbl(
    col.names = c("Dataset", "Sequences", "Min Length", "Max Length", "Total Ambiguities", "Ambiguity Characters"),
    caption = "Table 2: Summary of PRRSV-1 sequence datasets, quality metrics, and batch partitioning. Overview of the Cohort and Extended datasets for Whole Genome Sequences (WGS) and ORF5, detailing sequence counts, length ranges, and ambiguity codes.",
    booktabs = TRUE,
    align = c("l", "c", "c", "c", "c", "l")
  ) %>%
  kable_classic(full_width = FALSE, html_font = "Cambria") %>%
  
  # --- HIERARCHY ---
  add_indent(c(2, 3, 5, 6)) %>%
  
  # --- SEPARATORS ---
  # Add horizontal lines after each complete dataset block
  # Row 3 (End of WGS block), Row 6 (End of ORF5 block), 
  # Row 7 (Regional), Row 8 (Vaccines), Row 9 (WGS NCBI)
  row_spec(c(3, 6, 7, 8, 9), extra_css = "border-bottom: 1px solid #d3d3d3;") %>%
  
  # --- STYLING ---
  # Header style
  row_spec(0, bold = TRUE, extra_css = "border-bottom: 2px solid black;") %>%
  
  # Color the subsets slightly lighter for visual hierarchy
  row_spec(c(2, 3, 5, 6), color = "#404040") %>%
  
  # Bottom border for the whole table (End of Row 10)
  row_spec(nrow(data_val), extra_css = "border-bottom: 2px solid black;") %>%
  
  # Column widths
  column_spec(1, width = "15em") %>%
  column_spec(6, width = "12em") %>%
  
  save_kable(file = full_path, density = 600, zoom = 3)

message(paste("Validation table saved at:", full_path))