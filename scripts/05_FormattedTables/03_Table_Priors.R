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
table_name <- "Priors.png" 
full_path <- file.path(output_dir, table_name)

# ==========================================
# DATA PREPARATION
# ==========================================
# We manually set the repeated "Parameter" names to "" (empty string).
# This achieves the visual effect of merged cells but keeps the row structure intact for borders.

data_priors <- data.frame(
  Parameter = c(
    "Effective Reproduction Number", 
    "Infectious Period",
    "",                              # Blank for 2nd row of Infectious Period
    "Sampling Proportion",
    "",                              # Blank for 2nd row of Sampling Proportion
    "Clock Rate (Mean)",
    "",                              # Blank for 2nd row of Clock Rate (Mean)
    "Clock Rate (Std. Dev.)",
    "Tree Origin (tMRCA)",
    "",                              # Blank for 2nd row of Tree Origin
    "Change Times (Epochs)",         
    "Migration Rate",                # NEW PARAMETER
    "Substitution Model"
  ),
  Model = c(
    "Both",
    "BDSKYSA",
    "SMTBD",
    "BDSKYSA",
    "SMTBD",
    "BDSKYSA",
    "SMTBD",
    "Both",
    "BDSKYSA",
    "SMTBD",
    "BDSKYSA",                        
    "SMTBD",                         # Model for Migration Rate (Assumed SMTBD)
    "Both"
  ),
  Prior_Distribution = c(
    "Lognormal(0, 1.0)",
    "Lognormal(8.0, 0.1)",
    "Lognormal(8.0, 0.5)",
    "Beta(1.0, 1.0)",
    "Normal(0.33, 0.09)",
    "Lognormal(0.001, 1.25)",
    "Lognormal(0.0136, 0.1)",
    "Gamma(0.5396, 0.3819)",
    "Lognormal(1.20, 0.05)",
    "Lognormal(0.45, 0.1)",
    "Uniform(0, 1.2)",                
    "Exponential(1.0)",            # Prior for Migration Rate
    "bModelTest"
  ),
  Rationale = c(
    "Uninformative prior centered at 1.",
    "Tight prior (mean ~45 days, 95% CI ~37–65 days) based on experimental data for nursery-age infected pigs which made up the majority of the study population",
    "Broader prior (mean ~45 days, 95% CI ~20–140 days) to evaluate differences between infections in utero and after birth.",
    "Non-informative (uniform) prior allowing free estimation across the range 0 to 1.",
    "Confined to the 95% HPD interval derived from the Batch 1 posterior of the BDSKYSA analysis.",
    "Broad, uninformative distribution capturing the range for observed rates for PRRSV-1 centered around 0.001 substitutions per site per year.",
    "Confined by the posterior estimates from the BDSKYSA analysis due to weak local temporal signal in the Batch 1 only.",
    "Designed to induce an effective prior on the mean clock rate centered around 0.001 substitutions per site per year,",
    "Parameterized relative to the last Batch 3 sample and informed by root-to-tip regression while also accounting for the expected 1-to-8-week subclinical infection phase prior to the first clinical detection.",
    "Parameterized relative to the last Batch 1 sample and informed by root-to-tip regression while also accounting for the expected 1-to-8-week subclinical infection phase prior to the first clinical detection.",
    "Uninformative prior allowing free estimation of epoch configurations between the tMRCA and the end of Batch 3.", 
    "Uninformative prior that accommodates a broad range of biologically plausible rates and ensures robustness in the presence of limited signal for population structure",                             # Rationale for Migration Rate (Left Blank)
    "Co-estimated with default priors for invariable sites, rate heterogeneity, and for base frequencies."
  ),
  Reference = c(
    "-",
    "Charpin et al. (2012)",
    "Charpin et al. (2012), Butler et al., (2010; Rowland et al., (2003); R. W. Wills et al., (2003); You et al., (2022)",
    "-",
    "This study",
    "Boskova et al. (2018); Li et al. (2022), Parisio et al. (2024), Shin et al. (2022)",
    "This study",
    "Boskova et al. (2018)",
    "This study, Pedro Mil-Homens et al. (2024)",
    "This study, Pedro Mil-Homens et al. (2024)",
    "This study",                     
    "Seidel et al. (2024)",                             # Reference for Migration Rate (Left Blank)
    "Bouckaert & Drummond (2017)"
  )
)

# ==========================================
# TABLE GENERATION
# ==========================================
data_priors %>%
  kbl(
    format = "html", 
    col.names = c("Parameter", "Model", "Prior Distribution", "Rationale", "Reference"),
    caption = "Table 1: Prior distributions and rationale for the Bayesian phylodynamic models. Summary of the prior distributions applied to the BDSKYSA and SMTBD models, including reference sources",
    booktabs = TRUE,
    align = c("l", "c", "c", "l", "l") 
  ) %>%
  kable_styling(
    full_width = FALSE,
    position = "center"
  ) %>%
  # 1. HEADER LINE: Thick border (2px)
  row_spec(0, bold = TRUE, extra_css = "border-bottom: 2px solid black;") %>%
  
  # 2. SEPARATOR LINES: Solid border (1px) across the FULL row
  #    Applied to the last row of each parameter group.
  #    Updated indices: Added '12' for the new Migration Rate parameter
  row_spec(c(1, 3, 5, 7, 8, 10, 11, 12), extra_css = "border-bottom: 1px solid black;") %>%
  
  # 3. COLUMN SPECS
  column_spec(4, width = "25em") %>%  
  column_spec(5, width = "12em") %>%  
  
  save_kable(file = full_path, density = 600, zoom = 3)