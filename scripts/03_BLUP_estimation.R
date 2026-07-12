############################################################
##
## S3011 Perennial Sorghum GWAS Pipeline
##
## Script 03: BLUP Estimation
##
## Purpose:
##   - Fit a common linear mixed model for each trait
##   - Estimate accession-level best linear unbiased predictors
##     (BLUPs)
##   - Estimate broad-sense heritability from REML variance
##     components
##   - Generate the accession-level phenotype panel used in
##     downstream genome-wide analyses
##
## Model:
##
##   y = mu + Plot + Accession + residual
##
##   Implemented as:
##
##   y ~ 1 + (1 | Plot) + (1 | Accession)
##
##   Plot and Accession are fitted as random effects.
##   Models are fitted by restricted maximum likelihood (REML).
##
## Required preceding scripts:
##   00_setup.R
##   01_genotype_QC.R
##   02_aggregate_accession_genotypes.R
##
## Required inputs:
##   results/intermediate/geno_mat_imp_accessionLevel.rds
##   results/intermediate/pheno_raw.rds
##
## Outputs:
##   results/intermediate/blup_panel_all_traits.rds
##   results/intermediate/blup_models_all_traits.rds
##   results/tables/BLUP_panel_all_traits.csv
##   results/tables/Table_Heritability_BLUPmodel.csv
##   results/tables/VarianceComponents_BLUPmodel.csv
##
############################################################

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(tibble)
})

############################################################
## 1. Verify required pipeline objects
############################################################

required_objects <- c(
  "out_intermediate",
  "out_tables"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0L) {
  stop(
    "Required pipeline objects are missing: ",
    paste(missing_objects, collapse = ", "),
    ". Run 00_setup.R before running this script."
  )
}

############################################################
## 2. Load accession-level genotype and phenotype data
############################################################

geno_file <- file.path(
  out_intermediate,
  "geno_mat_imp_accessionLevel.rds"
)

pheno_file <- file.path(
  out_intermediate,
  "pheno_raw.rds"
)

if (!file.exists(geno_file)) {
  stop(
    "Accession-level genotype file not found: ",
    geno_file
  )
}

if (!file.exists(pheno_file)) {
  stop(
    "Phenotype file not found: ",
    pheno_file
  )
}

geno_mat_imp <- readRDS(geno_file)
pheno <- readRDS(pheno_file)

############################################################
## 3. Define grouping variables
############################################################

plot_col <- "Plot"
accession_col <- "Accession"

required_pheno_columns <- c(
  plot_col,
  accession_col
)

missing_pheno_columns <- setdiff(
  required_pheno_columns,
  colnames(pheno)
)

if (length(missing_pheno_columns) > 0L) {
  stop(
    "Required phenotype columns are missing: ",
    paste(
      missing_pheno_columns,
      collapse = ", "
    )
  )
}

pheno[[plot_col]] <- factor(
  pheno[[plot_col]]
)

pheno[[accession_col]] <- factor(
  pheno[[accession_col]]
)

############################################################
## 4. Define traits
############################################################

## Trait column names must match the phenotype data exactly.

traits_spec <- list(
  
  list(
    column = "Rhiz_num",
    type = "numeric",
    label = "Rhiz"
  ),
  
  list(
    column = "Regrowth",
    type = "binary",
    label = "Regrowth"
  ),
  
  list(
    column = "Plant.Height",
    type = "numeric",
    label = "PlantHeight"
  ),
  
  list(
    column = "Lateral.tillers",
    type = "numeric",
    label = "Tillers_LateralBranches"
  )
  
)

## Add HeadCompactness only when the column is present.
if ("Compact" %in% colnames(pheno)) {
  
  traits_spec <- append(
    traits_spec,
    list(
      list(
        column = "Compact",
        type = "numeric",
        label = "HeadCompactness"
      )
    )
  )
  
}

############################################################
## 5. Helper functions
############################################################

## Recode common binary values to 0/1.

recode_binary01 <- function(x) {
  
  x_clean <- trimws(
    tolower(
      as.character(x)
    )
  )
  
  dplyr::case_when(
    
    x_clean %in% c(
      "yes",
      "y",
      "1",
      "true",
      "present"
    ) ~ 1,
    
    x_clean %in% c(
      "no",
      "n",
      "0",
      "false",
      "absent"
    ) ~ 0,
    
    TRUE ~ NA_real_
    
  )
  
}

## Convert numeric traits safely.

sanitize_numeric <- function(x, trait_name) {
  
  x_character <- trimws(
    as.character(x)
  )
  
  x_character[
    x_character %in% c(
      "",
      "NA",
      "NaN",
      ".",
      "-"
    )
  ] <- NA_character_
  
  x_numeric <- suppressWarnings(
    as.numeric(x_character)
  )
  
  invalid_values <- unique(
    x_character[
      !is.na(x_character) &
        is.na(x_numeric)
    ]
  )
  
  if (length(invalid_values) > 0L) {
    warning(
      "Trait '",
      trait_name,
      "' contained non-numeric values that were converted to NA: ",
      paste(
        invalid_values,
        collapse = ", "
      )
    )
  }
  
  x_numeric
  
}

## Extract one variance component safely.

extract_variance <- function(
    variance_table,
    grouping_factor
) {
  
  variance_value <- variance_table$vcov[
    variance_table$grp == grouping_factor
  ]
  
  if (length(variance_value) == 0L) {
    return(0)
  }
  
  as.numeric(
    variance_value[1]
  )
  
}

############################################################
## 6. Fit one trait
############################################################

fit_trait <- function(
    data,
    trait_column,
    trait_type,
    trait_label,
    plot_column = "Plot",
    accession_column = "Accession"
) {
  
  if (!trait_column %in% colnames(data)) {
    stop(
      "Trait column not found in phenotype data: ",
      trait_column
    )
  }
  
  analysis_data <- data
  
  response_column <- paste0(
    "Y__",
    trait_label
  )
  
  if (trait_type == "binary") {
    
    analysis_data[[response_column]] <-
      recode_binary01(
        analysis_data[[trait_column]]
      )
    
  } else {
    
    analysis_data[[response_column]] <-
      sanitize_numeric(
        analysis_data[[trait_column]],
        trait_name = trait_column
      )
    
  }
  
  ## Retain complete records for the response and grouping factors.
  
  analysis_data <- analysis_data[
    !is.na(analysis_data[[response_column]]) &
      !is.na(analysis_data[[plot_column]]) &
      !is.na(analysis_data[[accession_column]]),
    ,
    drop = FALSE
  ]
  
  analysis_data[[plot_column]] <-
    droplevels(
      factor(
        analysis_data[[plot_column]]
      )
    )
  
  analysis_data[[accession_column]] <-
    droplevels(
      factor(
        analysis_data[[accession_column]]
      )
    )
  
  if (nrow(analysis_data) < 10L) {
    stop(
      "Too few non-missing observations for trait: ",
      trait_column
    )
  }
  
  if (
    nlevels(
      analysis_data[[accession_column]]
    ) < 2L
  ) {
    stop(
      "Fewer than two accessions remain for trait: ",
      trait_column
    )
  }
  
  ##########################################################
  ## Mixed model
  ##########################################################
  
  ## Plot identifiers are unique across the experiment,
  ## with each plot containing a single accession planting.
  ## Plot is therefore fitted as the field-level random
  ## effect, while Accession estimates accession-level
  ## genetic differences.
  
  model_formula <- as.formula(
    paste0(
      response_column,
      " ~ 1 + (1 | ",
      plot_column,
      ") + (1 | ",
      accession_column,
      ")"
    )
  )
  
  fitted_model <- lmer(
    model_formula,
    data = analysis_data,
    REML = TRUE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(
        maxfun = 200000
      )
    )
  )
  
  ##########################################################
  ## Variance components and heritability
  ##########################################################
  
  variance_table <- as.data.frame(
    VarCorr(
      fitted_model
    )
  )
  
  variance_accession <- extract_variance(
    variance_table,
    accession_column
  )
  
  variance_plot <- extract_variance(
    variance_table,
    plot_column
  )
  
  variance_residual <- sigma(
    fitted_model
  )^2
  
  total_variance <-
    variance_accession +
    variance_plot +
    variance_residual
  
  broad_sense_h2 <- if (
    is.finite(total_variance) &&
    total_variance > 0
  ) {
    
    variance_accession /
      total_variance
    
  } else {
    
    NA_real_
    
  }
  
  ##########################################################
  ## Accession-level BLUPs
  ##########################################################
  
  accession_blups <- ranef(
    fitted_model
  )[[accession_column]]
  
  accession_blups <- data.frame(
    
    Accession = rownames(
      accession_blups
    ),
    
    BLUP = as.numeric(
      accession_blups[[1]]
    ),
    
    stringsAsFactors = FALSE
    
  )
  
  colnames(accession_blups)[2] <- paste0(
    trait_label,
    "_BLUP"
  )
  
  ##########################################################
  ## Model diagnostics
  ##########################################################
  
  singular_model <- isSingular(
    fitted_model,
    tol = 1e-4
  )
  
  convergence_messages <-
    fitted_model@optinfo$conv$lme4$messages
  
  if (is.null(convergence_messages)) {
    convergence_messages <- NA_character_
  } else {
    convergence_messages <- paste(
      convergence_messages,
      collapse = "; "
    )
  }
  
  list(
    
    model = fitted_model,
    
    blups = accession_blups,
    
    variance_components = tibble(
      
      Trait = trait_label,
      
      Variance_Accession =
        variance_accession,
      
      Variance_Plot =
        variance_plot,
      
      Variance_Residual =
        variance_residual,
      
      H2 =
        broad_sense_h2,
      
      N_Observations =
        nrow(analysis_data),
      
      N_Accessions =
        nlevels(
          analysis_data[[accession_column]]
        ),
      
      N_Plots =
        nlevels(
          analysis_data[[plot_column]]
        ),
      
      Singular_Model =
        singular_model,
      
      Convergence_Message =
        convergence_messages
      
    )
    
  )
  
}

############################################################
## 7. Fit all traits
############################################################

blup_list <- list()
variance_component_list <- list()
model_list <- list()

for (trait_specification in traits_spec) {
  
  trait_result <- fit_trait(
    
    data = pheno,
    
    trait_column =
      trait_specification$column,
    
    trait_type =
      trait_specification$type,
    
    trait_label =
      trait_specification$label,
    
    plot_column =
      plot_col,
    
    accession_column =
      accession_col
    
  )
  
  blup_list[[
    trait_specification$label
  ]] <- trait_result$blups
  
  variance_component_list[[
    trait_specification$label
  ]] <- trait_result$variance_components
  
  model_list[[
    trait_specification$label
  ]] <- trait_result$model
  
  message(
    
    "Trait completed: ",
    trait_specification$label,
    
    " | H2 = ",
    round(
      trait_result$variance_components$H2,
      4
    ),
    
    " | observations = ",
    trait_result$variance_components$N_Observations,
    
    " | accessions = ",
    trait_result$variance_components$N_Accessions,
    
    " | plots = ",
    trait_result$variance_components$N_Plots
    
  )
  
}

############################################################
## 8. Merge trait-specific BLUPs
############################################################

## Use full joins so an accession is not discarded merely
## because one trait is missing.

blup_panel <- Reduce(
  
  function(x, y) {
    
    full_join(
      x,
      y,
      by = "Accession"
    )
    
  },
  
  blup_list
  
)

############################################################
## 9. Restrict the BLUP panel to genotyped accessions
############################################################

genotyped_accessions <- rownames(
  geno_mat_imp
)

blup_panel <- blup_panel %>%
  
  filter(
    Accession %in%
      genotyped_accessions
  ) %>%
  
  arrange(
    match(
      Accession,
      genotyped_accessions
    )
  )

if (
  anyDuplicated(
    blup_panel$Accession
  ) > 0L
) {
  stop(
    "Duplicated accessions were detected in the final BLUP panel."
  )
}

############################################################
## 10. Assemble variance-component and heritability tables
############################################################

variance_components_table <- bind_rows(
  variance_component_list
)

heritability_table <- variance_components_table %>%
  
  select(
    Trait,
    H2,
    N_Observations,
    N_Accessions,
    N_Plots,
    Singular_Model,
    Convergence_Message
  )

############################################################
## 11. Save outputs
############################################################

saveRDS(
  
  blup_panel,
  
  file.path(
    out_intermediate,
    "blup_panel_all_traits.rds"
  )
  
)

saveRDS(
  
  model_list,
  
  file.path(
    out_intermediate,
    "blup_models_all_traits.rds"
  )
  
)

write.csv(
  
  blup_panel,
  
  file.path(
    out_tables,
    "BLUP_panel_all_traits.csv"
  ),
  
  row.names = FALSE
  
)

write.csv(
  
  heritability_table,
  
  file.path(
    out_tables,
    "Table_Heritability_BLUPmodel.csv"
  ),
  
  row.names = FALSE
  
)

write.csv(
  
  variance_components_table,
  
  file.path(
    out_tables,
    "VarianceComponents_BLUPmodel.csv"
  ),
  
  row.names = FALSE
  
)

############################################################
## 12. Completion messages
############################################################

message(
  "03_BLUP_estimation.R completed successfully."
)

message(
  "BLUP panel saved to: ",
  file.path(
    out_intermediate,
    "blup_panel_all_traits.rds"
  )
)

message(
  "Readable BLUP table saved to: ",
  file.path(
    out_tables,
    "BLUP_panel_all_traits.csv"
  )
)

message(
  "Heritability table saved to: ",
  file.path(
    out_tables,
    "Table_Heritability_BLUPmodel.csv"
  )
)

message(
  "Variance-component table saved to: ",
  file.path(
    out_tables,
    "VarianceComponents_BLUPmodel.csv"
  )
)