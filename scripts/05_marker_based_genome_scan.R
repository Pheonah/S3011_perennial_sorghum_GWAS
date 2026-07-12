############################################################
##
## S3011 Perennial Sorghum GWAS Pipeline
##
## Script 05: Complementary Marker-Based Genome Scan
##
## Purpose:
##   - Conduct complementary single-marker regression scans
##     using accession-level mean allele dosage and accession-
##     level BLUPs
##   - Generate complete marker-scan summary statistics
##   - Identify the strongest marker on each chromosome
##   - Calculate Benjamini-Hochberg FDR and genomic inflation
##   - Generate Manhattan and Q-Q plots for all analyzed traits
##
## Analytical framework:
##
##   BLUP = beta0 + betaSNP(SNP dosage) + error
##
## This analysis is a marker-level genome scan and is not
## classical linkage-based interval mapping. It is intended to
## evaluate whether candidate regions prioritized by the
## exploratory PCA-adjusted GWAS are similarly prioritized by
## an alternative marker-level analytical framework.
##
## Required preceding scripts:
##   00_setup.R
##   01_genotype_QC.R
##   02_aggregate_accession_genotypes.R
##   03_BLUP_estimation.R
##   04_GWAS.R
##
## Required inputs:
##   results/intermediate/geno_mat_imp_accessionLevel.rds
##   results/intermediate/marker_metadata_accessionLevel.rds
##   results/intermediate/blup_panel_all_traits.rds
##
## Outputs:
##   results/intermediate/marker_scan_all_traits.rds
##   results/tables/marker_scan_all_results.csv
##   results/tables/marker_scan_peaks_per_chromosome.csv
##   results/tables/marker_scan_top20_per_trait.csv
##   results/tables/marker_scan_run_summary.csv
##   results/tables/marker_scan_QQ_statistics.csv
##   results/figures/marker_scan_<trait>_Manhattan.png
##   results/figures/marker_scan_<trait>_QQ.png
##
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

############################################################
## 1. Verify required pipeline objects
############################################################

required_objects <- c(
  "out_intermediate",
  "out_tables",
  "out_figures"
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
## 2. Define and verify input files
############################################################

genotype_file <- file.path(
  out_intermediate,
  "geno_mat_imp_accessionLevel.rds"
)

marker_file <- file.path(
  out_intermediate,
  "marker_metadata_accessionLevel.rds"
)

legacy_marker_file <- file.path(
  out_intermediate,
  "gt_filt_accessionLevelMarkers.rds"
)

blup_file <- file.path(
  out_intermediate,
  "blup_panel_all_traits.rds"
)

if (!file.exists(genotype_file)) {
  stop(
    "Accession-level genotype file not found: ",
    genotype_file
  )
}

if (!file.exists(blup_file)) {
  stop(
    "Accession-level BLUP panel not found: ",
    blup_file
  )
}

if (!file.exists(marker_file)) {
  
  if (file.exists(legacy_marker_file)) {
    
    warning(
      "Using legacy marker metadata file: ",
      legacy_marker_file
    )
    
    marker_file <- legacy_marker_file
    
  } else {
    
    stop(
      "Marker metadata file was not found. Expected either:\n",
      marker_file,
      "\nor\n",
      legacy_marker_file
    )
  }
}

############################################################
## 3. Load inputs
############################################################

geno_mat <- readRDS(
  genotype_file
)

marker_metadata <- readRDS(
  marker_file
)

blup_panel <- readRDS(
  blup_file
)

geno_mat <- as.matrix(
  geno_mat
)

storage.mode(geno_mat) <- "numeric"

if (is.null(rownames(geno_mat))) {
  stop(
    "The genotype matrix must contain accession IDs as row names."
  )
}

if (!"Accession" %in% colnames(blup_panel)) {
  stop(
    "The BLUP panel must contain an 'Accession' column."
  )
}

blup_panel$Accession <- as.character(
  blup_panel$Accession
)

############################################################
## 4. Standardize marker metadata
############################################################

chromosome_candidates <- c(
  "CHR",
  "CHROM",
  "X.CHROM",
  "Chromosome",
  "chromosome"
)

chromosome_column <- intersect(
  chromosome_candidates,
  colnames(marker_metadata)
)

if (length(chromosome_column) == 0L) {
  stop(
    "No chromosome column was found in marker metadata."
  )
}

chromosome_column <- chromosome_column[1]

position_candidates <- c(
  "BP",
  "POS",
  "Position",
  "position"
)

position_column <- intersect(
  position_candidates,
  colnames(marker_metadata)
)

if (length(position_column) == 0L) {
  stop(
    "No physical-position column was found in marker metadata."
  )
}

position_column <- position_column[1]

if ("SNP" %in% colnames(marker_metadata)) {
  
  marker_metadata$SNP <- as.character(
    marker_metadata$SNP
  )
  
} else {
  
  marker_metadata$SNP <- paste0(
    "S",
    marker_metadata[[chromosome_column]],
    "_",
    marker_metadata[[position_column]]
  )
}

marker_metadata$SNP <- make.unique(
  marker_metadata$SNP
)

marker_metadata$CHR <- as.integer(
  marker_metadata[[chromosome_column]]
)

marker_metadata$BP <- as.numeric(
  marker_metadata[[position_column]]
)

############################################################
## 5. Align marker metadata with genotype columns
############################################################

if (
  !is.null(colnames(geno_mat)) &&
  all(colnames(geno_mat) %in% marker_metadata$SNP)
) {
  
  marker_metadata <- marker_metadata[
    match(
      colnames(geno_mat),
      marker_metadata$SNP
    ),
    ,
    drop = FALSE
  ]
  
} else if (
  nrow(marker_metadata) == ncol(geno_mat)
) {
  
  warning(
    "SNP identifiers could not be matched explicitly. ",
    "Marker metadata are being aligned by row order."
  )
  
  colnames(geno_mat) <- marker_metadata$SNP
  
} else {
  
  stop(
    "Marker metadata rows do not match genotype columns. ",
    "Metadata rows = ",
    nrow(marker_metadata),
    "; genotype columns = ",
    ncol(geno_mat),
    "."
  )
}

if (
  !identical(
    colnames(geno_mat),
    marker_metadata$SNP
  )
) {
  stop(
    "Marker metadata and genotype columns are not aligned."
  )
}

if (
  anyNA(marker_metadata$CHR) ||
  anyNA(marker_metadata$BP)
) {
  stop(
    "Invalid chromosome or position values were found."
  )
}

############################################################
## 6. Align accessions
############################################################

common_accessions <- intersect(
  rownames(geno_mat),
  blup_panel$Accession
)

if (length(common_accessions) < 20L) {
  stop(
    "Too few accessions overlap between genotype and BLUP data: ",
    length(common_accessions)
  )
}

## Preserve genotype-matrix order.

common_accessions <- rownames(geno_mat)[
  rownames(geno_mat) %in%
    common_accessions
]

geno_mat <- geno_mat[
  common_accessions,
  ,
  drop = FALSE
]

blup_panel <- blup_panel[
  match(
    common_accessions,
    blup_panel$Accession
  ),
  ,
  drop = FALSE
]

if (
  !identical(
    rownames(geno_mat),
    blup_panel$Accession
  )
) {
  stop(
    "Genotype and BLUP accession orders could not be aligned."
  )
}

############################################################
## 7. Validate genotype matrix
############################################################

if (anyNA(geno_mat)) {
  stop(
    "Missing values remain in the accession-level genotype matrix."
  )
}

marker_variance <- apply(
  geno_mat,
  2,
  var,
  na.rm = TRUE
)

informative_markers <- is.finite(
  marker_variance
) & marker_variance > 0

if (!all(informative_markers)) {
  
  warning(
    "Removing ",
    sum(!informative_markers),
    " markers with zero or undefined variance."
  )
  
  geno_mat <- geno_mat[
    ,
    informative_markers,
    drop = FALSE
  ]
  
  marker_metadata <- marker_metadata[
    informative_markers,
    ,
    drop = FALSE
  ]
}

############################################################
## 8. Identify BLUP traits
############################################################

trait_columns <- grep(
  "_BLUP$",
  colnames(blup_panel),
  value = TRUE
)

if (length(trait_columns) == 0L) {
  stop(
    "No BLUP columns ending in '_BLUP' were found."
  )
}

message(
  "Marker-scan traits: ",
  paste(
    trait_columns,
    collapse = ", "
  )
)

############################################################
## 9. Analysis thresholds
############################################################

fdr_alpha <- 0.05

## These lines are displayed for visual reference.
## They are based on the number of markers tested.

bonferroni_alpha <- 0.05 / ncol(
  geno_mat
)

suggestive_alpha <- 1 / ncol(
  geno_mat
)

bonferroni_logp <- -log10(
  bonferroni_alpha
)

suggestive_logp <- -log10(
  suggestive_alpha
)

############################################################
## 10. Genomic inflation factor
############################################################

calculate_lambda_gc <- function(p_values) {
  
  valid_p <- p_values[
    is.finite(p_values) &
      !is.na(p_values) &
      p_values > 0 &
      p_values <= 1
  ]
  
  if (length(valid_p) < 100L) {
    return(NA_real_)
  }
  
  chi_square <- qchisq(
    1 - valid_p,
    df = 1
  )
  
  median(
    chi_square,
    na.rm = TRUE
  ) /
    qchisq(
      0.5,
      df = 1
    )
}

############################################################
## 11. Single-marker regression function
############################################################

run_marker_scan <- function(
    phenotype,
    genotype_matrix
) {
  
  run_one_marker <- function(marker_dosage) {
    
    marker_dosage <- as.numeric(
      marker_dosage
    )
    
    valid_records <- is.finite(
      phenotype
    ) & is.finite(
      marker_dosage
    )
    
    if (
      sum(valid_records) < 10L ||
      !is.finite(
        var(
          marker_dosage[valid_records]
        )
      ) ||
      var(
        marker_dosage[valid_records]
      ) <= 0
    ) {
      
      return(
        c(
          N = sum(valid_records),
          beta = NA_real_,
          SE = NA_real_,
          t_value = NA_real_,
          P = NA_real_,
          R2 = NA_real_
        )
      )
    }
    
    model_data <- data.frame(
      phenotype = phenotype[
        valid_records
      ],
      marker = marker_dosage[
        valid_records
      ]
    )
    
    fitted_model <- tryCatch(
      
      lm(
        phenotype ~ marker,
        data = model_data
      ),
      
      error = function(e) NULL
      
    )
    
    if (is.null(fitted_model)) {
      
      return(
        c(
          N = sum(valid_records),
          beta = NA_real_,
          SE = NA_real_,
          t_value = NA_real_,
          P = NA_real_,
          R2 = NA_real_
        )
      )
    }
    
    coefficient_table <- summary(
      fitted_model
    )$coefficients
    
    model_summary <- summary(
      fitted_model
    )
    
    c(
      N = sum(valid_records),
      
      beta = coefficient_table[
        "marker",
        "Estimate"
      ],
      
      SE = coefficient_table[
        "marker",
        "Std. Error"
      ],
      
      t_value = coefficient_table[
        "marker",
        "t value"
      ],
      
      P = coefficient_table[
        "marker",
        "Pr(>|t|)"
      ],
      
      R2 = model_summary$r.squared
    )
  }
  
  marker_results <- t(
    apply(
      genotype_matrix,
      2,
      run_one_marker
    )
  )
  
  as.data.frame(
    marker_results
  )
}

############################################################
## 12. Manhattan plot function
############################################################

draw_manhattan_plot <- function(
    scan_results,
    trait_label,
    output_file
) {
  
  plot_data <- scan_results %>%
    
    filter(
      is.finite(P),
      P > 0,
      P <= 1,
      !is.na(CHR),
      !is.na(BP)
    ) %>%
    
    arrange(
      CHR,
      BP
    )
  
  chromosome_information <- plot_data %>%
    
    group_by(CHR) %>%
    
    summarise(
      chromosome_length = max(
        BP,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    
    arrange(CHR) %>%
    
    mutate(
      offset = cumsum(
        lag(
          chromosome_length,
          default = 0
        )
      ),
      chromosome_center =
        offset +
        chromosome_length / 2
    )
  
  plot_data <- plot_data %>%
    
    left_join(
      chromosome_information,
      by = "CHR"
    ) %>%
    
    mutate(
      cumulative_position =
        BP +
        offset,
      chromosome_group =
        factor(
          CHR %% 2
        )
    )
  
  plot_object <- ggplot(
    plot_data,
    aes(
      x = cumulative_position,
      y = minus_log10_P
    )
  ) +
    
    geom_point(
      aes(
        group = chromosome_group
      ),
      size = 0.55,
      alpha = 0.70
    ) +
    
    geom_hline(
      yintercept = bonferroni_logp,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    geom_hline(
      yintercept = suggestive_logp,
      linetype = "dotted",
      linewidth = 0.5
    ) +
    
    scale_x_continuous(
      breaks =
        chromosome_information$chromosome_center,
      labels =
        chromosome_information$CHR,
      expand = expansion(
        mult = c(
          0.01,
          0.01
        )
      )
    ) +
    
    labs(
      title = paste0(
        "Complementary marker-based genome scan: ",
        trait_label
      ),
      x = "Chromosome",
      y = expression(
        -log[10](P)
      )
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      axis.text.x = element_text(
        size = 9
      )
    )
  
  ggsave(
    filename = output_file,
    plot = plot_object,
    width = 10,
    height = 5,
    dpi = 300
  )
}

############################################################
## 13. Q-Q plot function
############################################################

draw_qq_plot <- function(
    scan_results,
    trait_label,
    lambda_gc,
    output_file
) {
  
  valid_p <- scan_results$P[
    is.finite(scan_results$P) &
      scan_results$P > 0 &
      scan_results$P <= 1
  ]
  
  valid_p <- sort(
    valid_p
  )
  
  n_p <- length(
    valid_p
  )
  
  if (n_p < 10L) {
    stop(
      "Too few valid P-values to generate Q-Q plot for ",
      trait_label
    )
  }
  
  qq_data <- data.frame(
    Expected = -log10(
      ppoints(n_p)
    ),
    Observed = -log10(
      valid_p
    )
  )
  
  subtitle_text <- if (
    is.finite(lambda_gc)
  ) {
    
    paste0(
      "Genomic inflation factor, lambdaGC = ",
      format(
        lambda_gc,
        digits = 3
      )
    )
    
  } else {
    
    "Genomic inflation factor could not be estimated"
  }
  
  plot_object <- ggplot(
    qq_data,
    aes(
      x = Expected,
      y = Observed
    )
  ) +
    
    geom_point(
      size = 0.7,
      alpha = 0.65
    ) +
    
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    
    labs(
      title = paste0(
        "Q-Q plot: ",
        trait_label
      ),
      subtitle = subtitle_text,
      x = expression(
        Expected ~ -log[10](P)
      ),
      y = expression(
        Observed ~ -log[10](P)
      )
    ) +
    
    theme_classic(
      base_size = 12
    )
  
  ggsave(
    filename = output_file,
    plot = plot_object,
    width = 5.5,
    height = 5.5,
    dpi = 300
  )
}

############################################################
## 14. Run marker-based scans for all traits
############################################################

all_scan_results <- list()
run_summary_rows <- list()
qq_statistics_rows <- list()

for (trait_column in trait_columns) {
  
  message(
    "------------------------------------------------------------"
  )
  
  message(
    "Starting marker-based genome scan for: ",
    trait_column
  )
  
  phenotype_raw <- blup_panel[[
    trait_column
  ]]
  
  retained_accessions <- is.finite(
    as.numeric(
      phenotype_raw
    )
  )
  
  n_accessions <- sum(
    retained_accessions
  )
  
  if (n_accessions < 20L) {
    
    warning(
      "Skipping ",
      trait_column,
      " because only ",
      n_accessions,
      " accessions have valid BLUPs."
    )
    
    next
  }
  
  phenotype_vector <- as.numeric(
    phenotype_raw[
      retained_accessions
    ]
  )
  
  genotype_trait <- geno_mat[
    retained_accessions,
    ,
    drop = FALSE
  ]
  
  scan_statistics <- run_marker_scan(
    phenotype =
      phenotype_vector,
    genotype_matrix =
      genotype_trait
  )
  
  scan_results <- data.frame(
    Trait = trait_column,
    SNP = marker_metadata$SNP,
    CHR = marker_metadata$CHR,
    BP = marker_metadata$BP,
    N = scan_statistics$N,
    beta = scan_statistics$beta,
    SE = scan_statistics$SE,
    t_value = scan_statistics$t_value,
    R2 = scan_statistics$R2,
    P = scan_statistics$P,
    stringsAsFactors = FALSE
  )
  
  scan_results$FDR_BH <- p.adjust(
    scan_results$P,
    method = "BH"
  )
  
  scan_results$minus_log10_P <- -log10(
    scan_results$P
  )
  
  scan_results$Bonferroni_significant <- (
    scan_results$P <
      bonferroni_alpha
  )
  
  scan_results$Suggestive <- (
    scan_results$P <
      suggestive_alpha
  )
  
  scan_results$FDR_significant <- (
    scan_results$FDR_BH <
      fdr_alpha
  )
  
  scan_results <- scan_results[
    order(
      scan_results$CHR,
      scan_results$BP
    ),
    ,
    drop = FALSE
  ]
  
  lambda_gc <- calculate_lambda_gc(
    scan_results$P
  )
  
  minimum_p <- suppressWarnings(
    min(
      scan_results$P,
      na.rm = TRUE
    )
  )
  
  if (!is.finite(minimum_p)) {
    minimum_p <- NA_real_
  }
  
  n_bonferroni <- sum(
    scan_results$Bonferroni_significant,
    na.rm = TRUE
  )
  
  n_suggestive <- sum(
    scan_results$Suggestive,
    na.rm = TRUE
  )
  
  n_fdr <- sum(
    scan_results$FDR_significant,
    na.rm = TRUE
  )
  
  all_scan_results[[
    trait_column
  ]] <- scan_results
  
  ##########################################################
  ## Save trait-specific results
  ##########################################################
  
  write.csv(
    scan_results,
    file.path(
      out_tables,
      paste0(
        "marker_scan_",
        trait_column,
        "_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  saveRDS(
    scan_results,
    file.path(
      out_intermediate,
      paste0(
        "marker_scan_",
        trait_column,
        ".rds"
      )
    )
  )
  
  ##########################################################
  ## Generate figures
  ##########################################################
  
  trait_label <- gsub(
    "_BLUP$",
    "",
    trait_column
  )
  
  draw_manhattan_plot(
    scan_results =
      scan_results,
    trait_label =
      trait_label,
    output_file =
      file.path(
        out_figures,
        paste0(
          "marker_scan_",
          trait_column,
          "_Manhattan.png"
        )
      )
  )
  
  draw_qq_plot(
    scan_results =
      scan_results,
    trait_label =
      trait_label,
    lambda_gc =
      lambda_gc,
    output_file =
      file.path(
        out_figures,
        paste0(
          "marker_scan_",
          trait_column,
          "_QQ.png"
        )
      )
  )
  
  ##########################################################
  ## Save diagnostic summaries
  ##########################################################
  
  run_summary_rows[[
    trait_column
  ]] <- data.frame(
    Trait = trait_column,
    N_Accessions = n_accessions,
    N_Markers = nrow(
      scan_results
    ),
    Minimum_P = minimum_p,
    N_Bonferroni_significant =
      n_bonferroni,
    N_Suggestive =
      n_suggestive,
    N_FDR_BH_lt_0_05 =
      n_fdr,
    Lambda_GC =
      lambda_gc
  )
  
  qq_statistics_rows[[
    trait_column
  ]] <- data.frame(
    Trait = trait_column,
    Lambda_GC = lambda_gc,
    N_Accessions = n_accessions,
    N_Valid_Pvalues = sum(
      is.finite(
        scan_results$P
      ) &
        scan_results$P > 0 &
        scan_results$P <= 1
    )
  )
  
  message(
    "Completed ",
    trait_column,
    " | accessions = ",
    n_accessions,
    " | minimum P = ",
    format(
      minimum_p,
      scientific = TRUE,
      digits = 4
    ),
    " | lambdaGC = ",
    round(
      lambda_gc,
      3
    )
  )
}

############################################################
## 15. Combine complete scan results
############################################################

if (length(all_scan_results) == 0L) {
  stop(
    "No marker-based genome scans were completed."
  )
}

combined_scan_results <- bind_rows(
  all_scan_results
)

saveRDS(
  combined_scan_results,
  file.path(
    out_intermediate,
    "marker_scan_all_traits.rds"
  )
)

write.csv(
  combined_scan_results,
  file.path(
    out_tables,
    "marker_scan_all_results.csv"
  ),
  row.names = FALSE
)

############################################################
## 16. Identify strongest marker per chromosome
############################################################

chromosome_peaks <- combined_scan_results %>%
  
  filter(
    is.finite(P),
    !is.na(CHR),
    !is.na(BP)
  ) %>%
  
  group_by(
    Trait,
    CHR
  ) %>%
  
  slice_min(
    order_by = P,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup() %>%
  
  arrange(
    Trait,
    CHR
  ) %>%
  
  rename(
    Peak_SNP = SNP,
    Peak_BP = BP,
    Peak_beta = beta,
    Peak_R2 = R2,
    Peak_P = P,
    Peak_FDR_BH = FDR_BH,
    Peak_minus_log10_P =
      minus_log10_P
  )

write.csv(
  chromosome_peaks,
  file.path(
    out_tables,
    "marker_scan_peaks_per_chromosome.csv"
  ),
  row.names = FALSE
)

############################################################
## 17. Save the top 20 markers per trait
############################################################

top20_markers <- combined_scan_results %>%
  
  filter(
    is.finite(P)
  ) %>%
  
  group_by(
    Trait
  ) %>%
  
  slice_min(
    order_by = P,
    n = 20,
    with_ties = FALSE
  ) %>%
  
  ungroup() %>%
  
  arrange(
    Trait,
    P
  )

write.csv(
  top20_markers,
  file.path(
    out_tables,
    "marker_scan_top20_per_trait.csv"
  ),
  row.names = FALSE
)

############################################################
## 18. Save run summaries
############################################################

marker_scan_run_summary <- bind_rows(
  run_summary_rows
)

marker_scan_qq_statistics <- bind_rows(
  qq_statistics_rows
)

write.csv(
  marker_scan_run_summary,
  file.path(
    out_tables,
    "marker_scan_run_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  marker_scan_qq_statistics,
  file.path(
    out_tables,
    "marker_scan_QQ_statistics.csv"
  ),
  row.names = FALSE
)

############################################################
## 19. Save session information
############################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    out_intermediate,
    "sessionInfo_marker_based_genome_scan.txt"
  )
)

############################################################
## 20. Completion messages
############################################################

message(
  "05_marker_based_genome_scan.R completed successfully."
)

message(
  "Complete marker-scan results saved to: ",
  file.path(
    out_tables,
    "marker_scan_all_results.csv"
  )
)

message(
  "Chromosome-wise peaks saved to: ",
  file.path(
    out_tables,
    "marker_scan_peaks_per_chromosome.csv"
  )
)

message(
  "Figures saved in: ",
  out_figures
)