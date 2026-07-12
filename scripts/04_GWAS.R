############################################################
##
## S3011 Perennial Sorghum GWAS Pipeline
##
## Script 04: Genome-Wide Association Analysis
##
## Purpose:
##   - Align accession-level genotype and phenotype data
##   - Estimate population-structure covariates using PCA
##   - Conduct PCA-adjusted single-SNP association analyses
##   - Calculate SNP effects, P-values, BH-adjusted FDR,
##     proportion of variance explained (PVE), and lambdaGC
##   - Generate Manhattan and Q-Q plots for all BLUP traits
##   - Save complete association results and run summaries
##
## Association model:
##
##   BLUP = beta0 + betaSNP(SNP) +
##          beta1(PC1) + beta2(PC2) + beta3(PC3) + error
##
## Interpretation:
##   The analysis is treated as exploratory. Candidate intervals
##   are prioritized using a predefined suggestive threshold of
##   P < 1e-4. A nominal genome-wide reference line of P = 1e-5
##   is also displayed. Benjamini-Hochberg FDR values are reported
##   for every SNP.
##
## Required preceding scripts:
##   00_setup.R
##   01_genotype_QC.R
##   02_aggregate_accession_genotypes.R
##   03_BLUP_estimation.R
##
## Required inputs:
##   results/intermediate/geno_mat_imp_accessionLevel.rds
##   results/intermediate/marker_metadata_accessionLevel.rds
##   results/intermediate/blup_panel_all_traits.rds
##
## Outputs:
##   results/intermediate/res_GWAS_<trait>.rds
##   results/intermediate/GWAS_run_log.rds
##   results/tables/GWAS_<trait>_summary.csv
##   results/tables/GWAS_run_summary.csv
##   results/tables/GWAS_QQ_statistics.csv
##   results/tables/PCA_accession_scores.csv
##   results/tables/PCA_variance_explained.csv
##   results/figures/GWAS_<trait>_Manhattan.png
##   results/figures/GWAS_<trait>_QQ.png
##   results/figures/GWAS_<trait>_Manhattan_QQ.png
##
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(qqman)
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
## 2. Define input files
############################################################

genotype_file <- file.path(
  out_intermediate,
  "geno_mat_imp_accessionLevel.rds"
)

marker_file <- file.path(
  out_intermediate,
  "marker_metadata_accessionLevel.rds"
)

## Backward-compatible alternative from earlier Script 02 versions.
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
    "BLUP panel not found: ",
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

## Identify chromosome column.

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

## Identify physical-position column.

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

## Identify or construct SNP identifiers.

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
    "Genotype SNP names could not be matched explicitly. ",
    "Marker metadata are being aligned by row order."
  )
  
  colnames(geno_mat) <- marker_metadata$SNP
  
} else {
  
  stop(
    "Marker metadata rows do not match genotype matrix columns. ",
    "Metadata rows = ",
    nrow(marker_metadata),
    "; genotype SNPs = ",
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
    "Marker order does not match genotype-matrix column order."
  )
}

if (
  anyNA(marker_metadata$CHR) ||
  anyNA(marker_metadata$BP)
) {
  stop(
    "Missing or invalid chromosome/position values were found."
  )
}

############################################################
## 6. Align genotype and phenotype accessions
############################################################

common_accessions <- intersect(
  rownames(geno_mat),
  blup_panel$Accession
)

if (length(common_accessions) < 20L) {
  stop(
    "Too few overlapping accessions between genotype and BLUP data: ",
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
## 7. Check genotype matrix before PCA
############################################################

if (anyNA(geno_mat)) {
  stop(
    "Missing genotype values remain in the accession-level matrix."
  )
}

genotype_variance <- apply(
  geno_mat,
  2,
  var,
  na.rm = TRUE
)

informative_snp <- is.finite(
  genotype_variance
) & genotype_variance > 0

if (!all(informative_snp)) {
  
  warning(
    "Removing ",
    sum(!informative_snp),
    " zero-variance SNPs before PCA and GWAS."
  )
  
  geno_mat <- geno_mat[
    ,
    informative_snp,
    drop = FALSE
  ]
  
  marker_metadata <- marker_metadata[
    informative_snp,
    ,
    drop = FALSE
  ]
}

############################################################
## 8. Principal component analysis
############################################################

pca_result <- prcomp(
  geno_mat,
  center = TRUE,
  scale. = TRUE
)

if (ncol(pca_result$x) < 3L) {
  stop(
    "Fewer than three principal components were generated."
  )
}

pc_scores <- data.frame(
  Accession = rownames(geno_mat),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  stringsAsFactors = FALSE
)

pc_variance <- (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
) * 100

pc_variance_table <- data.frame(
  Principal_Component = paste0(
    "PC",
    seq_along(pc_variance)
  ),
  Variance_Explained_Percent = pc_variance
)

write.csv(
  pc_scores,
  file.path(
    out_tables,
    "PCA_accession_scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  pc_variance_table,
  file.path(
    out_tables,
    "PCA_variance_explained.csv"
  ),
  row.names = FALSE
)

phenotype_gwas <- blup_panel %>%
  left_join(
    pc_scores,
    by = "Accession"
  )

if (
  anyNA(
    phenotype_gwas[
      ,
      c("PC1", "PC2", "PC3")
    ]
  )
) {
  stop(
    "Missing PCA covariates were detected after alignment."
  )
}

############################################################
## 9. Identify BLUP traits
############################################################

trait_columns <- grep(
  "_BLUP$",
  colnames(phenotype_gwas),
  value = TRUE
)

if (length(trait_columns) == 0L) {
  stop(
    "No columns ending in '_BLUP' were found."
  )
}

message(
  "GWAS traits: ",
  paste(
    trait_columns,
    collapse = ", "
  )
)

############################################################
## 10. Analysis thresholds
############################################################

suggestive_p <- 1e-4
genomewide_reference_p <- 1e-5
fdr_alpha <- 0.05

suggestive_logp <- -log10(
  suggestive_p
)

genomewide_reference_logp <- -log10(
  genomewide_reference_p
)

############################################################
## 11. Helper: genomic inflation factor
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
## 12. Helper: Q-Q plot
############################################################

draw_qq_plot <- function(
    p_values,
    title,
    lambda_gc
) {
  
  valid_p <- p_values[
    is.finite(p_values) &
      !is.na(p_values) &
      p_values > 0 &
      p_values <= 1
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
      title
    )
  }
  
  expected <- -log10(
    ppoints(n_p)
  )
  
  observed <- -log10(
    valid_p
  )
  
  plot(
    expected,
    observed,
    pch = 16,
    cex = 0.55,
    xlab = expression(
      Expected ~ -log[10](P)
    ),
    ylab = expression(
      Observed ~ -log[10](P)
    ),
    main = title
  )
  
  abline(
    0,
    1,
    lty = 2,
    lwd = 1.2
  )
  
  if (is.finite(lambda_gc)) {
    
    legend(
      "topleft",
      legend = bquote(
        lambda[GC] == .(
          format(
            lambda_gc,
            digits = 3
          )
        )
      ),
      bty = "n"
    )
  }
}

############################################################
## 13. Helper: run one PCA-adjusted GWAS
############################################################

run_gwas_one_trait <- function(
    phenotype,
    genotype_matrix,
    pc_data
) {
  
  covariate_data <- data.frame(
    phenotype = as.numeric(
      phenotype
    ),
    PC1 = pc_data$PC1,
    PC2 = pc_data$PC2,
    PC3 = pc_data$PC3
  )
  
  ## Covariate-only null model used for PVE calculation.
  
  null_model <- lm(
    phenotype ~ PC1 + PC2 + PC3,
    data = covariate_data
  )
  
  residual_variance_null <- var(
    residuals(
      null_model
    ),
    na.rm = TRUE
  )
  
  run_single_snp <- function(snp_dosage) {
    
    snp_dosage <- as.numeric(
      snp_dosage
    )
    
    if (
      !is.finite(
        var(
          snp_dosage,
          na.rm = TRUE
        )
      ) ||
      var(
        snp_dosage,
        na.rm = TRUE
      ) <= 0
    ) {
      return(
        c(
          beta = NA_real_,
          SE = NA_real_,
          t_value = NA_real_,
          P = NA_real_,
          PVE = NA_real_
        )
      )
    }
    
    model_data <- data.frame(
      phenotype = covariate_data$phenotype,
      SNP = snp_dosage,
      PC1 = covariate_data$PC1,
      PC2 = covariate_data$PC2,
      PC3 = covariate_data$PC3
    )
    
    fitted_model <- tryCatch(
      lm(
        phenotype ~ SNP + PC1 + PC2 + PC3,
        data = model_data
      ),
      error = function(e) NULL
    )
    
    if (is.null(fitted_model)) {
      return(
        c(
          beta = NA_real_,
          SE = NA_real_,
          t_value = NA_real_,
          P = NA_real_,
          PVE = NA_real_
        )
      )
    }
    
    coefficient_table <- summary(
      fitted_model
    )$coefficients
    
    if (
      !"SNP" %in%
      rownames(coefficient_table)
    ) {
      return(
        c(
          beta = NA_real_,
          SE = NA_real_,
          t_value = NA_real_,
          P = NA_real_,
          PVE = NA_real_
        )
      )
    }
    
    residual_variance_full <- var(
      residuals(
        fitted_model
      ),
      na.rm = TRUE
    )
    
    pve <- if (
      is.finite(residual_variance_null) &&
      residual_variance_null > 0
    ) {
      
      100 * (
        residual_variance_null -
          residual_variance_full
      ) /
        residual_variance_null
      
    } else {
      
      NA_real_
    }
    
    ## Sampling variability can occasionally produce tiny negative
    ## values; these are set to zero for reporting.
    
    if (
      is.finite(pve) &&
      pve < 0 &&
      pve > -1e-8
    ) {
      pve <- 0
    }
    
    c(
      beta = coefficient_table[
        "SNP",
        "Estimate"
      ],
      SE = coefficient_table[
        "SNP",
        "Std. Error"
      ],
      t_value = coefficient_table[
        "SNP",
        "t value"
      ],
      P = coefficient_table[
        "SNP",
        "Pr(>|t|)"
      ],
      PVE = pve
    )
  }
  
  gwas_results <- t(
    apply(
      genotype_matrix,
      2,
      run_single_snp
    )
  )
  
  as.data.frame(
    gwas_results
  )
}

############################################################
## 14. Run GWAS for all traits
############################################################

gwas_run_log <- list()
gwas_summary_rows <- list()
qq_statistics_rows <- list()

for (trait_column in trait_columns) {
  
  message(
    "------------------------------------------------------------"
  )
  
  message(
    "Starting GWAS for trait: ",
    trait_column
  )
  
  phenotype_raw <- phenotype_gwas[[
    trait_column
  ]]
  
  retained_accessions <- !is.na(
    phenotype_raw
  ) & is.finite(
    as.numeric(
      phenotype_raw
    )
  )
  
  n_used <- sum(
    retained_accessions
  )
  
  if (n_used < 30L) {
    
    warning(
      "Skipping ",
      trait_column,
      " because only ",
      n_used,
      " accessions have non-missing BLUPs."
    )
    
    gwas_run_log[[
      trait_column
    ]] <- list(
      skipped = TRUE,
      n = n_used
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
  
  pc_trait <- phenotype_gwas[
    retained_accessions,
    c("PC1", "PC2", "PC3"),
    drop = FALSE
  ]
  
  association_statistics <- run_gwas_one_trait(
    phenotype = phenotype_vector,
    genotype_matrix = genotype_trait,
    pc_data = pc_trait
  )
  
  association_results <- data.frame(
    SNP = marker_metadata$SNP,
    CHR = marker_metadata$CHR,
    BP = marker_metadata$BP,
    beta = association_statistics$beta,
    SE = association_statistics$SE,
    t_value = association_statistics$t_value,
    P = association_statistics$P,
    PVE_percent = association_statistics$PVE,
    stringsAsFactors = FALSE
  )
  
  association_results$FDR_BH <- p.adjust(
    association_results$P,
    method = "BH"
  )
  
  association_results$minus_log10_P <- -log10(
    association_results$P
  )
  
  association_results$Suggestive <- (
    association_results$P <
      suggestive_p
  )
  
  association_results$FDR_significant <- (
    association_results$FDR_BH <
      fdr_alpha
  )
  
  association_results <- association_results[
    order(
      association_results$CHR,
      association_results$BP
    ),
    ,
    drop = FALSE
  ]
  
  lambda_gc <- calculate_lambda_gc(
    association_results$P
  )
  
  minimum_p <- suppressWarnings(
    min(
      association_results$P,
      na.rm = TRUE
    )
  )
  
  if (!is.finite(minimum_p)) {
    minimum_p <- NA_real_
  }
  
  n_suggestive <- sum(
    association_results$Suggestive,
    na.rm = TRUE
  )
  
  n_fdr_significant <- sum(
    association_results$FDR_significant,
    na.rm = TRUE
  )
  
  ##########################################################
  ## Save complete association outputs
  ##########################################################
  
  output_rds <- file.path(
    out_intermediate,
    paste0(
      "res_GWAS_",
      trait_column,
      ".rds"
    )
  )
  
  output_csv <- file.path(
    out_tables,
    paste0(
      "GWAS_",
      trait_column,
      "_summary.csv"
    )
  )
  
  saveRDS(
    association_results,
    output_rds
  )
  
  write.csv(
    association_results,
    output_csv,
    row.names = FALSE
  )
  
  ##########################################################
  ## Prepare plotting data
  ##########################################################
  
  plot_data <- association_results[
    is.finite(association_results$P) &
      association_results$P > 0 &
      association_results$P <= 1,
    c("SNP", "CHR", "BP", "P"),
    drop = FALSE
  ]
  
  maximum_logp <- max(
    -log10(
      plot_data$P
    ),
    na.rm = TRUE
  )
  
  y_axis_max <- max(
    maximum_logp + 0.5,
    genomewide_reference_logp + 0.5
  )
  
  plot_title <- gsub(
    "_BLUP$",
    "",
    trait_column
  )
  
  ##########################################################
  ## Combined Manhattan and Q-Q figure
  ##########################################################
  
  combined_figure <- file.path(
    out_figures,
    paste0(
      "GWAS_",
      trait_column,
      "_Manhattan_QQ.png"
    )
  )
  
  png(
    combined_figure,
    width = 7,
    height = 8,
    units = "in",
    res = 300
  )
  
  par(
    mfrow = c(2, 1),
    mar = c(4.5, 4.5, 3.2, 1.0)
  )
  
  manhattan(
    plot_data,
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    main = paste0(
      plot_title,
      ": PCA-adjusted GWAS"
    ),
    ylim = c(
      0,
      y_axis_max
    ),
    genomewideline =
      genomewide_reference_logp,
    suggestiveline =
      suggestive_logp,
    cex = 0.55,
    cex.axis = 0.8
  )
  
  draw_qq_plot(
    p_values = plot_data$P,
    title = paste0(
      plot_title,
      ": Q-Q plot"
    ),
    lambda_gc = lambda_gc
  )
  
  dev.off()
  
  ##########################################################
  ## Standalone Manhattan figure
  ##########################################################
  
  manhattan_figure <- file.path(
    out_figures,
    paste0(
      "GWAS_",
      trait_column,
      "_Manhattan.png"
    )
  )
  
  png(
    manhattan_figure,
    width = 9,
    height = 4.5,
    units = "in",
    res = 300
  )
  
  par(
    mfrow = c(1, 1),
    mar = c(4.5, 4.5, 3.2, 1.0)
  )
  
  manhattan(
    plot_data,
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    main = paste0(
      plot_title,
      ": PCA-adjusted GWAS"
    ),
    ylim = c(
      0,
      y_axis_max
    ),
    genomewideline =
      genomewide_reference_logp,
    suggestiveline =
      suggestive_logp,
    cex = 0.55,
    cex.axis = 0.8
  )
  
  dev.off()
  
  ##########################################################
  ## Standalone Q-Q figure
  ##########################################################
  
  qq_figure <- file.path(
    out_figures,
    paste0(
      "GWAS_",
      trait_column,
      "_QQ.png"
    )
  )
  
  png(
    qq_figure,
    width = 5.5,
    height = 5.5,
    units = "in",
    res = 300
  )
  
  par(
    mfrow = c(1, 1),
    mar = c(4.5, 4.5, 3.2, 1.0)
  )
  
  draw_qq_plot(
    p_values = plot_data$P,
    title = paste0(
      plot_title,
      ": Q-Q plot"
    ),
    lambda_gc = lambda_gc
  )
  
  dev.off()
  
  ##########################################################
  ## Save run diagnostics
  ##########################################################
  
  gwas_run_log[[
    trait_column
  ]] <- list(
    skipped = FALSE,
    n = n_used,
    lambda_GC = lambda_gc,
    minimum_P = minimum_p,
    number_suggestive = n_suggestive,
    number_FDR_significant = n_fdr_significant,
    output_rds = output_rds,
    output_csv = output_csv
  )
  
  gwas_summary_rows[[
    trait_column
  ]] <- data.frame(
    Trait = trait_column,
    N_Accessions = n_used,
    N_SNPs = nrow(association_results),
    Minimum_P = minimum_p,
    N_Suggestive_P_lt_1e_4 =
      n_suggestive,
    N_FDR_BH_lt_0_05 =
      n_fdr_significant,
    Lambda_GC = lambda_gc
  )
  
  qq_statistics_rows[[
    trait_column
  ]] <- data.frame(
    Trait = trait_column,
    Lambda_GC = lambda_gc,
    N_Accessions = n_used,
    N_Valid_Pvalues = sum(
      is.finite(
        association_results$P
      ) &
        association_results$P > 0 &
        association_results$P <= 1
    )
  )
  
  message(
    "Completed ",
    trait_column,
    " | n = ",
    n_used,
    " | lambdaGC = ",
    round(
      lambda_gc,
      3
    ),
    " | minimum P = ",
    format(
      minimum_p,
      scientific = TRUE,
      digits = 4
    ),
    " | suggestive SNPs = ",
    n_suggestive,
    " | FDR-significant SNPs = ",
    n_fdr_significant
  )
}

############################################################
## 15. Save run summaries
############################################################

saveRDS(
  gwas_run_log,
  file.path(
    out_intermediate,
    "GWAS_run_log.rds"
  )
)

if (length(gwas_summary_rows) > 0L) {
  
  gwas_run_summary <- bind_rows(
    gwas_summary_rows
  )
  
  write.csv(
    gwas_run_summary,
    file.path(
      out_tables,
      "GWAS_run_summary.csv"
    ),
    row.names = FALSE
  )
}

if (length(qq_statistics_rows) > 0L) {
  
  qq_statistics <- bind_rows(
    qq_statistics_rows
  )
  
  write.csv(
    qq_statistics,
    file.path(
      out_tables,
      "GWAS_QQ_statistics.csv"
    ),
    row.names = FALSE
  )
}

############################################################
## 16. Completion messages
############################################################

message(
  "04_GWAS.R completed successfully."
)

message(
  "GWAS tables saved in: ",
  out_tables
)

message(
  "GWAS figures saved in: ",
  out_figures
)

message(
  "GWAS intermediate objects saved in: ",
  out_intermediate
)