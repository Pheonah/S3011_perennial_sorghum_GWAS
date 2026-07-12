############################################################
##
## S3011 Perennial Sorghum GWAS Pipeline
##
## Script 08: Local Ancestry Proxy Analysis
##
## Purpose:
##   - Orient accession-level SNP dosage so that increasing
##     dosage corresponds to the allele associated with
##     increased rhizome BLUP
##   - Summarize oriented dosage in non-overlapping 500-kb
##     genomic windows
##   - Use window-level mean dosage as an allele-dosage-derived
##     proxy for local ancestry
##   - Compare ancestry-proxy profiles between accessions with
##     high and low rhizome expression
##   - Test associations between window-level ancestry proxies
##     and rhizome BLUPs using Spearman correlation
##   - Overlay exploratory GWAS candidate intervals on the
##     ancestry-proxy profiles
##   - Save reproducible tables, intermediate objects, and figures
##
## Important interpretation:
##   This procedure does not reconstruct phased local ancestry.
##   SNP dosage is oriented according to the direction of the
##   rhizome GWAS effect and summarized within physical windows.
##   The resulting values should therefore be interpreted as
##   allele-dosage-derived ancestry proxies rather than direct
##   haplotypic assignments to Sorghum halepense or S. bicolor.
##
## Required preceding scripts:
##   00_setup.R
##   01_genotype_QC.R
##   02_aggregate_accession_genotypes.R
##   03_BLUP_estimation.R
##   04_GWAS.R
##   05_marker_based_genome_scan.R
##   06_GWAS_marker_overlap.R
##   07_candidate_gene_annotation.R
##
## Required inputs:
##   results/intermediate/
##     geno_mat_imp_accessionLevel.rds
##     marker_metadata_accessionLevel.rds
##     blup_panel_all_traits.rds
##     res_GWAS_Rhiz_BLUP.rds
##
##   results/tables/
##     GWAS_candidate_intervals_250kb.csv
##
## Outputs:
##   results/intermediate/
##     oriented_dosage_accessionLevel.rds
##     ancestry_proxy_500kb_matrix.rds
##     ancestry_proxy_500kb_long.rds
##     ancestry_proxy_window_definitions.rds
##     ancestry_proxy_high_low_difference.rds
##     sessionInfo_local_ancestry_proxy.txt
##
##   results/tables/
##     allele_dosage_orientation_rhizome.csv
##     ancestry_proxy_window_summary.csv
##     ancestry_proxy_high_low_difference.csv
##     ancestry_proxy_top_difference_windows.csv
##     ancestry_proxy_spearman_by_window.csv
##     ancestry_proxy_candidate_interval_summary.csv
##     rhizome_group_thresholds.csv
##
##   results/figures/
##     Fig_AncestryProxy_HighVsLow_Rhizome.png
##     Fig_AncestryProxy_with_GWASIntervals.png
##     Fig_AncestryProxy_Spearman_Rhizome.png
##
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(tibble)
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
## 2. Analysis parameters
############################################################

window_size_bp <- 500000L

rhizome_blup_column <- "Rhiz_BLUP"

minimum_markers_per_window <- 1L

gwas_interval_file <- file.path(
  out_tables,
  "GWAS_candidate_intervals_250kb.csv"
)

############################################################
## 3. Define required input files
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

rhizome_gwas_file <- file.path(
  out_intermediate,
  "res_GWAS_Rhiz_BLUP.rds"
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

if (!file.exists(rhizome_gwas_file)) {
  stop(
    "Rhizome GWAS result file not found: ",
    rhizome_gwas_file
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
## 4. Load inputs
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

rhizome_gwas <- readRDS(
  rhizome_gwas_file
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

if (!rhizome_blup_column %in% colnames(blup_panel)) {
  stop(
    "The BLUP panel does not contain the required rhizome column: ",
    rhizome_blup_column
  )
}

blup_panel$Accession <- as.character(
  blup_panel$Accession
)

############################################################
## 5. Standardize marker metadata
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

if (
  anyNA(marker_metadata$CHR) ||
  anyNA(marker_metadata$BP)
) {
  stop(
    "Marker metadata contain invalid chromosome or position values."
  )
}

############################################################
## 6. Align marker metadata with genotype columns
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

############################################################
## 7. Standardize and align rhizome GWAS results
############################################################

required_gwas_columns <- c(
  "SNP",
  "CHR",
  "BP",
  "beta",
  "P"
)

missing_gwas_columns <- setdiff(
  required_gwas_columns,
  colnames(rhizome_gwas)
)

if (length(missing_gwas_columns) > 0L) {
  stop(
    "Rhizome GWAS results are missing required columns: ",
    paste(
      missing_gwas_columns,
      collapse = ", "
    )
  )
}

rhizome_gwas <- rhizome_gwas %>%
  mutate(
    SNP = as.character(SNP),
    CHR = as.integer(CHR),
    BP = as.numeric(BP),
    beta = as.numeric(beta),
    P = as.numeric(P)
  )

if (!all(marker_metadata$SNP %in% rhizome_gwas$SNP)) {
  stop(
    "Some genotype markers are missing from the rhizome GWAS results."
  )
}

rhizome_gwas <- rhizome_gwas[
  match(
    marker_metadata$SNP,
    rhizome_gwas$SNP
  ),
  ,
  drop = FALSE
]

if (
  !identical(
    marker_metadata$SNP,
    rhizome_gwas$SNP
  )
) {
  stop(
    "Rhizome GWAS results could not be aligned to marker metadata."
  )
}

############################################################
## 8. Align genotype and BLUP accessions
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

if (anyNA(geno_mat)) {
  stop(
    "Missing values remain in the accession-level genotype matrix."
  )
}

############################################################
## 9. Orient SNP dosage by rhizome-effect direction
############################################################

## SNPs with negative rhizome GWAS effects are transformed as
## 2 - dosage. After orientation, increasing dosage corresponds
## to the allele associated with increased rhizome BLUP.

rhizome_beta <- rhizome_gwas$beta

flip_marker <- is.finite(
  rhizome_beta
) & rhizome_beta < 0

oriented_dosage <- geno_mat

oriented_dosage[
  ,
  flip_marker
] <- 2 - oriented_dosage[
  ,
  flip_marker,
  drop = FALSE
]

saveRDS(
  oriented_dosage,
  file.path(
    out_intermediate,
    "oriented_dosage_accessionLevel.rds"
  )
)

############################################################
## 10. Save SNP orientation and allele-frequency table
############################################################

mean_oriented_dosage <- colMeans(
  oriented_dosage,
  na.rm = TRUE
)

oriented_allele_frequency <-
  mean_oriented_dosage / 2

orientation_table <- data.frame(
  SNP = marker_metadata$SNP,
  CHR = marker_metadata$CHR,
  BP = marker_metadata$BP,
  Rhizome_GWAS_beta = rhizome_beta,
  Rhizome_GWAS_P = rhizome_gwas$P,
  Dosage_Flipped = flip_marker,
  Mean_Oriented_Dosage =
    mean_oriented_dosage,
  Frequency_Oriented_Allele =
    oriented_allele_frequency,
  Frequency_Alternative_Orientation =
    1 - oriented_allele_frequency,
  stringsAsFactors = FALSE
)

write_csv(
  orientation_table,
  file.path(
    out_tables,
    "allele_dosage_orientation_rhizome.csv"
  )
)

############################################################
## 11. Define non-overlapping 500-kb windows
############################################################

window_definitions <- bind_rows(
  lapply(
    sort(
      unique(
        marker_metadata$CHR
      )
    ),
    function(chromosome_id) {
      
      chromosome_positions <-
        marker_metadata$BP[
          marker_metadata$CHR ==
            chromosome_id
        ]
      
      chromosome_maximum <- max(
        chromosome_positions,
        na.rm = TRUE
      )
      
      window_starts <- seq(
        from = 0,
        to = chromosome_maximum,
        by = window_size_bp
      )
      
      tibble(
        CHR = chromosome_id,
        Window = seq_along(
          window_starts
        ),
        Window_Start = window_starts,
        Window_End = window_starts +
          window_size_bp -
          1,
        Window_Midpoint_BP =
          window_starts +
          window_size_bp / 2,
        Window_Midpoint_Mb =
          (
            window_starts +
              window_size_bp / 2
          ) / 1e6
      )
    }
  )
)

window_definitions$Window_ID <- paste0(
  "Chr",
  window_definitions$CHR,
  "_W",
  window_definitions$Window
)

############################################################
## 12. Assign markers to windows
############################################################

marker_window_map <- marker_metadata %>%
  transmute(
    SNP,
    CHR,
    BP,
    Window = floor(
      BP /
        window_size_bp
    ) + 1L
  ) %>%
  left_join(
    window_definitions,
    by = c(
      "CHR",
      "Window"
    )
  )

if (anyNA(marker_window_map$Window_ID)) {
  stop(
    "Some SNPs could not be assigned to a genomic window."
  )
}

markers_per_window <- marker_window_map %>%
  count(
    CHR,
    Window,
    Window_ID,
    Window_Start,
    Window_End,
    Window_Midpoint_BP,
    Window_Midpoint_Mb,
    name = "N_Markers"
  )

############################################################
## 13. Calculate accession-by-window ancestry-proxy matrix
############################################################

ancestry_proxy_matrix <- matrix(
  NA_real_,
  nrow = nrow(
    oriented_dosage
  ),
  ncol = nrow(
    window_definitions
  ),
  dimnames = list(
    rownames(
      oriented_dosage
    ),
    window_definitions$Window_ID
  )
)

for (
  window_index in
  seq_len(
    nrow(
      window_definitions
    )
  )
) {
  
  window_row <- window_definitions[
    window_index,
    ,
    drop = FALSE
  ]
  
  marker_indices <- which(
    marker_metadata$CHR ==
      window_row$CHR &
      marker_metadata$BP >=
      window_row$Window_Start &
      marker_metadata$BP <=
      window_row$Window_End
  )
  
  if (
    length(marker_indices) <
    minimum_markers_per_window
  ) {
    next
  }
  
  ancestry_proxy_matrix[
    ,
    window_index
  ] <- rowMeans(
    oriented_dosage[
      ,
      marker_indices,
      drop = FALSE
    ],
    na.rm = TRUE
  )
}

## Convert NaN values produced by completely missing windows to NA.

ancestry_proxy_matrix[
  is.nan(
    ancestry_proxy_matrix
  )
] <- NA_real_

saveRDS(
  ancestry_proxy_matrix,
  file.path(
    out_intermediate,
    "ancestry_proxy_500kb_matrix.rds"
  )
)

saveRDS(
  window_definitions,
  file.path(
    out_intermediate,
    "ancestry_proxy_window_definitions.rds"
  )
)

############################################################
## 14. Convert window matrix to long format
############################################################

ancestry_proxy_long <- as.data.frame(
  ancestry_proxy_matrix
) %>%
  rownames_to_column(
    "Accession"
  ) %>%
  pivot_longer(
    cols = -Accession,
    names_to = "Window_ID",
    values_to = "Ancestry_Proxy_Dosage"
  ) %>%
  left_join(
    window_definitions,
    by = "Window_ID"
  ) %>%
  left_join(
    markers_per_window,
    by = c(
      "CHR",
      "Window",
      "Window_ID",
      "Window_Start",
      "Window_End",
      "Window_Midpoint_BP",
      "Window_Midpoint_Mb"
    )
  )

saveRDS(
  ancestry_proxy_long,
  file.path(
    out_intermediate,
    "ancestry_proxy_500kb_long.rds"
  )
)

############################################################
## 15. Save genome-wide window summaries
############################################################

window_summary <- ancestry_proxy_long %>%
  group_by(
    CHR,
    Window,
    Window_ID,
    Window_Start,
    Window_End,
    Window_Midpoint_BP,
    Window_Midpoint_Mb,
    N_Markers
  ) %>%
  summarise(
    N_Accessions_With_Data = sum(
      is.finite(
        Ancestry_Proxy_Dosage
      )
    ),
    Mean_Ancestry_Proxy_Dosage = mean(
      Ancestry_Proxy_Dosage,
      na.rm = TRUE
    ),
    SD_Ancestry_Proxy_Dosage = sd(
      Ancestry_Proxy_Dosage,
      na.rm = TRUE
    ),
    Minimum_Ancestry_Proxy_Dosage = min(
      Ancestry_Proxy_Dosage,
      na.rm = TRUE
    ),
    Maximum_Ancestry_Proxy_Dosage = max(
      Ancestry_Proxy_Dosage,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    across(
      c(
        Mean_Ancestry_Proxy_Dosage,
        SD_Ancestry_Proxy_Dosage,
        Minimum_Ancestry_Proxy_Dosage,
        Maximum_Ancestry_Proxy_Dosage
      ),
      ~ifelse(
        is.nan(.x) |
          is.infinite(.x),
        NA_real_,
        .x
      )
    )
  ) %>%
  arrange(
    CHR,
    Window
  )

write_csv(
  window_summary,
  file.path(
    out_tables,
    "ancestry_proxy_window_summary.csv"
  )
)

############################################################
## 16. Define high- and low-rhizome accession groups
############################################################

rhizome_data <- blup_panel %>%
  transmute(
    Accession,
    Rhizome_BLUP =
      as.numeric(
        .data[[
          rhizome_blup_column
        ]]
      )
  ) %>%
  filter(
    is.finite(
      Rhizome_BLUP
    )
  )

rhizome_q25 <- quantile(
  rhizome_data$Rhizome_BLUP,
  probabilities = 0.25,
  na.rm = TRUE,
  names = FALSE
)

rhizome_median <- median(
  rhizome_data$Rhizome_BLUP,
  na.rm = TRUE
)

rhizome_q75 <- quantile(
  rhizome_data$Rhizome_BLUP,
  probabilities = 0.75,
  na.rm = TRUE,
  names = FALSE
)

rhizome_data <- rhizome_data %>%
  mutate(
    Rhizome_Group = case_when(
      Rhizome_BLUP <= rhizome_q25 ~
        "Low",
      Rhizome_BLUP >= rhizome_q75 ~
        "High",
      TRUE ~
        "Intermediate"
    )
  )

group_threshold_table <- tibble(
  Trait = rhizome_blup_column,
  Lower_Quartile_Threshold = rhizome_q25,
  Median = rhizome_median,
  Upper_Quartile_Threshold = rhizome_q75,
  N_Low = sum(
    rhizome_data$Rhizome_Group ==
      "Low"
  ),
  N_Intermediate = sum(
    rhizome_data$Rhizome_Group ==
      "Intermediate"
  ),
  N_High = sum(
    rhizome_data$Rhizome_Group ==
      "High"
  )
)

write_csv(
  group_threshold_table,
  file.path(
    out_tables,
    "rhizome_group_thresholds.csv"
  )
)

############################################################
## 17. Calculate high-minus-low ancestry-proxy differences
############################################################

ancestry_with_groups <- ancestry_proxy_long %>%
  left_join(
    rhizome_data,
    by = "Accession"
  ) %>%
  filter(
    Rhizome_Group %in%
      c(
        "High",
        "Low"
      ),
    is.finite(
      Ancestry_Proxy_Dosage
    )
  )

high_low_summary <- ancestry_with_groups %>%
  group_by(
    CHR,
    Window,
    Window_ID,
    Window_Start,
    Window_End,
    Window_Midpoint_BP,
    Window_Midpoint_Mb,
    N_Markers,
    Rhizome_Group
  ) %>%
  summarise(
    Mean_Ancestry_Proxy_Dosage = mean(
      Ancestry_Proxy_Dosage,
      na.rm = TRUE
    ),
    SD_Ancestry_Proxy_Dosage = sd(
      Ancestry_Proxy_Dosage,
      na.rm = TRUE
    ),
    N_Accessions = n_distinct(
      Accession
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Rhizome_Group,
    values_from = c(
      Mean_Ancestry_Proxy_Dosage,
      SD_Ancestry_Proxy_Dosage,
      N_Accessions
    )
  ) %>%
  mutate(
    High_minus_Low_Difference =
      Mean_Ancestry_Proxy_Dosage_High -
      Mean_Ancestry_Proxy_Dosage_Low,
    Absolute_High_minus_Low_Difference =
      abs(
        High_minus_Low_Difference
      )
  ) %>%
  arrange(
    CHR,
    Window
  )

saveRDS(
  high_low_summary,
  file.path(
    out_intermediate,
    "ancestry_proxy_high_low_difference.rds"
  )
)

write_csv(
  high_low_summary,
  file.path(
    out_tables,
    "ancestry_proxy_high_low_difference.csv"
  )
)

############################################################
## 18. Identify strongest ancestry-proxy contrasts
############################################################

top_difference_windows <- high_low_summary %>%
  filter(
    is.finite(
      Absolute_High_minus_Low_Difference
    )
  ) %>%
  group_by(
    CHR
  ) %>%
  slice_max(
    order_by =
      Absolute_High_minus_Low_Difference,
    n = 3,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    CHR,
    desc(
      Absolute_High_minus_Low_Difference
    )
  )

write_csv(
  top_difference_windows,
  file.path(
    out_tables,
    "ancestry_proxy_top_difference_windows.csv"
  )
)

############################################################
## 19. Spearman correlation with rhizome BLUP
############################################################

calculate_window_spearman <- function(
    window_data
) {
  
  valid_data <- window_data %>%
    filter(
      is.finite(
        Ancestry_Proxy_Dosage
      ),
      is.finite(
        Rhizome_BLUP
      )
    )
  
  if (
    nrow(valid_data) < 10L ||
    length(
      unique(
        valid_data$Ancestry_Proxy_Dosage
      )
    ) < 2L ||
    length(
      unique(
        valid_data$Rhizome_BLUP
      )
    ) < 2L
  ) {
    
    return(
      tibble(
        N_Accessions = nrow(
          valid_data
        ),
        Spearman_rho = NA_real_,
        P = NA_real_
      )
    )
  }
  
  correlation_test <- suppressWarnings(
    cor.test(
      valid_data$Ancestry_Proxy_Dosage,
      valid_data$Rhizome_BLUP,
      method = "spearman",
      exact = FALSE
    )
  )
  
  tibble(
    N_Accessions = nrow(
      valid_data
    ),
    Spearman_rho = unname(
      correlation_test$estimate
    ),
    P = correlation_test$p.value
  )
}

spearman_results <- ancestry_proxy_long %>%
  left_join(
    rhizome_data %>%
      select(
        Accession,
        Rhizome_BLUP
      ),
    by = "Accession"
  ) %>%
  group_by(
    CHR,
    Window,
    Window_ID,
    Window_Start,
    Window_End,
    Window_Midpoint_BP,
    Window_Midpoint_Mb,
    N_Markers
  ) %>%
  group_modify(
    ~calculate_window_spearman(.x)
  ) %>%
  ungroup() %>%
  mutate(
    FDR_BH = p.adjust(
      P,
      method = "BH"
    ),
    minus_log10_P = -log10(
      P
    )
  ) %>%
  arrange(
    CHR,
    Window
  )

write_csv(
  spearman_results,
  file.path(
    out_tables,
    "ancestry_proxy_spearman_by_window.csv"
  )
)

############################################################
## 20. Load rhizome GWAS candidate intervals
############################################################

if (file.exists(gwas_interval_file)) {
  
  gwas_candidate_intervals <- read_csv(
    gwas_interval_file,
    show_col_types = FALSE
  ) %>%
    mutate(
      Trait = as.character(
        Trait
      ),
      CHR = as.integer(
        CHR
      ),
      GWAS_Lead_BP = as.numeric(
        GWAS_Lead_BP
      ),
      Window_Start = as.numeric(
        Window_Start
      ),
      Window_End = as.numeric(
        Window_End
      )
    ) %>%
    filter(
      str_detect(
        Trait,
        regex(
          "rhiz",
          ignore_case = TRUE
        )
      )
    )
  
} else {
  
  warning(
    "GWAS candidate interval file was not found: ",
    gwas_interval_file,
    ". Candidate-interval summaries and overlays will use ",
    "the three manuscript rhizome lead SNPs."
  )
  
  gwas_candidate_intervals <- tibble(
    Trait = "Rhiz_BLUP",
    CHR = c(
      1L,
      3L,
      4L
    ),
    GWAS_Interval_ID = c(
      "Rhiz_BLUP_Chr1_62967991",
      "Rhiz_BLUP_Chr3_62120700",
      "Rhiz_BLUP_Chr4_8603936"
    ),
    GWAS_Lead_SNP = c(
      "S1_62967991",
      "S3_62120700",
      "S4_8603936"
    ),
    GWAS_Lead_BP = c(
      62967991,
      62120700,
      8603936
    ),
    GWAS_P = c(
      8.0e-6,
      1.0e-4,
      1.3e-5
    ),
    Window_Start = pmax(
      1,
      GWAS_Lead_BP -
        250000
    ),
    Window_End =
      GWAS_Lead_BP +
      250000
  )
}

############################################################
## 21. Summarize ancestry proxies in candidate intervals
############################################################

candidate_interval_window_map <- bind_rows(
  lapply(
    seq_len(
      nrow(
        gwas_candidate_intervals
      )
    ),
    function(interval_index) {
      
      interval_row <-
        gwas_candidate_intervals[
          interval_index,
          ,
          drop = FALSE
        ]
      
      window_definitions %>%
        filter(
          CHR ==
            interval_row$CHR,
          Window_End >=
            interval_row$Window_Start,
          Window_Start <=
            interval_row$Window_End
        ) %>%
        mutate(
          GWAS_Interval_ID =
            interval_row$GWAS_Interval_ID,
          GWAS_Lead_SNP =
            interval_row$GWAS_Lead_SNP,
          GWAS_Lead_BP =
            interval_row$GWAS_Lead_BP,
          GWAS_P =
            interval_row$GWAS_P
        )
    }
  )
)

candidate_interval_summary <-
  candidate_interval_window_map %>%
  left_join(
    high_low_summary,
    by = c(
      "CHR",
      "Window",
      "Window_ID",
      "Window_Start",
      "Window_End",
      "Window_Midpoint_BP",
      "Window_Midpoint_Mb"
    )
  ) %>%
  left_join(
    spearman_results %>%
      select(
        CHR,
        Window,
        Window_ID,
        Spearman_rho,
        P,
        FDR_BH
      ),
    by = c(
      "CHR",
      "Window",
      "Window_ID"
    )
  ) %>%
  rename(
    Spearman_P = P,
    Spearman_FDR_BH = FDR_BH
  ) %>%
  arrange(
    CHR,
    GWAS_Lead_BP,
    Window
  )

write_csv(
  candidate_interval_summary,
  file.path(
    out_tables,
    "ancestry_proxy_candidate_interval_summary.csv"
  )
)

############################################################
## 22. Figure: high-minus-low ancestry-proxy difference
############################################################

figure_high_low <- ggplot(
  high_low_summary,
  aes(
    x = Window_Midpoint_Mb,
    y = High_minus_Low_Difference
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  geom_line(
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    ~CHR,
    scales = "free_x",
    ncol = 2
  ) +
  
  labs(
    title =
      "Allele-dosage-derived ancestry-proxy difference",
    subtitle =
      "Mean oriented dosage in high-rhizome accessions minus low-rhizome accessions",
    x =
      "Genomic position (Mb)",
    y =
      "High minus low mean ancestry-proxy dosage"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    out_figures,
    "Fig_AncestryProxy_HighVsLow_Rhizome.png"
  ),
  plot = figure_high_low,
  width = 8,
  height = 6,
  dpi = 300
)

############################################################
## 23. Figure: ancestry proxy with GWAS intervals
############################################################

overlay_intervals <- gwas_candidate_intervals %>%
  transmute(
    CHR,
    xmin =
      Window_Start /
      1e6,
    xmax =
      Window_End /
      1e6,
    Lead_Mb =
      GWAS_Lead_BP /
      1e6,
    GWAS_Lead_SNP
  )

figure_overlay <- ggplot(
  high_low_summary,
  aes(
    x = Window_Midpoint_Mb,
    y = High_minus_Low_Difference
  )
) +
  
  geom_rect(
    data = overlay_intervals,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = -Inf,
      ymax = Inf
    ),
    inherit.aes = FALSE,
    alpha = 0.15
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  geom_line(
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  
  geom_vline(
    data = overlay_intervals,
    aes(
      xintercept = Lead_Mb
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.6
  ) +
  
  facet_wrap(
    ~CHR,
    scales = "free_x",
    ncol = 2
  ) +
  
  labs(
    title =
      "Ancestry-proxy contrasts at exploratory GWAS candidate intervals",
    subtitle =
      "Shaded regions indicate +/-250-kb candidate intervals; dotted lines indicate GWAS lead SNPs",
    x =
      "Genomic position (Mb)",
    y =
      "High minus low mean ancestry-proxy dosage"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    out_figures,
    "Fig_AncestryProxy_with_GWASIntervals.png"
  ),
  plot = figure_overlay,
  width = 8,
  height = 6,
  dpi = 300
)

############################################################
## 24. Figure: Spearman correlation by genomic window
############################################################

figure_spearman <- ggplot(
  spearman_results,
  aes(
    x = Window_Midpoint_Mb,
    y = Spearman_rho
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  geom_line(
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    ~CHR,
    scales = "free_x",
    ncol = 2
  ) +
  
  labs(
    title =
      "Association between ancestry-proxy dosage and rhizome BLUP",
    subtitle =
      "Spearman correlation calculated independently within each 500-kb genomic window",
    x =
      "Genomic position (Mb)",
    y =
      "Spearman rho"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    out_figures,
    "Fig_AncestryProxy_Spearman_Rhizome.png"
  ),
  plot = figure_spearman,
  width = 8,
  height = 6,
  dpi = 300
)

############################################################
## 25. Save session information
############################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    out_intermediate,
    "sessionInfo_local_ancestry_proxy.txt"
  )
)

############################################################
## 26. Completion messages
############################################################

message(
  "08_local_ancestry_analysis.R completed successfully."
)

message(
  "Accessions analyzed: ",
  nrow(
    oriented_dosage
  )
)

message(
  "Markers analyzed: ",
  ncol(
    oriented_dosage
  )
)

message(
  "Genomic windows defined: ",
  nrow(
    window_definitions
  )
)

message(
  "High-rhizome accessions: ",
  group_threshold_table$N_High
)

message(
  "Low-rhizome accessions: ",
  group_threshold_table$N_Low
)

message(
  "Tables saved in: ",
  out_tables
)

message(
  "Figures saved in: ",
  out_figures
)

message(
  "Intermediate objects saved in: ",
  out_intermediate
)