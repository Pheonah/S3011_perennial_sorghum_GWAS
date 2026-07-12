############################################################
##
## S3011 Perennial Sorghum GWAS Pipeline
##
## Script 06: GWAS–Marker-Scan Overlap Analysis
##
## Purpose:
##   - Identify candidate intervals from the exploratory GWAS
##   - Identify chromosome-wise peaks from the complementary
##     marker-based genome scan
##   - Compare marker-scan peaks with GWAS candidate intervals
##     using a +/-250-kb physical window
##   - Calculate exact SNP matches and physical distances
##   - Summarize concordance between analytical frameworks
##   - Generate a figure displaying overlapping intervals
##
## Analytical interpretation:
##   Concordance between the exploratory GWAS and complementary
##   marker-based genome scan is interpreted as complementary
##   methodological support for prioritizing candidate intervals.
##   Because both analyses use the same population and genotype
##   dataset, overlap is not interpreted as independent validation
##   of causal loci.
##
## Required preceding scripts:
##   00_setup.R
##   01_genotype_QC.R
##   02_aggregate_accession_genotypes.R
##   03_BLUP_estimation.R
##   04_GWAS.R
##   05_marker_based_genome_scan.R
##
## Required inputs:
##   results/intermediate/res_GWAS_<trait>.rds
##   results/intermediate/marker_scan_all_traits.rds
##
## Outputs:
##   results/tables/GWAS_candidate_intervals_250kb.csv
##   results/tables/marker_scan_chromosome_peaks.csv
##   results/tables/GWAS_marker_overlap_250kb.csv
##   results/tables/GWAS_marker_overlap_summary.csv
##   results/tables/GWAS_marker_exact_matches.csv
##   results/figures/GWAS_marker_overlap_250kb.png
##   results/intermediate/GWAS_marker_overlap_250kb.rds
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
## 2. Analysis parameters
############################################################

## Predefined suggestive GWAS threshold used in the manuscript.

gwas_p_threshold <- 1e-4

## Physical window used to define candidate intervals.

overlap_window_bp <- 250000L

############################################################
## 3. Locate GWAS and marker-scan files
############################################################

gwas_files <- list.files(
  path = out_intermediate,
  pattern = "^res_GWAS_.*_BLUP\\.rds$",
  full.names = TRUE
)

if (length(gwas_files) == 0L) {
  stop(
    "No GWAS result files were found in: ",
    out_intermediate,
    "\nExpected files named res_GWAS_<trait>_BLUP.rds."
  )
}

marker_scan_file <- file.path(
  out_intermediate,
  "marker_scan_all_traits.rds"
)

if (!file.exists(marker_scan_file)) {
  stop(
    "Combined marker-scan file not found: ",
    marker_scan_file,
    "\nRun 05_marker_based_genome_scan.R first."
  )
}

############################################################
## 4. Helper functions
############################################################

## Extract the trait name from a GWAS filename.

extract_trait_name <- function(file_path) {
  
  file_name <- basename(
    file_path
  )
  
  sub(
    "^res_GWAS_(.*)\\.rds$",
    "\\1",
    file_name
  )
}

## Standardize GWAS result columns.

standardize_gwas_columns <- function(
    gwas_data,
    trait_name
) {
  
  required_columns <- c(
    "SNP",
    "CHR",
    "BP",
    "P"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(gwas_data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "GWAS file for trait ",
      trait_name,
      " is missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  gwas_data %>%
    
    transmute(
      Trait = trait_name,
      SNP = as.character(SNP),
      CHR = as.integer(CHR),
      BP = as.numeric(BP),
      P = as.numeric(P),
      beta = if (
        "beta" %in% colnames(gwas_data)
      ) {
        as.numeric(beta)
      } else {
        NA_real_
      },
      FDR_BH = if (
        "FDR_BH" %in% colnames(gwas_data)
      ) {
        as.numeric(FDR_BH)
      } else {
        p.adjust(
          as.numeric(P),
          method = "BH"
        )
      }
    )
}

## Declump GWAS signals by retaining the most significant SNP
## and removing other suggestive SNPs within +/-window_bp.

declump_gwas_signals <- function(
    gwas_data,
    p_threshold,
    window_bp
) {
  
  suggestive_signals <- gwas_data %>%
    
    filter(
      is.finite(P),
      P > 0,
      P < p_threshold,
      !is.na(CHR),
      !is.na(BP)
    ) %>%
    
    arrange(
      CHR,
      P,
      BP
    )
  
  if (nrow(suggestive_signals) == 0L) {
    return(
      suggestive_signals
    )
  }
  
  retained_leads <- list()
  
  for (
    chromosome_id in
    sort(
      unique(
        suggestive_signals$CHR
      )
    )
  ) {
    
    chromosome_data <- suggestive_signals %>%
      
      filter(
        CHR == chromosome_id
      ) %>%
      
      arrange(
        P,
        BP
      )
    
    chromosome_leads <- list()
    
    while (nrow(chromosome_data) > 0L) {
      
      lead_snp <- chromosome_data[
        1,
        ,
        drop = FALSE
      ]
      
      chromosome_leads[[
        length(chromosome_leads) + 1L
      ]] <- lead_snp
      
      lead_position <- lead_snp$BP[1]
      
      chromosome_data <- chromosome_data %>%
        
        filter(
          abs(
            BP -
              lead_position
          ) >
            window_bp
        )
    }
    
    retained_leads[[
      as.character(
        chromosome_id
      )
    ]] <- bind_rows(
      chromosome_leads
    )
  }
  
  bind_rows(
    retained_leads
  )
}

############################################################
## 5. Load and combine GWAS results
############################################################

gwas_results_list <- lapply(
  gwas_files,
  function(gwas_file) {
    
    trait_name <- extract_trait_name(
      gwas_file
    )
    
    gwas_data <- readRDS(
      gwas_file
    )
    
    standardize_gwas_columns(
      gwas_data =
        gwas_data,
      trait_name =
        trait_name
    )
  }
)

gwas_results <- bind_rows(
  gwas_results_list
)

message(
  "GWAS traits detected: ",
  paste(
    sort(
      unique(
        gwas_results$Trait
      )
    ),
    collapse = ", "
  )
)

############################################################
## 6. Define GWAS candidate intervals
############################################################

gwas_candidate_leads <- bind_rows(
  lapply(
    split(
      gwas_results,
      gwas_results$Trait
    ),
    function(trait_data) {
      
      declump_gwas_signals(
        gwas_data =
          trait_data,
        p_threshold =
          gwas_p_threshold,
        window_bp =
          overlap_window_bp
      )
    }
  )
)

if (nrow(gwas_candidate_leads) == 0L) {
  stop(
    "No GWAS candidate signals were detected at P < ",
    gwas_p_threshold,
    "."
  )
}

gwas_candidate_intervals <- gwas_candidate_leads %>%
  
  mutate(
    GWAS_Lead_SNP = SNP,
    GWAS_Lead_BP = BP,
    GWAS_P = P,
    GWAS_beta = beta,
    GWAS_FDR_BH = FDR_BH,
    Window_Start = pmax(
      1,
      BP -
        overlap_window_bp
    ),
    Window_End =
      BP +
      overlap_window_bp,
    GWAS_Interval_ID = paste0(
      Trait,
      "_Chr",
      CHR,
      "_",
      BP
    )
  ) %>%
  
  select(
    Trait,
    CHR,
    GWAS_Interval_ID,
    GWAS_Lead_SNP,
    GWAS_Lead_BP,
    GWAS_P,
    GWAS_beta,
    GWAS_FDR_BH,
    Window_Start,
    Window_End
  ) %>%
  
  arrange(
    Trait,
    CHR,
    GWAS_Lead_BP
  )

write.csv(
  gwas_candidate_intervals,
  file.path(
    out_tables,
    "GWAS_candidate_intervals_250kb.csv"
  ),
  row.names = FALSE
)

############################################################
## 7. Load marker-based genome-scan results
############################################################

marker_scan_results <- readRDS(
  marker_scan_file
)

required_marker_columns <- c(
  "Trait",
  "SNP",
  "CHR",
  "BP",
  "P"
)

missing_marker_columns <- setdiff(
  required_marker_columns,
  colnames(marker_scan_results)
)

if (length(missing_marker_columns) > 0L) {
  stop(
    "Marker-scan results are missing required columns: ",
    paste(
      missing_marker_columns,
      collapse = ", "
    )
  )
}

marker_scan_results <- marker_scan_results %>%
  
  mutate(
    Trait = as.character(Trait),
    SNP = as.character(SNP),
    CHR = as.integer(CHR),
    BP = as.numeric(BP),
    P = as.numeric(P),
    FDR_BH = if (
      "FDR_BH" %in%
      colnames(marker_scan_results)
    ) {
      as.numeric(FDR_BH)
    } else {
      ave(
        P,
        Trait,
        FUN = function(x) {
          p.adjust(
            x,
            method = "BH"
          )
        }
      )
    }
  )

############################################################
## 8. Identify strongest marker-scan peak per chromosome
############################################################

marker_chromosome_peaks <- marker_scan_results %>%
  
  filter(
    is.finite(P),
    P > 0,
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
  
  transmute(
    Trait,
    CHR,
    Marker_Peak_SNP = SNP,
    Marker_Peak_BP = BP,
    Marker_P = P,
    Marker_beta = if (
      "beta" %in%
      colnames(marker_scan_results)
    ) {
      beta
    } else {
      NA_real_
    },
    Marker_R2 = if (
      "R2" %in%
      colnames(marker_scan_results)
    ) {
      R2
    } else {
      NA_real_
    },
    Marker_FDR_BH = FDR_BH
  ) %>%
  
  arrange(
    Trait,
    CHR
  )

write.csv(
  marker_chromosome_peaks,
  file.path(
    out_tables,
    "marker_scan_chromosome_peaks.csv"
  ),
  row.names = FALSE
)

############################################################
## 9. Compare GWAS intervals with marker-scan peaks
############################################################

overlap_results <- gwas_candidate_intervals %>%
  
  inner_join(
    marker_chromosome_peaks,
    by = c(
      "Trait",
      "CHR"
    )
  ) %>%
  
  mutate(
    Distance_bp = abs(
      Marker_Peak_BP -
        GWAS_Lead_BP
    ),
    Within_250kb = (
      Marker_Peak_BP >=
        Window_Start &
        Marker_Peak_BP <=
        Window_End
    ),
    Exact_Position_Match = (
      Marker_Peak_BP ==
        GWAS_Lead_BP
    ),
    Exact_SNP_Match = (
      Marker_Peak_SNP ==
        GWAS_Lead_SNP
    )
  ) %>%
  
  arrange(
    Trait,
    CHR,
    Distance_bp
  )

write.csv(
  overlap_results,
  file.path(
    out_tables,
    "GWAS_marker_overlap_250kb.csv"
  ),
  row.names = FALSE
)

saveRDS(
  overlap_results,
  file.path(
    out_intermediate,
    "GWAS_marker_overlap_250kb.rds"
  )
)

############################################################
## 10. Extract overlapping and exact-match results
############################################################

overlapping_intervals <- overlap_results %>%
  
  filter(
    Within_250kb
  )

exact_matches <- overlap_results %>%
  
  filter(
    Exact_Position_Match |
      Exact_SNP_Match
  )

write.csv(
  exact_matches,
  file.path(
    out_tables,
    "GWAS_marker_exact_matches.csv"
  ),
  row.names = FALSE
)

############################################################
## 11. Summarize overlap by trait
############################################################

gwas_interval_counts <- gwas_candidate_intervals %>%
  
  count(
    Trait,
    name = "N_GWAS_candidate_intervals"
  )

marker_peak_counts <- marker_chromosome_peaks %>%
  
  count(
    Trait,
    name = "N_marker_chromosome_peaks"
  )

overlap_counts <- overlap_results %>%
  
  group_by(
    Trait
  ) %>%
  
  summarise(
    N_GWAS_intervals_compared =
      n_distinct(
        GWAS_Interval_ID
      ),
    N_intervals_with_marker_peak_250kb =
      n_distinct(
        GWAS_Interval_ID[
          Within_250kb
        ]
      ),
    N_exact_position_matches =
      sum(
        Exact_Position_Match,
        na.rm = TRUE
      ),
    N_exact_SNP_matches =
      sum(
        Exact_SNP_Match,
        na.rm = TRUE
      ),
    Minimum_distance_bp =
      if (
        any(
          is.finite(
            Distance_bp
          )
        )
      ) {
        min(
          Distance_bp,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    .groups = "drop"
  )

overlap_summary <- gwas_interval_counts %>%
  
  full_join(
    marker_peak_counts,
    by = "Trait"
  ) %>%
  
  full_join(
    overlap_counts,
    by = "Trait"
  ) %>%
  
  mutate(
    across(
      where(is.numeric),
      ~replace(
        .x,
        is.na(.x),
        0
      )
    ),
    Percent_GWAS_intervals_overlapping = ifelse(
      N_GWAS_candidate_intervals > 0,
      round(
        100 *
          N_intervals_with_marker_peak_250kb /
          N_GWAS_candidate_intervals,
        1
      ),
      NA_real_
    )
  ) %>%
  
  arrange(
    Trait
  )

write.csv(
  overlap_summary,
  file.path(
    out_tables,
    "GWAS_marker_overlap_summary.csv"
  ),
  row.names = FALSE
)

############################################################
## 12. Generate overlap figure
############################################################

figure_data <- overlap_results %>%
  
  filter(
    Within_250kb
  ) %>%
  
  mutate(
    Trait_Label = gsub(
      "_BLUP$",
      "",
      Trait
    ),
    Interval_Label = paste0(
      "Chr",
      CHR,
      ": ",
      GWAS_Lead_SNP
    ),
    Distance_kb =
      Distance_bp /
      1000
  )

if (nrow(figure_data) > 0L) {
  
  overlap_figure <- ggplot(
    figure_data,
    aes(
      x = Distance_kb,
      y = reorder(
        Interval_Label,
        Distance_kb
      )
    )
  ) +
    
    geom_point(
      size = 3
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    facet_wrap(
      ~Trait_Label,
      scales = "free_y"
    ) +
    
    labs(
      title =
        "Concordance between exploratory GWAS intervals and marker-scan peaks",
      subtitle =
        "Points show physical distance between the GWAS lead SNP and chromosome-wise marker-scan peak",
      x =
        "Distance between lead markers (kb)",
      y =
        "GWAS candidate interval"
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
      "GWAS_marker_overlap_250kb.png"
    ),
    plot = overlap_figure,
    width = 9,
    height = 6,
    dpi = 300
  )
  
} else {
  
  warning(
    "No marker-scan peaks occurred within +/-250 kb of ",
    "GWAS candidate intervals; overlap figure was not generated."
  )
}

############################################################
## 13. Save session information
############################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    out_intermediate,
    "sessionInfo_GWAS_marker_overlap.txt"
  )
)

############################################################
## 14. Completion messages
############################################################

message(
  "06_GWAS_marker_overlap.R completed successfully."
)

message(
  "GWAS candidate intervals: ",
  nrow(
    gwas_candidate_intervals
  )
)

message(
  "Marker-scan chromosome peaks: ",
  nrow(
    marker_chromosome_peaks
  )
)

message(
  "Overlapping intervals within +/-",
  overlap_window_bp / 1000,
  " kb: ",
  nrow(
    overlapping_intervals
  )
)

message(
  "Exact SNP matches: ",
  nrow(
    exact_matches
  )
)

message(
  "Tables saved in: ",
  out_tables
)

message(
  "Figures saved in: ",
  out_figures
)