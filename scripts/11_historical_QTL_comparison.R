############################################################
##
## S3011 Perennial Sorghum GWAS Pipeline
##
## Script 11: Historical QTL Comparison
##
## Purpose:
##   - Compare S3011 rhizome-associated candidate intervals with
##     previously reported perenniality QTL intervals
##   - Quantify interval overlap and physical distance
##   - Generate a reproducible correspondence table for the
##     manuscript and supplementary materials
##   - Produce a chromosome-level comparison figure
##
## Interpretation:
##   Positional correspondence is descriptive. Overlap or proximity
##   between the present candidate intervals and historical QTL does
##   not establish identity of causal genes or independent validation.
##
## Required preceding scripts:
##   00_setup.R
##   01_genotype_QC.R
##   02_aggregate_accession_genotypes.R
##   03_BLUP_estimation.R
##   04_GWAS.R
##   05_marker_based_genome_scan.R
##   06_GWAS_marker_overlap.R
##
## Required input:
##   results/tables/GWAS_candidate_intervals_250kb.csv
##
## Curated historical QTL:
##   Historical intervals are encoded below from the sources used
##   in the manuscript. Coordinates must be checked against the
##   cited publications and confirmed to refer to a compatible
##   sorghum reference assembly before final publication.
##
## Outputs:
##   results/tables/
##     historical_perenniality_QTL_curated.csv
##     S3011_historical_QTL_pairwise_comparison.csv
##     S3011_historical_QTL_nearest_match.csv
##     S3011_historical_QTL_summary.csv
##
##   results/figures/
##     Fig_Historical_QTL_comparison.png
##
##   results/intermediate/
##     historical_QTL_comparison.rds
##     sessionInfo_historical_QTL_comparison.txt
##
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
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

## Distance used to describe a present-day candidate interval
## as being near a historical interval when they do not overlap.

proximity_tolerance_bp <- 2000000

############################################################
## 3. Define and verify input file
############################################################

candidate_interval_file <- file.path(
  out_tables,
  "GWAS_candidate_intervals_250kb.csv"
)

if (!file.exists(candidate_interval_file)) {
  stop(
    "GWAS candidate-interval file not found: ",
    candidate_interval_file,
    "\nRun 06_GWAS_marker_overlap.R first."
  )
}

############################################################
## 4. Load S3011 rhizome-associated candidate intervals
############################################################

candidate_intervals <- read_csv(
  candidate_interval_file,
  show_col_types = FALSE
)

required_interval_columns <- c(
  "Trait",
  "CHR",
  "GWAS_Interval_ID",
  "GWAS_Lead_SNP",
  "GWAS_Lead_BP",
  "GWAS_P",
  "Window_Start",
  "Window_End"
)

missing_interval_columns <- setdiff(
  required_interval_columns,
  colnames(candidate_intervals)
)

if (length(missing_interval_columns) > 0L) {
  stop(
    "Candidate-interval table is missing required columns: ",
    paste(
      missing_interval_columns,
      collapse = ", "
    )
  )
}

s3011_intervals <- candidate_intervals %>%
  mutate(
    Trait = as.character(Trait),
    CHR = as.integer(CHR),
    GWAS_Interval_ID = as.character(GWAS_Interval_ID),
    GWAS_Lead_SNP = as.character(GWAS_Lead_SNP),
    GWAS_Lead_BP = as.numeric(GWAS_Lead_BP),
    GWAS_P = as.numeric(GWAS_P),
    Window_Start = as.numeric(Window_Start),
    Window_End = as.numeric(Window_End)
  ) %>%
  filter(
    str_detect(
      Trait,
      regex(
        "rhiz",
        ignore_case = TRUE
      )
    )
  ) %>%
  distinct(
    GWAS_Interval_ID,
    .keep_all = TRUE
  ) %>%
  transmute(
    Trait,
    CHR,
    S3011_Interval_ID = GWAS_Interval_ID,
    S3011_Lead_SNP = GWAS_Lead_SNP,
    S3011_Lead_BP = GWAS_Lead_BP,
    S3011_GWAS_P = GWAS_P,
    S3011_Start_BP = Window_Start,
    S3011_End_BP = Window_End,
    S3011_Lead_Mb = GWAS_Lead_BP / 1e6,
    S3011_Start_Mb = Window_Start / 1e6,
    S3011_End_Mb = Window_End / 1e6
  ) %>%
  arrange(
    CHR,
    S3011_Lead_BP
  )

if (nrow(s3011_intervals) == 0L) {
  stop(
    "No rhizome-associated S3011 candidate intervals were found."
  )
}

############################################################
## 5. Curated historical perenniality QTL intervals
############################################################

## Trait abbreviations:
##   RG = regrowth
##   RN = rhizome number
##   RD = rhizome diameter
##
## IMPORTANT:
##   Verify interval coordinates and reference-assembly
##   compatibility against the cited publications before final
##   manuscript submission.

historical_qtl <- tribble(
  ~CHR, ~Historical_QTL, ~Historical_Trait, ~Start_Mb, ~End_Mb, ~Source,
  1L, "qRG.1C.H6.1", "RG", 66.9, 69.1, "Kong et al. 2017",
  1L, "qRG.1.H4.1",  "RG", 71.0, 71.0, "Kong et al. 2017",
  3L, "qRG.3.H6.1",  "RG", 58.0, 60.0, "Kong et al. 2017",
  3L, "qRN.3.H6.1",  "RN", 58.0, 60.0, "Kong et al. 2017",
  3L, "qRD.3.H6.1",  "RD", 58.0, 60.0, "Kong et al. 2017",
  4L, "qRG.4B.H6.1", "RG",  5.9,  8.5, "Kong et al. 2017",
  4L, "qRG.4.H6.1",  "RG",  8.9,  8.9, "Kong et al. 2017"
) %>%
  mutate(
    CHR = as.integer(CHR),
    Start_Mb = as.numeric(Start_Mb),
    End_Mb = as.numeric(End_Mb),
    Historical_Start_BP = Start_Mb * 1e6,
    Historical_End_BP = End_Mb * 1e6,
    Historical_Midpoint_BP =
      (
        Historical_Start_BP +
          Historical_End_BP
      ) / 2,
    Historical_Midpoint_Mb =
      (
        Start_Mb +
          End_Mb
      ) / 2,
    Historical_Interval_ID = paste0(
      Historical_QTL,
      "_Chr",
      CHR
    )
  ) %>%
  arrange(
    CHR,
    Historical_Start_BP,
    Historical_End_BP
  )

write_csv(
  historical_qtl,
  file.path(
    out_tables,
    "historical_perenniality_QTL_curated.csv"
  )
)

############################################################
## 6. Helper functions for interval comparison
############################################################

calculate_overlap_bp <- function(
  start_a,
  end_a,
  start_b,
  end_b
) {

  pmax(
    0,
    pmin(
      end_a,
      end_b
    ) -
      pmax(
        start_a,
        start_b
      ) +
      1
  )
}

calculate_signed_distance_bp <- function(
  start_a,
  end_a,
  start_b,
  end_b
) {

  case_when(
    end_a < start_b ~
      end_a -
        start_b,

    start_a > end_b ~
      start_a -
        end_b,

    TRUE ~
      0
  )
}

classify_relationship <- function(
  overlap_bp,
  signed_distance_bp,
  tolerance_bp
) {

  case_when(
    overlap_bp > 0 ~
      "Overlapping intervals",

    abs(
      signed_distance_bp
    ) <= tolerance_bp &
      signed_distance_bp < 0 ~
      "Upstream and within proximity tolerance",

    abs(
      signed_distance_bp
    ) <= tolerance_bp &
      signed_distance_bp > 0 ~
      "Downstream and within proximity tolerance",

    signed_distance_bp < 0 ~
      "Upstream of historical interval",

    signed_distance_bp > 0 ~
      "Downstream of historical interval",

    TRUE ~
      "Undetermined"
  )
}

############################################################
## 7. Pairwise S3011-historical interval comparison
############################################################

pairwise_comparison <- s3011_intervals %>%
  inner_join(
    historical_qtl,
    by = "CHR",
    relationship = "many-to-many"
  ) %>%
  mutate(
    Overlap_BP = calculate_overlap_bp(
      S3011_Start_BP,
      S3011_End_BP,
      Historical_Start_BP,
      Historical_End_BP
    ),

    Overlap_Mb =
      Overlap_BP /
      1e6,

    Signed_Distance_BP =
      calculate_signed_distance_bp(
        S3011_Start_BP,
        S3011_End_BP,
        Historical_Start_BP,
        Historical_End_BP
      ),

    Absolute_Distance_BP =
      abs(
        Signed_Distance_BP
      ),

    Absolute_Distance_Mb =
      Absolute_Distance_BP /
      1e6,

    Lead_to_Historical_Midpoint_BP =
      S3011_Lead_BP -
      Historical_Midpoint_BP,

    Absolute_Lead_to_Historical_Midpoint_BP =
      abs(
        Lead_to_Historical_Midpoint_BP
      ),

    S3011_Interval_Width_BP =
      S3011_End_BP -
      S3011_Start_BP +
      1,

    Historical_Interval_Width_BP =
      Historical_End_BP -
      Historical_Start_BP +
      1,

    Proportion_S3011_Interval_Overlapped =
      ifelse(
        S3011_Interval_Width_BP > 0,
        Overlap_BP /
          S3011_Interval_Width_BP,
        NA_real_
      ),

    Proportion_Historical_Interval_Overlapped =
      ifelse(
        Historical_Interval_Width_BP > 0,
        Overlap_BP /
          Historical_Interval_Width_BP,
        NA_real_
      ),

    Relationship = classify_relationship(
      overlap_bp =
        Overlap_BP,
      signed_distance_bp =
        Signed_Distance_BP,
      tolerance_bp =
        proximity_tolerance_bp
    ),

    Exact_Lead_Inside_Historical_Interval =
      S3011_Lead_BP >=
      Historical_Start_BP &
      S3011_Lead_BP <=
      Historical_End_BP,

    Within_Proximity_Tolerance =
      Overlap_BP > 0 |
      Absolute_Distance_BP <=
      proximity_tolerance_bp
  ) %>%
  arrange(
    CHR,
    S3011_Lead_BP,
    Absolute_Distance_BP,
    Historical_Start_BP
  )

write_csv(
  pairwise_comparison,
  file.path(
    out_tables,
    "S3011_historical_QTL_pairwise_comparison.csv"
  )
)

############################################################
## 8. Identify nearest historical interval per S3011 region
############################################################

nearest_match <- pairwise_comparison %>%
  group_by(
    S3011_Interval_ID
  ) %>%
  arrange(
    desc(
      Overlap_BP
    ),
    Absolute_Distance_BP,
    Absolute_Lead_to_Historical_Midpoint_BP
  ) %>%
  slice(
    1
  ) %>%
  ungroup() %>%
  mutate(
    Nearest_Match_Rank = 1L
  ) %>%
  arrange(
    CHR,
    S3011_Lead_BP
  )

write_csv(
  nearest_match,
  file.path(
    out_tables,
    "S3011_historical_QTL_nearest_match.csv"
  )
)

############################################################
## 9. Construct manuscript-style summary table
############################################################

summary_table <- pairwise_comparison %>%
  group_by(
    CHR,
    S3011_Interval_ID,
    S3011_Lead_SNP,
    S3011_Lead_BP,
    S3011_Lead_Mb,
    S3011_Start_BP,
    S3011_End_BP,
    S3011_Start_Mb,
    S3011_End_Mb
  ) %>%
  arrange(
    desc(
      Overlap_BP
    ),
    Absolute_Distance_BP,
    Historical_Start_BP
  ) %>%
  summarise(
    N_Historical_Intervals_On_Chromosome =
      n(),

    N_Overlapping_Historical_Intervals =
      sum(
        Overlap_BP > 0
      ),

    N_Proximal_Historical_Intervals =
      sum(
        Within_Proximity_Tolerance
      ),

    Nearest_Historical_QTL =
      first(
        Historical_QTL
      ),

    Nearest_Historical_Trait =
      first(
        Historical_Trait
      ),

    Nearest_Historical_Source =
      first(
        Source
      ),

    Nearest_Historical_Interval_Mb =
      paste0(
        format(
          first(
            Start_Mb
          ),
          trim = TRUE,
          nsmall = 1
        ),
        "-",
        format(
          first(
            End_Mb
          ),
          trim = TRUE,
          nsmall = 1
        )
      ),

    Nearest_Relationship =
      first(
        Relationship
      ),

    Nearest_Distance_Mb =
      first(
        Absolute_Distance_Mb
      ),

    Lead_Inside_Nearest_Historical_Interval =
      first(
        Exact_Lead_Inside_Historical_Interval
      ),

    Historical_QTLs_On_Chromosome =
      paste(
        paste0(
          Historical_QTL,
          " (",
          Historical_Trait,
          "; ",
          format(
            Start_Mb,
            trim = TRUE,
            nsmall = 1
          ),
          "-",
          format(
            End_Mb,
            trim = TRUE,
            nsmall = 1
          ),
          " Mb)"
        ),
        collapse = "; "
      ),

    .groups =
      "drop"
  ) %>%
  mutate(
    S3011_Lead_Label = paste0(
      S3011_Lead_SNP,
      " (",
      format(
        S3011_Lead_Mb,
        digits = 4,
        trim = TRUE
      ),
      " Mb)"
    ),

    S3011_Interval_Label = paste0(
      format(
        S3011_Start_Mb,
        digits = 4,
        trim = TRUE
      ),
      "-",
      format(
        S3011_End_Mb,
        digits = 4,
        trim = TRUE
      ),
      " Mb"
    ),

    Correspondence_Interpretation = case_when(
      N_Overlapping_Historical_Intervals > 0 ~
        "S3011 candidate interval overlaps at least one reported perenniality interval",

      N_Proximal_Historical_Intervals > 0 ~
        paste0(
          "S3011 candidate interval is within ",
          proximity_tolerance_bp / 1e6,
          " Mb of a reported perenniality interval"
        ),

      TRUE ~
        "No close positional correspondence under the specified tolerance"
    )
  ) %>%
  select(
    CHR,
    S3011_Interval_ID,
    S3011_Lead_Label,
    S3011_Interval_Label,
    N_Overlapping_Historical_Intervals,
    Nearest_Historical_QTL,
    Nearest_Historical_Trait,
    Nearest_Historical_Interval_Mb,
    Nearest_Relationship,
    Nearest_Distance_Mb,
    Lead_Inside_Nearest_Historical_Interval,
    Historical_QTLs_On_Chromosome,
    Correspondence_Interpretation
  ) %>%
  arrange(
    CHR
  )

write_csv(
  summary_table,
  file.path(
    out_tables,
    "S3011_historical_QTL_summary.csv"
  )
)

############################################################
## 10. Save complete comparison object
############################################################

comparison_object <- list(
  parameters = list(
    proximity_tolerance_bp =
      proximity_tolerance_bp
  ),
  s3011_intervals =
    s3011_intervals,
  historical_qtl =
    historical_qtl,
  pairwise_comparison =
    pairwise_comparison,
  nearest_match =
    nearest_match,
  summary_table =
    summary_table
)

saveRDS(
  comparison_object,
  file.path(
    out_intermediate,
    "historical_QTL_comparison.rds"
  )
)

############################################################
## 11. Prepare chromosome-level comparison figure
############################################################

s3011_plot_data <- s3011_intervals %>%
  transmute(
    CHR,
    Feature_Type =
      "S3011 candidate interval",
    Feature_Label =
      S3011_Lead_SNP,
    Start_Mb =
      S3011_Start_Mb,
    End_Mb =
      S3011_End_Mb,
    Midpoint_Mb =
      S3011_Lead_Mb
  )

historical_plot_data <- historical_qtl %>%
  transmute(
    CHR,
    Feature_Type =
      "Historical perenniality interval",
    Feature_Label =
      paste0(
        Historical_QTL,
        " (",
        Historical_Trait,
        ")"
      ),
    Start_Mb,
    End_Mb,
    Midpoint_Mb =
      Historical_Midpoint_Mb
  )

comparison_plot_data <- bind_rows(
  s3011_plot_data,
  historical_plot_data
) %>%
  group_by(
    CHR
  ) %>%
  arrange(
    Feature_Type,
    Start_Mb,
    End_Mb,
    .by_group = TRUE
  ) %>%
  mutate(
    Feature_Row =
      row_number()
  ) %>%
  ungroup()

comparison_plot <- ggplot(
  comparison_plot_data,
  aes(
    y =
      Feature_Row
  )
) +

  geom_segment(
    aes(
      x =
        Start_Mb,
      xend =
        End_Mb,
      yend =
        Feature_Row,
      linetype =
        Feature_Type
    ),
    linewidth =
      1.2
  ) +

  geom_point(
    aes(
      x =
        Midpoint_Mb,
      shape =
        Feature_Type
    ),
    size =
      2.5
  ) +

  geom_text(
    aes(
      x =
        End_Mb,
      label =
        Feature_Label
    ),
    hjust =
      -0.05,
    size =
      3
  ) +

  facet_wrap(
    ~CHR,
    scales =
      "free",
    ncol =
      1
  ) +

  scale_y_continuous(
    breaks =
      NULL
  ) +

  labs(
    title =
      "Positional comparison of S3011 candidate intervals and historical perenniality QTL",
    subtitle =
      "Line segments show reported physical intervals; points indicate interval midpoint or S3011 lead-SNP position",
    x =
      "Physical position (Mb)",
    y =
      NULL,
    linetype =
      "Feature type",
    shape =
      "Feature type"
  ) +

  coord_cartesian(
    clip =
      "off"
  ) +

  theme_classic(
    base_size =
      12
  ) +

  theme(
    plot.title =
      element_text(
        face =
          "bold"
      ),
    strip.text =
      element_text(
        face =
          "bold"
      ),
    plot.margin =
      margin(
        5.5,
        120,
        5.5,
        5.5
      ),
    legend.position =
      "bottom"
  )

ggsave(
  filename =
    file.path(
      out_figures,
      "Fig_Historical_QTL_comparison.png"
    ),
  plot =
    comparison_plot,
  width =
    10,
  height =
    8,
  dpi =
    300
)

############################################################
## 12. Save session information
############################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    out_intermediate,
    "sessionInfo_historical_QTL_comparison.txt"
  )
)

############################################################
## 13. Completion messages
############################################################

message(
  "11_historical_QTL_comparison.R completed successfully."
)

message(
  "S3011 candidate intervals analyzed: ",
  nrow(
    s3011_intervals
  )
)

message(
  "Historical intervals included: ",
  nrow(
    historical_qtl
  )
)

message(
  "Pairwise comparisons generated: ",
  nrow(
    pairwise_comparison
  )
)

message(
  "Summary table saved to: ",
  file.path(
    out_tables,
    "S3011_historical_QTL_summary.csv"
  )
)

message(
  "Comparison figure saved to: ",
  file.path(
    out_figures,
    "Fig_Historical_QTL_comparison.png"
  )
)
