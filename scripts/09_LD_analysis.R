############################################################
##
## S3011 Perennial Sorghum GWAS Pipeline
##
## Script 09: Linkage Disequilibrium Analysis
##
## Purpose:
##   - Characterize local linkage disequilibrium (LD) within the
##     rhizome-associated candidate intervals identified by the
##     exploratory GWAS
##   - Calculate pairwise SNP LD as squared Pearson correlation
##     (r^2) using accession-level mean allele dosage
##   - Calculate LD between each marker and the GWAS lead SNP
##   - Summarize LD-supported regions at r^2 >= 0.50 and
##     r^2 >= 0.70
##   - Generate lead-marker LD profiles and interval-specific
##     pairwise LD heatmaps
##
## Interpretation:
##   The S3011 population is an advanced diploid-extraction panel
##   with family structure and accession-level mean allele dosage.
##   LD summaries are therefore used to describe local marker
##   correlation around prioritized candidate intervals, not to
##   define formal linkage-map QTL confidence intervals.
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
##   08_local_ancestry_analysis.R
##
## Required inputs:
##   results/intermediate/
##     geno_mat_imp_accessionLevel.rds
##     marker_metadata_accessionLevel.rds
##
##   results/tables/
##     GWAS_candidate_intervals_250kb.csv
##
## Outputs:
##   results/intermediate/
##     LD_matrices_rhizome_intervals.rds
##     sessionInfo_LD_analysis.txt
##
##   results/tables/
##     LD_marker_to_lead_all_intervals.csv
##     LD_summary_r2_0.50_rhizome_intervals.csv
##     LD_summary_r2_0.70_rhizome_intervals.csv
##     LD_interval_marker_counts.csv
##
##   results/figures/
##     Fig_LD_to_lead_rhizome_intervals.png
##     Fig_LD_heatmap_<interval>.png
##
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
})

required_objects <- c("out_intermediate", "out_tables", "out_figures")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]
if (length(missing_objects) > 0L) {
  stop("Required pipeline objects are missing: ", paste(missing_objects, collapse = ", "),
       ". Run 00_setup.R before running this script.")
}

primary_r2_threshold <- 0.50
stringent_r2_threshold <- 0.70
minimum_markers_per_interval <- 2L

genotype_file <- file.path(out_intermediate, "geno_mat_imp_accessionLevel.rds")
marker_file <- file.path(out_intermediate, "marker_metadata_accessionLevel.rds")
legacy_marker_file <- file.path(out_intermediate, "gt_filt_accessionLevelMarkers.rds")
candidate_interval_file <- file.path(out_tables, "GWAS_candidate_intervals_250kb.csv")

if (!file.exists(genotype_file)) stop("Accession-level genotype file not found: ", genotype_file)
if (!file.exists(candidate_interval_file)) {
  stop("GWAS candidate-interval file not found: ", candidate_interval_file,
       "\nRun 06_GWAS_marker_overlap.R first.")
}
if (!file.exists(marker_file)) {
  if (file.exists(legacy_marker_file)) {
    warning("Using legacy marker metadata file: ", legacy_marker_file)
    marker_file <- legacy_marker_file
  } else {
    stop("Marker metadata file was not found. Expected either:\n", marker_file,
         "\nor\n", legacy_marker_file)
  }
}

geno_mat <- as.matrix(readRDS(genotype_file))
storage.mode(geno_mat) <- "numeric"
marker_metadata <- readRDS(marker_file)
candidate_intervals <- read_csv(candidate_interval_file, show_col_types = FALSE)
if (is.null(rownames(geno_mat))) stop("The genotype matrix must contain accession IDs as row names.")

chromosome_candidates <- c("CHR", "CHROM", "X.CHROM", "Chromosome", "chromosome")
chromosome_column <- intersect(chromosome_candidates, colnames(marker_metadata))
if (length(chromosome_column) == 0L) stop("No chromosome column was found in marker metadata.")
chromosome_column <- chromosome_column[1]
position_candidates <- c("BP", "POS", "Position", "position")
position_column <- intersect(position_candidates, colnames(marker_metadata))
if (length(position_column) == 0L) stop("No physical-position column was found in marker metadata.")
position_column <- position_column[1]

if ("SNP" %in% colnames(marker_metadata)) {
  marker_metadata$SNP <- as.character(marker_metadata$SNP)
} else {
  marker_metadata$SNP <- paste0("S", marker_metadata[[chromosome_column]], "_", marker_metadata[[position_column]])
}
marker_metadata$SNP <- make.unique(marker_metadata$SNP)
marker_metadata$CHR <- as.integer(marker_metadata[[chromosome_column]])
marker_metadata$BP <- as.numeric(marker_metadata[[position_column]])
if (anyNA(marker_metadata$CHR) || anyNA(marker_metadata$BP)) {
  stop("Marker metadata contain invalid chromosome or position values.")
}

if (!is.null(colnames(geno_mat)) && all(colnames(geno_mat) %in% marker_metadata$SNP)) {
  marker_metadata <- marker_metadata[match(colnames(geno_mat), marker_metadata$SNP), , drop = FALSE]
} else if (nrow(marker_metadata) == ncol(geno_mat)) {
  warning("SNP identifiers could not be matched explicitly. Marker metadata are being aligned by row order.")
  colnames(geno_mat) <- marker_metadata$SNP
} else {
  stop("Marker metadata rows do not match genotype columns. Metadata rows = ", nrow(marker_metadata),
       "; genotype columns = ", ncol(geno_mat), ".")
}
if (!identical(colnames(geno_mat), marker_metadata$SNP)) stop("Marker metadata and genotype columns are not aligned.")
if (anyNA(geno_mat)) stop("Missing values remain in the accession-level genotype matrix.")

required_interval_columns <- c("Trait", "CHR", "GWAS_Interval_ID", "GWAS_Lead_SNP",
                               "GWAS_Lead_BP", "GWAS_P", "Window_Start", "Window_End")
missing_interval_columns <- setdiff(required_interval_columns, colnames(candidate_intervals))
if (length(missing_interval_columns) > 0L) {
  stop("Candidate-interval table is missing required columns: ",
       paste(missing_interval_columns, collapse = ", "))
}

rhizome_intervals <- candidate_intervals %>%
  mutate(
    Trait = as.character(Trait), CHR = as.integer(CHR),
    GWAS_Interval_ID = as.character(GWAS_Interval_ID),
    GWAS_Lead_SNP = as.character(GWAS_Lead_SNP),
    GWAS_Lead_BP = as.numeric(GWAS_Lead_BP), GWAS_P = as.numeric(GWAS_P),
    Window_Start = as.numeric(Window_Start), Window_End = as.numeric(Window_End)
  ) %>%
  filter(str_detect(Trait, regex("rhiz", ignore_case = TRUE))) %>%
  distinct(Trait, CHR, GWAS_Interval_ID, GWAS_Lead_SNP, GWAS_Lead_BP,
           GWAS_P, Window_Start, Window_End, .keep_all = TRUE) %>%
  arrange(CHR, GWAS_Lead_BP)
if (nrow(rhizome_intervals) == 0L) stop("No rhizome-associated candidate intervals were found.")

safe_interval_label <- function(x) {
  x %>% str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>% str_replace_all("^_|_$", "")
}

summarize_ld_threshold <- function(lead_ld_table, threshold) {
  threshold_hits <- lead_ld_table %>% filter(is.finite(R2_to_Lead), R2_to_Lead >= threshold)
  if (nrow(threshold_hits) == 0L) {
    return(tibble(
      Trait = unique(lead_ld_table$Trait)[1], CHR = unique(lead_ld_table$CHR)[1],
      Candidate_Interval_ID = unique(lead_ld_table$Candidate_Interval_ID)[1],
      Lead_SNP = unique(lead_ld_table$Lead_SNP)[1], Lead_BP = unique(lead_ld_table$Lead_BP)[1],
      R2_Threshold = threshold, N_Markers_In_Interval = nrow(lead_ld_table),
      N_Markers_At_Or_Above_Threshold = 0L, LD_Block_Start_BP = NA_real_,
      LD_Block_End_BP = NA_real_, LD_Block_Width_BP = NA_real_,
      Minimum_R2_Among_Selected = NA_real_, Median_R2_Among_Selected = NA_real_,
      Maximum_R2_Among_Selected = NA_real_
    ))
  }
  tibble(
    Trait = unique(lead_ld_table$Trait)[1], CHR = unique(lead_ld_table$CHR)[1],
    Candidate_Interval_ID = unique(lead_ld_table$Candidate_Interval_ID)[1],
    Lead_SNP = unique(lead_ld_table$Lead_SNP)[1], Lead_BP = unique(lead_ld_table$Lead_BP)[1],
    R2_Threshold = threshold, N_Markers_In_Interval = nrow(lead_ld_table),
    N_Markers_At_Or_Above_Threshold = nrow(threshold_hits),
    LD_Block_Start_BP = min(threshold_hits$BP, na.rm = TRUE),
    LD_Block_End_BP = max(threshold_hits$BP, na.rm = TRUE),
    LD_Block_Width_BP = max(threshold_hits$BP, na.rm = TRUE) - min(threshold_hits$BP, na.rm = TRUE),
    Minimum_R2_Among_Selected = min(threshold_hits$R2_to_Lead, na.rm = TRUE),
    Median_R2_Among_Selected = median(threshold_hits$R2_to_Lead, na.rm = TRUE),
    Maximum_R2_Among_Selected = max(threshold_hits$R2_to_Lead, na.rm = TRUE)
  )
}

ld_matrix_list <- list()
lead_ld_list <- list()
marker_count_list <- list()
summary_r2_050_list <- list()
summary_r2_070_list <- list()

for (interval_index in seq_len(nrow(rhizome_intervals))) {
  interval_row <- rhizome_intervals[interval_index, , drop = FALSE]
  interval_markers <- marker_metadata %>%
    filter(CHR == interval_row$CHR, BP >= interval_row$Window_Start, BP <= interval_row$Window_End) %>%
    arrange(BP, SNP)
  n_markers <- nrow(interval_markers)
  message("Interval ", interval_row$GWAS_Interval_ID, ": ", n_markers, " markers detected.")

  marker_count_list[[interval_row$GWAS_Interval_ID]] <- tibble(
    Trait = interval_row$Trait, CHR = interval_row$CHR,
    Candidate_Interval_ID = interval_row$GWAS_Interval_ID,
    Lead_SNP = interval_row$GWAS_Lead_SNP, Lead_BP = interval_row$GWAS_Lead_BP,
    Interval_Start_BP = interval_row$Window_Start, Interval_End_BP = interval_row$Window_End,
    N_Markers_In_Interval = n_markers,
    Lead_SNP_Present = interval_row$GWAS_Lead_SNP %in% interval_markers$SNP
  )

  if (n_markers < minimum_markers_per_interval) {
    warning("Skipping interval ", interval_row$GWAS_Interval_ID,
            " because fewer than ", minimum_markers_per_interval, " markers were present.")
    next
  }
  if (!interval_row$GWAS_Lead_SNP %in% interval_markers$SNP) {
    stop("Lead SNP ", interval_row$GWAS_Lead_SNP, " was not found in interval ",
         interval_row$GWAS_Interval_ID, ".")
  }

  genotype_interval <- geno_mat[, interval_markers$SNP, drop = FALSE]
  marker_variance <- apply(genotype_interval, 2, var, na.rm = TRUE)
  variable_markers <- is.finite(marker_variance) & marker_variance > 0
  if (!all(variable_markers)) {
    warning("Removing ", sum(!variable_markers), " zero-variance markers from interval ",
            interval_row$GWAS_Interval_ID, ".")
    genotype_interval <- genotype_interval[, variable_markers, drop = FALSE]
    interval_markers <- interval_markers[variable_markers, , drop = FALSE]
  }
  if (ncol(genotype_interval) < minimum_markers_per_interval) {
    warning("Skipping interval ", interval_row$GWAS_Interval_ID,
            " after zero-variance filtering.")
    next
  }
  if (!interval_row$GWAS_Lead_SNP %in% colnames(genotype_interval)) {
    stop("The lead SNP became unavailable after marker filtering for ",
         interval_row$GWAS_Interval_ID, ".")
  }

  correlation_matrix <- cor(genotype_interval, use = "pairwise.complete.obs", method = "pearson")
  ld_matrix <- correlation_matrix^2

  ld_matrix_list[[interval_row$GWAS_Interval_ID]] <- list(
    Trait = interval_row$Trait, CHR = interval_row$CHR,
    Candidate_Interval_ID = interval_row$GWAS_Interval_ID,
    Lead_SNP = interval_row$GWAS_Lead_SNP, Lead_BP = interval_row$GWAS_Lead_BP,
    Marker_Metadata = interval_markers, LD_R2_Matrix = ld_matrix
  )

  lead_snp <- interval_row$GWAS_Lead_SNP
  lead_r2 <- as.numeric(ld_matrix[lead_snp, colnames(ld_matrix)])
  marker_bp <- interval_markers$BP[match(colnames(ld_matrix), interval_markers$SNP)]

  lead_ld_table <- tibble(
    Trait = interval_row$Trait, CHR = interval_row$CHR,
    Candidate_Interval_ID = interval_row$GWAS_Interval_ID,
    Lead_SNP = lead_snp, Lead_BP = interval_row$GWAS_Lead_BP,
    SNP = colnames(ld_matrix), BP = marker_bp,
    Distance_From_Lead_BP = marker_bp - interval_row$GWAS_Lead_BP,
    Absolute_Distance_From_Lead_BP = abs(marker_bp - interval_row$GWAS_Lead_BP),
    R2_to_Lead = lead_r2,
    At_Or_Above_R2_0_50 = lead_r2 >= primary_r2_threshold,
    At_Or_Above_R2_0_70 = lead_r2 >= stringent_r2_threshold
  ) %>% arrange(BP, SNP)

  lead_ld_list[[interval_row$GWAS_Interval_ID]] <- lead_ld_table
  summary_r2_050_list[[interval_row$GWAS_Interval_ID]] <-
    summarize_ld_threshold(lead_ld_table, primary_r2_threshold)
  summary_r2_070_list[[interval_row$GWAS_Interval_ID]] <-
    summarize_ld_threshold(lead_ld_table, stringent_r2_threshold)

  ld_long <- as.data.frame(as.table(ld_matrix), stringsAsFactors = FALSE)
  colnames(ld_long) <- c("SNP1", "SNP2", "R2")
  marker_order <- interval_markers %>% arrange(BP, SNP) %>% pull(SNP)
  ld_long <- ld_long %>%
    mutate(
      SNP1 = factor(SNP1, levels = marker_order),
      SNP2 = factor(SNP2, levels = rev(marker_order))
    )

  heatmap_plot <- ggplot(ld_long, aes(x = SNP1, y = SNP2, fill = R2)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1), name = expression(r^2)) +
    labs(
      title = paste0("Pairwise LD within ", interval_row$GWAS_Interval_ID),
      subtitle = paste0("Lead SNP: ", interval_row$GWAS_Lead_SNP),
      x = "Marker", y = "Marker"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6)
    )

  safe_label <- safe_interval_label(interval_row$GWAS_Interval_ID)
  ggsave(
    filename = file.path(out_figures, paste0("Fig_LD_heatmap_", safe_label, ".png")),
    plot = heatmap_plot, width = 8, height = 7, dpi = 300
  )
}

interval_marker_counts <- bind_rows(marker_count_list)
write_csv(interval_marker_counts, file.path(out_tables, "LD_interval_marker_counts.csv"))
if (length(lead_ld_list) == 0L) stop("No candidate intervals contained sufficient markers for LD analysis.")

lead_ld_all <- bind_rows(lead_ld_list)
ld_summary_r2_050 <- bind_rows(summary_r2_050_list)
ld_summary_r2_070 <- bind_rows(summary_r2_070_list)

write_csv(lead_ld_all, file.path(out_tables, "LD_marker_to_lead_all_intervals.csv"))
write_csv(ld_summary_r2_050, file.path(out_tables, "LD_summary_r2_0.50_rhizome_intervals.csv"))
write_csv(ld_summary_r2_070, file.path(out_tables, "LD_summary_r2_0.70_rhizome_intervals.csv"))
saveRDS(ld_matrix_list, file.path(out_intermediate, "LD_matrices_rhizome_intervals.rds"))

lead_profile_plot <- ggplot(lead_ld_all, aes(x = BP / 1e6, y = R2_to_Lead)) +
  geom_point(size = 1.8, alpha = 0.75) +
  geom_line(linewidth = 0.35, alpha = 0.70) +
  geom_hline(yintercept = primary_r2_threshold, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = stringent_r2_threshold, linetype = "dotted", linewidth = 0.5) +
  geom_vline(aes(xintercept = Lead_BP / 1e6), linetype = "dotdash", linewidth = 0.45) +
  facet_wrap(~Candidate_Interval_ID, scales = "free_x", ncol = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Linkage disequilibrium with rhizome GWAS lead SNPs",
    subtitle = "Dashed and dotted reference lines indicate r² = 0.50 and r² = 0.70, respectively",
    x = "Physical position (Mb)", y = expression(r^2 ~ "to lead SNP")
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))

ggsave(
  filename = file.path(out_figures, "Fig_LD_to_lead_rhizome_intervals.png"),
  plot = lead_profile_plot, width = 8, height = 7, dpi = 300
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_intermediate, "sessionInfo_LD_analysis.txt")
)

message("09_LD_analysis.R completed successfully.")
message("Rhizome candidate intervals requested: ", nrow(rhizome_intervals))
message("Intervals successfully analyzed: ", length(ld_matrix_list))
message("Lead-marker LD table saved to: ",
        file.path(out_tables, "LD_marker_to_lead_all_intervals.csv"))
message("LD summaries saved in: ", out_tables)
message("LD figures saved in: ", out_figures)
message("LD matrices saved to: ",
        file.path(out_intermediate, "LD_matrices_rhizome_intervals.rds"))
