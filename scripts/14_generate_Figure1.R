############################################################
##
## S3011 Perennial Sorghum GWAS Pipeline
##
## Generate Figure 1: Phenotypic Distributions
##
## Panel A:
##   Rhizome-number scores for the 606 genotyped individuals
##   retained after genotype quality control.
##
## Panel B:
##   Accession-level BLUPs for overwinter regrowth and rhizome
##   number for the 124 genotyped accessions used in downstream
##   genome-wide analyses.
##
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

############################################################
## 1. Project paths
############################################################

## Edit this path only if the S3011 project is stored elsewhere.
root_s3011 <- "E:/VINCENT Momocs/NEWRQTLMAPPING/S3011"

out_results      <- file.path(root_s3011, "results")
out_figures      <- file.path(out_results, "figures")
out_intermediate <- file.path(out_results, "intermediate")

dir.create(out_figures, recursive = TRUE, showWarnings = FALSE)

############################################################
## 2. Required input files
############################################################

pheno_rds <- file.path(
  out_intermediate,
  "pheno_raw.rds"
)

plant_genotype_qc_rds <- file.path(
  out_intermediate,
  "gt_filt_markerPlantQC.rds"
)

accession_genotype_rds <- file.path(
  out_intermediate,
  "geno_mat_imp_accessionLevel.rds"
)

rhiz_regrowth_blup_rds <- file.path(
  out_intermediate,
  "blup_panel_rhiz_reg.rds"
)

required_files <- c(
  pheno_rds,
  plant_genotype_qc_rds,
  accession_genotype_rds,
  rhiz_regrowth_blup_rds
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "The following required files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

############################################################
## 3. Load data
############################################################

pheno <- readRDS(pheno_rds)
gt_plant_qc <- readRDS(plant_genotype_qc_rds)
geno_accession <- readRDS(accession_genotype_rds)
blup_panel <- readRDS(rhiz_regrowth_blup_rds)

############################################################
## 4. Validate required columns and identifiers
############################################################

required_pheno_columns <- c("index2", "Rhiz_num")
missing_pheno_columns <- setdiff(required_pheno_columns, names(pheno))

if (length(missing_pheno_columns) > 0L) {
  stop(
    "Phenotype data are missing required columns: ",
    paste(missing_pheno_columns, collapse = ", ")
  )
}

required_blup_columns <- c(
  "Accession",
  "Rhiz_BLUP",
  "Regrowth_BLUP"
)

missing_blup_columns <- setdiff(required_blup_columns, names(blup_panel))

if (length(missing_blup_columns) > 0L) {
  stop(
    "BLUP data are missing required columns: ",
    paste(missing_blup_columns, collapse = ", ")
  )
}

if (ncol(gt_plant_qc) < 4L) {
  stop(
    "The plant-level genotype QC object does not contain genotype columns."
  )
}

if (is.null(rownames(geno_accession))) {
  stop(
    "The accession-level genotype matrix must contain accession IDs as row names."
  )
}

############################################################
## 5. Prepare Panel A data: 606 genotyped individuals
############################################################

genotyped_plant_ids <- colnames(gt_plant_qc)[4:ncol(gt_plant_qc)]

pheno_606 <- pheno %>%
  mutate(
    index2 = as.character(index2),
    Rhiz_num = as.numeric(as.character(Rhiz_num))
  ) %>%
  filter(
    index2 %in% genotyped_plant_ids,
    is.finite(Rhiz_num)
  )

n_panel_a <- nrow(pheno_606)

if (n_panel_a != 606L) {
  stop(
    "Panel A should contain 606 genotyped individuals, but ",
    n_panel_a,
    " phenotype records were matched. Check plant identifiers."
  )
}

pheno_606_plot <- pheno_606 %>%
  count(
    Rhiz_num,
    name = "Frequency"
  ) %>%
  arrange(Rhiz_num)

if (sum(pheno_606_plot$Frequency) != 606L) {
  stop("Panel A frequencies do not sum to 606.")
}

############################################################
## 6. Prepare Panel B data: 124 genotyped accessions
############################################################

genotyped_accessions <- rownames(geno_accession)

if (length(genotyped_accessions) != 124L) {
  stop(
    "Expected 124 genotyped accessions, but found ",
    length(genotyped_accessions),
    "."
  )
}

blup_panel_124 <- blup_panel %>%
  mutate(
    Accession = as.character(Accession)
  ) %>%
  filter(
    Accession %in% genotyped_accessions
  ) %>%
  distinct(
    Accession,
    .keep_all = TRUE
  ) %>%
  arrange(
    match(Accession, genotyped_accessions)
  )

n_panel_b_accessions <- n_distinct(blup_panel_124$Accession)

if (n_panel_b_accessions != 124L) {
  missing_accessions <- setdiff(
    genotyped_accessions,
    blup_panel_124$Accession
  )
  
  stop(
    "Panel B should contain 124 genotyped accessions, but ",
    n_panel_b_accessions,
    " were matched.\nMissing accessions:\n",
    paste(missing_accessions, collapse = "\n")
  )
}

blup_long <- blup_panel_124 %>%
  select(
    Accession,
    Regrowth_BLUP,
    Rhiz_BLUP
  ) %>%
  pivot_longer(
    cols = c(
      Regrowth_BLUP,
      Rhiz_BLUP
    ),
    names_to = "Trait",
    values_to = "BLUP"
  ) %>%
  mutate(
    Trait = recode(
      Trait,
      Regrowth_BLUP = "Overwinter regrowth",
      Rhiz_BLUP = "Rhizome number"
    ),
    Trait = factor(
      Trait,
      levels = c(
        "Overwinter regrowth",
        "Rhizome number"
      )
    )
  ) %>%
  filter(
    is.finite(BLUP)
  )

trait_counts <- table(blup_long$Trait)

if (any(trait_counts != 124L)) {
  stop(
    "Each Panel B trait should contain 124 BLUP values.\nObserved counts:\n",
    paste(
      names(trait_counts),
      trait_counts,
      sep = ": ",
      collapse = "\n"
    )
  )
}

############################################################
## 7. Panel A: rhizome-number distribution
############################################################

panel_a <- ggplot(
  pheno_606_plot,
  aes(
    x = Rhiz_num,
    y = Frequency
  )
) +
  geom_col(
    width = 0.10,
    fill = "#5B3517",
    color = "black",
    linewidth = 0.45
  ) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1, 2, 3),
    labels = c("0", "0.5", "1", "2", "3"),
    limits = c(-0.2, 3.2)
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.08)
    )
  ) +
  labs(
    title = "A",
    x = "Rhizome number score",
    y = "Number of genotyped individuals"
  ) +
  theme_classic(
    base_size = 15
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0
    ),
    axis.title = element_text(
      size = 15
    ),
    axis.text = element_text(
      size = 13,
      color = "black"
    ),
    axis.line = element_line(
      linewidth = 0.65
    ),
    axis.ticks = element_line(
      linewidth = 0.65
    ),
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )

############################################################
## 8. Panel B: accession-level BLUP distributions
############################################################

panel_b <- ggplot(
  blup_long,
  aes(
    x = Trait,
    y = BLUP,
    fill = Trait
  )
) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    width = 0.55,
    alpha = 0.75,
    color = "black",
    linewidth = 0.45
  ) +
  geom_jitter(
    width = 0.07,
    height = 0,
    size = 1.35,
    alpha = 0.50,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Overwinter regrowth" = "#0E5448",
      "Rhizome number" = "#8A5527"
    )
  ) +
  coord_cartesian(
    ylim = c(-0.40, 0.40)
  ) +
  scale_y_continuous(
    breaks = c(-0.25, 0, 0.25)
  ) +
  labs(
    title = "B",
    x = NULL,
    y = "Accession-level BLUP"
  ) +
  guides(
    fill = "none"
  ) +
  theme_classic(
    base_size = 15
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0
    ),
    axis.title.y = element_text(
      size = 15
    ),
    axis.text.x = element_text(
      size = 11,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    axis.line = element_line(
      linewidth = 0.6
    ),
    axis.ticks = element_line(
      linewidth = 0.6
    ),
    plot.margin = margin(
      10,
      10,
      10,
      15
    )
  )

############################################################
## 9. Combine panels
############################################################

figure_1 <- panel_a +
  panel_b +
  plot_layout(
    widths = c(1, 1.15)
  )

############################################################
## 10. Save Figure 1
############################################################

png_file <- file.path(
  out_figures,
  "Fig_1_Phenotype_Distributions.png"
)

pdf_file <- file.path(
  out_figures,
  "Fig_1_Phenotype_Distributions.pdf"
)

ggsave(
  filename = png_file,
  plot = figure_1,
  width = 12,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

pdf_device <- if (capabilities("cairo")) {
  cairo_pdf
} else {
  "pdf"
}

ggsave(
  filename = pdf_file,
  plot = figure_1,
  width = 12,
  height = 6,
  units = "in",
  device = pdf_device,
  bg = "white"
)

############################################################
## 11. Save figure data and report checks
############################################################

write.csv(
  pheno_606_plot,
  file.path(
    out_figures,
    "Fig_1_PanelA_Data.csv"
  ),
  row.names = FALSE
)

write.csv(
  blup_long,
  file.path(
    out_figures,
    "Fig_1_PanelB_Data.csv"
  ),
  row.names = FALSE
)

message("Figure 1 generated successfully.")
message("Panel A individuals: ", n_panel_a)
message("Panel B accessions: ", n_panel_b_accessions)
message("PNG: ", png_file)
message("PDF: ", pdf_file)