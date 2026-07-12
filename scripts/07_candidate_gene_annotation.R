############################################################
##
## S3011 Perennial Sorghum GWAS Pipeline
##
## Script 07: Candidate-Gene Annotation
##
## Purpose:
##   - Load rhizome-associated candidate intervals identified
##     by the exploratory GWAS
##   - Extract genes overlapping each +/-250-kb interval from
##     the Sorghum bicolor v5.1 reference annotation
##   - Join genomic coordinates to Phytozome P14 functional
##     annotations
##   - Calculate each gene's distance from the GWAS lead SNP
##   - Assign broad, annotation-based functional categories
##   - Save complete and prioritized candidate-gene tables
##
## Interpretation:
##   Genes identified by this script are positional candidates.
##   Functional categories are based on existing annotations and
##   do not establish causal relationships with rhizome development.
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
## Required analysis input:
##   results/tables/GWAS_candidate_intervals_250kb.csv
##
## Required external annotation files:
##   annotation/Sbicolor_730_v5.1.gene.gff3
##   annotation/Sbicolor_730_v5.1.P14.annotation_info.txt
##
## Outputs:
##   results/tables/
##     rhizome_candidate_intervals_used_for_annotation.csv
##     candidate_genes_per_interval_250kb_v5.1_all.csv
##     candidate_genes_per_interval_250kb_v5.1_annotated.csv
##     candidate_genes_prioritized_250kb_v5.1.csv
##     candidate_gene_functional_category_summary.csv
##
##   results/intermediate/
##     candidate_genes_annotated.rds
##     sessionInfo_candidate_gene_annotation.txt
##
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
})

############################################################
## 1. Verify required pipeline objects
############################################################

required_objects <- c(
  "root_s3011",
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
## 2. Define input paths
############################################################

candidate_interval_file <- file.path(
  out_tables,
  "GWAS_candidate_intervals_250kb.csv"
)

## Annotation files are stored in a project-level
## annotation directory by default.
annotation_dir <- file.path(
  root_s3011,
  "annotation"
)

gff_file <- file.path(
  annotation_dir,
  "Sbicolor_730_v5.1.gene.gff3"
)

functional_annotation_file <- file.path(
  annotation_dir,
  "Sbicolor_730_v5.1.P14.annotation_info.txt"
)

required_files <- c(
  candidate_interval_file,
  gff_file,
  functional_annotation_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "The following required files were not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    ),
    "\n\nPlace the Sorghum v5.1 annotation files in:\n",
    annotation_dir
  )
}

############################################################
## 3. Load and select rhizome-associated candidate intervals
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

candidate_intervals <- candidate_intervals %>%
  mutate(
    Trait = as.character(Trait),
    CHR = as.integer(CHR),
    GWAS_Lead_SNP = as.character(GWAS_Lead_SNP),
    GWAS_Lead_BP = as.numeric(GWAS_Lead_BP),
    GWAS_P = as.numeric(GWAS_P),
    Window_Start = as.numeric(Window_Start),
    Window_End = as.numeric(Window_End)
  )

## Select the rhizome-number trait while allowing reasonable
## variations in the trait label used by preceding scripts.

rhizome_trait_pattern <- regex(
  "rhiz",
  ignore_case = TRUE
)

rhizome_intervals <- candidate_intervals %>%
  filter(
    str_detect(
      Trait,
      rhizome_trait_pattern
    )
  ) %>%
  distinct(
    Trait,
    CHR,
    GWAS_Interval_ID,
    GWAS_Lead_SNP,
    GWAS_Lead_BP,
    GWAS_P,
    Window_Start,
    Window_End,
    .keep_all = TRUE
  ) %>%
  arrange(
    CHR,
    GWAS_Lead_BP
  )

if (nrow(rhizome_intervals) == 0L) {
  stop(
    "No rhizome-associated candidate intervals were found in: ",
    candidate_interval_file
  )
}

write_csv(
  rhizome_intervals,
  file.path(
    out_tables,
    "rhizome_candidate_intervals_used_for_annotation.csv"
  )
)

message(
  "Rhizome-associated candidate intervals loaded: ",
  nrow(rhizome_intervals)
)

############################################################
## 4. Load Sorghum bicolor v5.1 gene coordinates
############################################################

gff <- read_tsv(
  gff_file,
  comment = "#",
  col_names = FALSE,
  show_col_types = FALSE,
  progress = FALSE
)

if (ncol(gff) < 9L) {
  stop(
    "The GFF3 file contains fewer than nine columns: ",
    gff_file
  )
}

gff <- gff[, 1:9]

colnames(gff) <- c(
  "Reference_Chromosome",
  "Source",
  "Feature_Type",
  "Gene_Start",
  "Gene_End",
  "Score",
  "Strand",
  "Phase",
  "Attributes"
)

############################################################
## 5. Helper functions for GFF3 parsing
############################################################

extract_gff_attribute <- function(
    attributes,
    attribute_name
) {
  
  pattern <- paste0(
    "(?:^|;)",
    attribute_name,
    "=([^;]+)"
  )
  
  extracted <- str_match(
    attributes,
    pattern
  )[, 2]
  
  extracted <- str_replace_all(
    extracted,
    "%20",
    " "
  )
  
  extracted
}

standardize_chromosome <- function(x) {
  
  x_character <- as.character(x)
  
  chromosome_number <- str_extract(
    x_character,
    "[0-9]+"
  )
  
  suppressWarnings(
    as.integer(chromosome_number)
  )
}

############################################################
## 6. Extract gene records from the GFF3 file
############################################################

genes <- gff %>%
  filter(
    Feature_Type == "gene"
  ) %>%
  mutate(
    CHR = standardize_chromosome(
      Reference_Chromosome
    ),
    Gene_Start = as.numeric(
      Gene_Start
    ),
    Gene_End = as.numeric(
      Gene_End
    ),
    Gene_ID = extract_gff_attribute(
      Attributes,
      "ID"
    ),
    GFF_Name = extract_gff_attribute(
      Attributes,
      "Name"
    ),
    GFF_Description = extract_gff_attribute(
      Attributes,
      "description"
    )
  ) %>%
  filter(
    !is.na(CHR),
    is.finite(Gene_Start),
    is.finite(Gene_End),
    !is.na(Gene_ID)
  ) %>%
  select(
    CHR,
    Reference_Chromosome,
    Gene_Start,
    Gene_End,
    Strand,
    Gene_ID,
    GFF_Name,
    GFF_Description
  ) %>%
  distinct()

if (nrow(genes) == 0L) {
  stop(
    "No valid gene records could be extracted from the GFF3 file."
  )
}

############################################################
## 7. Load Phytozome P14 functional annotations
############################################################

functional_annotation <- read_tsv(
  functional_annotation_file,
  show_col_types = FALSE,
  progress = FALSE
)

if (!"locusName" %in% colnames(functional_annotation)) {
  stop(
    "The P14 annotation file must contain a 'locusName' column."
  )
}

functional_annotation <- functional_annotation %>%
  mutate(
    Gene_Core = str_extract(
      locusName,
      "Sobic\\.[0-9]{3}G[0-9]{6}"
    )
  )

if (
  all(
    is.na(
      functional_annotation$Gene_Core
    )
  )
) {
  stop(
    "No canonical Sobic gene identifiers could be parsed ",
    "from the P14 locusName column."
  )
}

############################################################
## 8. Ensure optional annotation fields exist
############################################################

optional_annotation_columns <- c(
  "Pfam",
  "GO",
  "Best-hit-arabi-defline",
  "Best-hit-rice-defline"
)

for (
  optional_column in
  optional_annotation_columns
) {
  
  if (
    !optional_column %in%
    colnames(functional_annotation)
  ) {
    
    functional_annotation[[
      optional_column
    ]] <- NA_character_
    
    warning(
      "Optional annotation column was not found and ",
      "was filled with NA: ",
      optional_column
    )
  }
}

############################################################
## 9. Extract genes overlapping each candidate interval
############################################################

extract_interval_genes <- function(
    interval_row,
    gene_table
) {
  
  chromosome_id <- as.integer(
    interval_row$CHR
  )
  
  interval_start <- as.numeric(
    interval_row$Window_Start
  )
  
  interval_end <- as.numeric(
    interval_row$Window_End
  )
  
  lead_position <- as.numeric(
    interval_row$GWAS_Lead_BP
  )
  
  interval_genes <- gene_table %>%
    filter(
      CHR == chromosome_id,
      Gene_Start <= interval_end,
      Gene_End >= interval_start
    )
  
  if (nrow(interval_genes) == 0L) {
    warning(
      "No genes were found in interval ",
      interval_row$GWAS_Interval_ID
    )
    
    return(
      tibble()
    )
  }
  
  interval_genes %>%
    mutate(
      Trait = interval_row$Trait,
      Candidate_Interval_ID =
        interval_row$GWAS_Interval_ID,
      GWAS_Lead_SNP =
        interval_row$GWAS_Lead_SNP,
      GWAS_Lead_BP =
        lead_position,
      GWAS_P =
        interval_row$GWAS_P,
      Interval_Start =
        interval_start,
      Interval_End =
        interval_end,
      Gene_Midpoint_BP =
        (
          Gene_Start +
            Gene_End
        ) / 2,
      Distance_Gene_Midpoint_to_Lead_BP =
        abs(
          Gene_Midpoint_BP -
            lead_position
        ),
      Distance_Gene_Midpoint_to_Lead_kb =
        Distance_Gene_Midpoint_to_Lead_BP /
        1000,
      Lead_SNP_Within_Gene = (
        lead_position >= Gene_Start &
          lead_position <= Gene_End
      )
    ) %>%
    select(
      Trait,
      CHR,
      Candidate_Interval_ID,
      GWAS_Lead_SNP,
      GWAS_Lead_BP,
      GWAS_P,
      Interval_Start,
      Interval_End,
      Gene_ID,
      Gene_Start,
      Gene_End,
      Gene_Midpoint_BP,
      Strand,
      Distance_Gene_Midpoint_to_Lead_BP,
      Distance_Gene_Midpoint_to_Lead_kb,
      Lead_SNP_Within_Gene,
      GFF_Name,
      GFF_Description,
      Reference_Chromosome
    )
}

candidate_genes <- map_dfr(
  seq_len(
    nrow(
      rhizome_intervals
    )
  ),
  function(interval_index) {
    
    extract_interval_genes(
      interval_row =
        rhizome_intervals[
          interval_index,
          ,
          drop = FALSE
        ],
      gene_table =
        genes
    )
  }
)

if (nrow(candidate_genes) == 0L) {
  stop(
    "No genes overlapped the rhizome-associated candidate intervals."
  )
}

############################################################
## 10. Standardize gene identifiers and join annotations
############################################################

candidate_genes <- candidate_genes %>%
  mutate(
    Gene_Core = str_extract(
      Gene_ID,
      "Sobic\\.[0-9]{3}G[0-9]{6}"
    )
  )

candidate_genes_annotated <- candidate_genes %>%
  left_join(
    functional_annotation,
    by = "Gene_Core",
    suffix = c(
      "",
      "_P14"
    )
  )

############################################################
## 11. Build a combined annotation text field
############################################################

annotation_text <- function(...) {
  
  fields <- list(...)
  
  fields <- lapply(
    fields,
    function(x) {
      
      x <- as.character(x)
      
      x[
        is.na(x)
      ] <- ""
      
      x
    }
  )
  
  do.call(
    paste,
    c(
      fields,
      sep = " | "
    )
  )
}

candidate_genes_annotated <-
  candidate_genes_annotated %>%
  mutate(
    Combined_Annotation_Text =
      annotation_text(
        GFF_Description,
        Pfam,
        GO,
        `Best-hit-arabi-defline`,
        `Best-hit-rice-defline`
      )
  )

############################################################
## 12. Assign broad annotation-based functional categories
############################################################

candidate_genes_classified <-
  candidate_genes_annotated %>%
  mutate(
    Functional_Category = case_when(
      
      str_detect(
        Combined_Annotation_Text,
        regex(
          "\\bTCP\\b|TCP transcription factor",
          ignore_case = TRUE
        )
      ) ~ "TCP transcription factor",
      
      str_detect(
        Combined_Annotation_Text,
        regex(
          "\\bAP2\\b|\\bERF\\b|ethylene response factor",
          ignore_case = TRUE
        )
      ) ~ "AP2/ERF transcription factor",
      
      str_detect(
        Combined_Annotation_Text,
        regex(
          "\\bWRKY\\b",
          ignore_case = TRUE
        )
      ) ~ "WRKY transcription factor",
      
      str_detect(
        Combined_Annotation_Text,
        regex(
          "cytokinin|response regulator|histidine phosphotransfer",
          ignore_case = TRUE
        )
      ) ~ "Cytokinin signaling",
      
      str_detect(
        Combined_Annotation_Text,
        regex(
          "auxin|AUX/IAA|auxin response factor|\\bARF\\b",
          ignore_case = TRUE
        )
      ) ~ "Auxin signaling",
      
      str_detect(
        Combined_Annotation_Text,
        regex(
          "abscisic acid|\\bABA\\b",
          ignore_case = TRUE
        )
      ) ~ "Abscisic-acid signaling",
      
      str_detect(
        Combined_Annotation_Text,
        regex(
          "meristem|axillary bud|bud outgrowth|tiller|branching",
          ignore_case = TRUE
        )
      ) ~ "Meristem or branching regulation",
      
      str_detect(
        Combined_Annotation_Text,
        regex(
          "starch|sucrose|carbohydrate|sugar transporter|glycolysis",
          ignore_case = TRUE
        )
      ) ~ "Carbohydrate metabolism or transport",
      
      str_detect(
        Combined_Annotation_Text,
        regex(
          "root|rhizome|underground stem",
          ignore_case = TRUE
        )
      ) ~ "Root or underground-organ development",
      
      str_detect(
        Combined_Annotation_Text,
        regex(
          "cold|freezing|drought|stress response|oxidative stress",
          ignore_case = TRUE
        )
      ) ~ "Abiotic-stress response",
      
      TRUE ~ "Other or insufficiently characterized"
    ),
    
    Functionally_Prioritized = (
      Functional_Category !=
        "Other or insufficiently characterized"
    )
  )

############################################################
## 13. Calculate an annotation-priority rank
############################################################

candidate_genes_classified <-
  candidate_genes_classified %>%
  mutate(
    Annotation_Priority_Score =
      as.integer(
        Functionally_Prioritized
      ) +
      as.integer(
        Lead_SNP_Within_Gene
      ) +
      as.integer(
        Distance_Gene_Midpoint_to_Lead_kb <= 50
      ),
    
    Annotation_Priority_Class =
      case_when(
        Annotation_Priority_Score >= 3 ~
          "Highest positional-functional priority",
        
        Annotation_Priority_Score == 2 ~
          "High positional-functional priority",
        
        Annotation_Priority_Score == 1 ~
          "Moderate positional or functional priority",
        
        TRUE ~
          "Unprioritized positional candidate"
      )
  )

############################################################
## 14. Construct complete and annotated tables
############################################################

candidate_genes_all <- candidate_genes_classified %>%
  arrange(
    CHR,
    GWAS_Lead_BP,
    Distance_Gene_Midpoint_to_Lead_BP,
    Gene_Start
  )

candidate_genes_with_annotation <-
  candidate_genes_all %>%
  filter(
    !is.na(locusName) |
      !is.na(Pfam) |
      !is.na(GO) |
      !is.na(`Best-hit-arabi-defline`) |
      !is.na(`Best-hit-rice-defline`) |
      !is.na(GFF_Description)
  )

prioritized_candidate_genes <-
  candidate_genes_all %>%
  filter(
    Functionally_Prioritized |
      Lead_SNP_Within_Gene |
      Distance_Gene_Midpoint_to_Lead_kb <= 50
  ) %>%
  arrange(
    desc(
      Annotation_Priority_Score
    ),
    CHR,
    Distance_Gene_Midpoint_to_Lead_BP
  )

############################################################
## 15. Save candidate-gene tables
############################################################

write_csv(
  candidate_genes_all,
  file.path(
    out_tables,
    "candidate_genes_per_interval_250kb_v5.1_all.csv"
  )
)

write_csv(
  candidate_genes_with_annotation,
  file.path(
    out_tables,
    "candidate_genes_per_interval_250kb_v5.1_annotated.csv"
  )
)

write_csv(
  prioritized_candidate_genes,
  file.path(
    out_tables,
    "candidate_genes_prioritized_250kb_v5.1.csv"
  )
)

saveRDS(
  candidate_genes_all,
  file.path(
    out_intermediate,
    "candidate_genes_annotated.rds"
  )
)

############################################################
## 16. Summarize functional categories by interval
############################################################

functional_category_summary <-
  candidate_genes_all %>%
  count(
    CHR,
    Candidate_Interval_ID,
    Functional_Category,
    name = "N_Genes"
  ) %>%
  arrange(
    CHR,
    Candidate_Interval_ID,
    desc(
      N_Genes
    )
  )

write_csv(
  functional_category_summary,
  file.path(
    out_tables,
    "candidate_gene_functional_category_summary.csv"
  )
)

############################################################
## 17. Save interval-level gene counts
############################################################

interval_gene_counts <-
  rhizome_intervals %>%
  select(
    Trait,
    CHR,
    GWAS_Interval_ID,
    GWAS_Lead_SNP,
    GWAS_Lead_BP,
    Window_Start,
    Window_End
  ) %>%
  left_join(
    candidate_genes_all %>%
      count(
        Candidate_Interval_ID,
        name = "N_Overlapping_Genes"
      ),
    by = c(
      "GWAS_Interval_ID" =
        "Candidate_Interval_ID"
    )
  ) %>%
  mutate(
    N_Overlapping_Genes =
      replace_na(
        N_Overlapping_Genes,
        0L
      )
  )

write_csv(
  interval_gene_counts,
  file.path(
    out_tables,
    "candidate_gene_counts_per_interval.csv"
  )
)

############################################################
## 18. Save session information
############################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    out_intermediate,
    "sessionInfo_candidate_gene_annotation.txt"
  )
)

############################################################
## 19. Completion messages
############################################################

message(
  "07_candidate_gene_annotation.R completed successfully."
)

message(
  "Rhizome-associated candidate intervals analyzed: ",
  nrow(
    rhizome_intervals
  )
)

message(
  "Total positional candidate genes identified: ",
  nrow(
    candidate_genes_all
  )
)

message(
  "Candidate genes with functional annotations: ",
  nrow(
    candidate_genes_with_annotation
  )
)

message(
  "Prioritized positional-functional candidates: ",
  nrow(
    prioritized_candidate_genes
  )
)

message(
  "Candidate-gene tables saved in: ",
  out_tables
)