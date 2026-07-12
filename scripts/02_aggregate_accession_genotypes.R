############################################################
##
## S3011 Perennial Sorghum GWAS Pipeline
##
## Script 02: Accession-Level Genotype Aggregation
##
## Purpose:
##   - Load plant-level genotypes retained after quality control
##   - Match genotyped plants to their accession identifiers
##   - Calculate mean allele dosage for each accession and SNP
##   - Impute remaining missing values using SNP-wise means
##   - Remove SNPs with zero or undefined variance
##   - Save accession-level genotype data and marker metadata
##
## Manuscript:
##   Nabukalu et al.
##   Genetic Architecture of Rhizome Development in
##   Diploid Perennial Sorghum Reveals Candidate Regions
##   Supported by Local Ancestry
##
## Required preceding scripts:
##   00_setup.R
##   01_genotype_QC.R
##
############################################################

############################################################
## 1. Load quality-filtered inputs
############################################################

gt_filt <- readRDS(
  file.path(
    out_intermediate,
    "gt_filt_markerPlantQC.rds"
  )
)

pheno <- readRDS(
  file.path(
    out_intermediate,
    "pheno_raw.rds"
  )
)

############################################################
## 2. Validate required columns
############################################################

required_pheno_columns <- c(
  "index2",
  "Accession"
)

missing_pheno_columns <- setdiff(
  required_pheno_columns,
  colnames(pheno)
)

if (length(missing_pheno_columns) > 0L) {
  stop(
    "Missing required phenotype columns: ",
    paste(
      missing_pheno_columns,
      collapse = ", "
    )
  )
}

if (ncol(gt_filt) < 4L) {
  stop(
    "The genotype table must contain at least three ",
    "marker-metadata columns followed by plant genotypes."
  )
}

############################################################
## 3. Define marker metadata and SNP identifiers
############################################################

## The first three columns contain marker metadata.
marker_metadata <- gt_filt[, 1:3, drop = FALSE]

## Prefer an existing marker-ID column when available.
possible_id_columns <- c(
  "SNP",
  "Marker",
  "ID",
  "rs",
  "Name"
)

id_column <- intersect(
  possible_id_columns,
  colnames(marker_metadata)
)

if (length(id_column) > 0L) {
  
  snp_ids <- as.character(
    marker_metadata[[id_column[1]]]
  )
  
} else if (
  all(
    c("X.CHROM", "POS") %in%
    colnames(marker_metadata)
  )
) {
  
  snp_ids <- paste0(
    "S",
    marker_metadata$X.CHROM,
    "_",
    marker_metadata$POS
  )
  
} else {
  
  snp_ids <- paste0(
    "SNP_",
    seq_len(nrow(marker_metadata))
  )
  
}

## Ensure SNP identifiers are unique.
snp_ids <- make.unique(snp_ids)

marker_metadata$SNP <- snp_ids

############################################################
## 4. Build the plant-level genotype matrix
############################################################

## Genotype columns begin in column 4.
plant_genotype_columns <- 4:ncol(gt_filt)

## Transpose so that:
##   rows    = individual plants
##   columns = SNP markers
geno_mat_plants <- t(
  as.matrix(
    gt_filt[, plant_genotype_columns, drop = FALSE]
  )
)

storage.mode(geno_mat_plants) <- "numeric"

plant_ids <- rownames(geno_mat_plants)

if (is.null(plant_ids) || any(plant_ids == "")) {
  stop(
    "Plant identifiers are missing from genotype-column names."
  )
}

if (anyDuplicated(plant_ids) > 0L) {
  stop(
    "Duplicated plant identifiers were detected ",
    "in the genotype matrix."
  )
}

colnames(geno_mat_plants) <- snp_ids

############################################################
## 5. Construct and validate the plant-to-accession map
############################################################

map_plant_acc <- pheno[
  ,
  required_pheno_columns,
  drop = FALSE
]

colnames(map_plant_acc) <- c(
  "PlantID",
  "Accession"
)

map_plant_acc$PlantID <- as.character(
  map_plant_acc$PlantID
)

map_plant_acc$Accession <- as.character(
  map_plant_acc$Accession
)

## Remove incomplete mappings.
map_plant_acc <- map_plant_acc[
  !is.na(map_plant_acc$PlantID) &
    map_plant_acc$PlantID != "" &
    !is.na(map_plant_acc$Accession) &
    map_plant_acc$Accession != "",
  ,
  drop = FALSE
]

## Remove exact duplicate mappings.
map_plant_acc <- unique(map_plant_acc)

## Check whether any plant is assigned to multiple accessions.
mapping_counts <- aggregate(
  Accession ~ PlantID,
  data = map_plant_acc,
  FUN = function(x) length(unique(x))
)

conflicting_plants <- mapping_counts$PlantID[
  mapping_counts$Accession > 1L
]

if (length(conflicting_plants) > 0L) {
  stop(
    "The following plants map to more than one accession: ",
    paste(
      conflicting_plants,
      collapse = ", "
    )
  )
}

## Retain mappings only for plants present in the genotype matrix.
map_plant_acc <- map_plant_acc[
  map_plant_acc$PlantID %in% plant_ids,
  ,
  drop = FALSE
]

plants_without_accession <- setdiff(
  plant_ids,
  map_plant_acc$PlantID
)

if (length(plants_without_accession) > 0L) {
  warning(
    length(plants_without_accession),
    " genotyped plants lacked an accession mapping ",
    "and were excluded from accession-level aggregation."
  )
}

## Restrict the genotype matrix to mapped plants.
mapped_plant_ids <- intersect(
  plant_ids,
  map_plant_acc$PlantID
)

geno_mat_plants <- geno_mat_plants[
  mapped_plant_ids,
  ,
  drop = FALSE
]

map_plant_acc <- map_plant_acc[
  match(
    mapped_plant_ids,
    map_plant_acc$PlantID
  ),
  ,
  drop = FALSE
]

if (!identical(
  rownames(geno_mat_plants),
  map_plant_acc$PlantID
)) {
  stop(
    "Plant order differs between the genotype matrix ",
    "and the accession mapping table."
  )
}

############################################################
## 6. Aggregate plant genotypes to accession-level means
############################################################

accession_ids <- sort(
  unique(map_plant_acc$Accession)
)

geno_mat_accession <- matrix(
  NA_real_,
  nrow = length(accession_ids),
  ncol = ncol(geno_mat_plants),
  dimnames = list(
    accession_ids,
    colnames(geno_mat_plants)
  )
)

for (accession_id in accession_ids) {
  
  accession_plants <- map_plant_acc$PlantID[
    map_plant_acc$Accession == accession_id
  ]
  
  accession_dosage <- colMeans(
    geno_mat_plants[
      accession_plants,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  ## colMeans() returns NaN when all plants have missing
  ## calls for a given SNP.
  accession_dosage[
    is.nan(accession_dosage)
  ] <- NA_real_
  
  geno_mat_accession[
    accession_id,
  ] <- accession_dosage
}

############################################################
## 7. Summarize genotyped plants per accession
############################################################

plants_per_accession <- as.data.frame(
  table(map_plant_acc$Accession),
  stringsAsFactors = FALSE
)

colnames(plants_per_accession) <- c(
  "Accession",
  "Number_genotyped_plants"
)

plants_per_accession <- plants_per_accession[
  order(plants_per_accession$Accession),
  ,
  drop = FALSE
]

write.csv(
  plants_per_accession,
  file.path(
    out_tables,
    "genotyped_plants_per_accession.csv"
  ),
  row.names = FALSE
)

############################################################
## 8. Impute remaining missing accession-level values
############################################################

geno_mat_imputed <- geno_mat_accession

for (j in seq_len(ncol(geno_mat_imputed))) {
  
  missing_index <- is.na(
    geno_mat_imputed[, j]
  )
  
  if (any(missing_index)) {
    
    snp_mean <- mean(
      geno_mat_imputed[, j],
      na.rm = TRUE
    )
    
    if (is.finite(snp_mean)) {
      geno_mat_imputed[
        missing_index,
        j
      ] <- snp_mean
    }
  }
}

############################################################
## 9. Remove zero-variance or undefined SNPs
############################################################

snp_variance <- apply(
  geno_mat_imputed,
  2,
  var,
  na.rm = TRUE
)

retain_snp <- is.finite(snp_variance) &
  snp_variance > 0

n_removed <- sum(!retain_snp)

if (n_removed > 0L) {
  message(
    "Removing ",
    n_removed,
    " SNPs with zero or undefined accession-level variance."
  )
}

geno_mat_imputed <- geno_mat_imputed[
  ,
  retain_snp,
  drop = FALSE
]

marker_metadata_filtered <- marker_metadata[
  match(
    colnames(geno_mat_imputed),
    marker_metadata$SNP
  ),
  ,
  drop = FALSE
]

if (anyNA(marker_metadata_filtered$SNP)) {
  stop(
    "Marker metadata could not be matched to all retained SNPs."
  )
}

############################################################
## 10. Save accession-level outputs
############################################################

saveRDS(
  geno_mat_imputed,
  file.path(
    out_intermediate,
    "geno_mat_imp_accessionLevel.rds"
  )
)

saveRDS(
  marker_metadata_filtered,
  file.path(
    out_intermediate,
    "marker_metadata_accessionLevel.rds"
  )
)

write.csv(
  geno_mat_imputed,
  file.path(
    out_tables,
    "accession_level_mean_allele_dosage.csv"
  ),
  row.names = TRUE
)

write.csv(
  marker_metadata_filtered,
  file.path(
    out_tables,
    "retained_marker_metadata.csv"
  ),
  row.names = FALSE
)

############################################################
## 11. Save aggregation summary
############################################################

aggregation_summary <- data.frame(
  Metric = c(
    "Mapped genotyped plants",
    "Accessions retained",
    "SNPs before accession-level variance filtering",
    "SNPs after accession-level variance filtering",
    "SNPs removed for zero or undefined variance"
  ),
  Value = c(
    nrow(geno_mat_plants),
    nrow(geno_mat_imputed),
    ncol(geno_mat_accession),
    ncol(geno_mat_imputed),
    n_removed
  )
)

write.csv(
  aggregation_summary,
  file.path(
    out_tables,
    "accession_genotype_aggregation_summary.csv"
  ),
  row.names = FALSE
)

############################################################
## 12. Diagnostics
############################################################

message(
  "02_aggregate_accession_genotypes.R completed successfully."
)

message(
  "Mapped genotyped plants: ",
  nrow(geno_mat_plants)
)

message(
  "Accessions retained: ",
  nrow(geno_mat_imputed)
)

message(
  "SNPs retained: ",
  ncol(geno_mat_imputed)
)