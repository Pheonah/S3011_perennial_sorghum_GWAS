############################################################
##
## S3011 Perennial Sorghum GWAS Pipeline
##
## Script 01: Plant-Level Genotype Quality Control
##
## Purpose:
##   - Import genotype and phenotype datasets
##   - Remove individuals with excessive missing genotype calls
##   - Remove low-quality SNP markers
##   - Save quality-filtered datasets for downstream analyses
##
## Manuscript:
##   Nabukalu et al.
##   Genetic Architecture of Rhizome Development in
##   Diploid Perennial Sorghum Reveals Candidate Regions
##   Supported by Local Ancestry
##
############################################################

############################################################
## Read input data
############################################################

gt <- read.csv(geno_file,
               header = TRUE,
               na.strings = "NaN")

pheno <- read.csv(pheno_file,
                  header = TRUE)

############################################################
## Plant-level genotype quality control
############################################################

## Columns 4:n contain individual plant genotype calls

geno_cols <- 4:ncol(gt)

## Remove plants with >50% missing genotype calls

na_per_ind <- apply(
  gt[, geno_cols],
  2,
  function(x) sum(is.na(x))
)

max_na <- 0.50 * nrow(gt)

keep_ind <- names(na_per_ind)[
  na_per_ind <= max_na
]

gt_ind <- gt[
  ,
  c(1:3,
    match(keep_ind,
          colnames(gt)))
]

############################################################
## Marker-level quality control
############################################################

## Criteria:
##   Call rate ≥80%
##   Minor allele frequency (MAF) ≥0.05

geno_only <- gt_ind[, 4:ncol(gt_ind)]

gtable <- data.frame(
  CHROM = gt_ind$X.CHROM,
  POS   = gt_ind$POS
)

gtable$count0 <- apply(
  geno_only,
  1,
  function(x) sum(x == 0,
                  na.rm = TRUE)
)

gtable$count1 <- apply(
  geno_only,
  1,
  function(x) sum(x == 1,
                  na.rm = TRUE)
)

gtable$count2 <- apply(
  geno_only,
  1,
  function(x) sum(x == 2,
                  na.rm = TRUE)
)

gtable$n_obs <-
  gtable$count0 +
  gtable$count1 +
  gtable$count2

n_ind <- ncol(geno_only)

gtable$MAF <-
  pmin(
    gtable$count1 +
      2 * gtable$count2,
    gtable$count1 +
      2 * gtable$count0
  ) /
  (2 * gtable$n_obs)

keep_mrk <- with(
  gtable,
  n_obs >= 0.80 * n_ind &
    MAF >= 0.05
)

gt_filt <- gt_ind[
  keep_mrk,
]

############################################################
## Save intermediate files
############################################################

saveRDS(
  gt_filt,
  file = file.path(
    out_intermediate,
    "gt_filt_markerPlantQC.rds"
  )
)

saveRDS(
  pheno,
  file = file.path(
    out_intermediate,
    "pheno_raw.rds"
  )
)

############################################################
## Save QC summary
############################################################

qc_summary <- data.frame(
  
  Plants_before = length(geno_cols),
  
  Plants_after = ncol(gt_filt) - 3,
  
  Markers_before = nrow(gt),
  
  Markers_after = nrow(gt_filt)
  
)

write.csv(
  qc_summary,
  file.path(
    out_tables,
    "QC_summary.csv"
  ),
  row.names = FALSE
)

############################################################
## Diagnostics
############################################################

message("01_genotype_QC.R completed successfully.")

message("Plants retained: ",
        ncol(gt_filt) - 3)

message("Markers retained: ",
        nrow(gt_filt))