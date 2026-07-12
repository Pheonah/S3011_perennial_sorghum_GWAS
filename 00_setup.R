############################################################
##
## S3011 Perennial Sorghum GWAS Pipeline
##
## Script 00: Project Initialization
##
## Purpose:
##   - Load required R packages
##   - Define project directories and file paths
##   - Create output folders
##   - Verify required input files
##   - Record session information for reproducibility
##
## Manuscript:
##   Nabukalu et al.
##   Genetic Architecture of Rhizome Development in
##   Diploid Perennial Sorghum Reveals Candidate Regions
##   Supported by Local Ancestry
##
############################################################

############################################################
## Required packages
############################################################

library(dplyr)
library(tidyr)
library(data.table)
library(lme4)
library(GAPIT3)
library(qqman)
library(ggplot2)

############################################################
## Project directory
############################################################

## Set the project root directory before running the pipeline.
## Example:
## root_s3011 <- "D:/VINCENT Momocs/NEWRQTLMAPPING/S3011"

setwd(root_s3011)

############################################################
## Input files
############################################################

## Define these objects before running the pipeline.
##
## geno_file  <- file.path(root_s3011,
##                         "data",
##                         "S3011_geno_GBS_filtered.csv")
##
## pheno_file <- file.path(root_s3011,
##                         "data",
##                         "S3011_rhizome_regrowth_pheno.csv")

############################################################
## Output directories
############################################################

out_results      <- file.path(root_s3011, "results")
out_figures      <- file.path(out_results, "figures")
out_tables       <- file.path(out_results, "tables")
out_intermediate <- file.path(out_results, "intermediate")

dir.create(out_results,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(out_figures,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(out_tables,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(out_intermediate,
           recursive = TRUE,
           showWarnings = FALSE)

############################################################
## Diagnostics
############################################################

message("00_setup.R: Environment initialized.")
message("Working directory: ", getwd())

message("Genotype file: ", geno_file,
        ifelse(file.exists(geno_file),
               " [FOUND]",
               " [MISSING]"))

message("Phenotype file: ", pheno_file,
        ifelse(file.exists(pheno_file),
               " [FOUND]",
               " [MISSING]"))

message("Intermediate directory: ", out_intermediate,
        ifelse(dir.exists(out_intermediate),
               " [FOUND]",
               " [MISSING]"))

############################################################
## Stop if required files are missing
############################################################

if (!file.exists(geno_file) || !file.exists(pheno_file)) {
  stop(
    "One or more required input files are missing.\n",
    "Please verify the paths to 'geno_file' and 'pheno_file'."
  )
}

############################################################
## Save session information
############################################################

writeLines(
  capture.output(sessionInfo()),
  file.path(out_results, "sessionInfo.txt")
)

message("Project initialization completed successfully.")