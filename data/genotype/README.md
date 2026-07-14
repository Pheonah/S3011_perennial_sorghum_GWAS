# Genotype Data

This directory contains the genotype data used for the S3011 perennial sorghum genome-wide association analyses.

## Files

### gt_s3011_mask6.csv

Plant-level filtered genotype matrix used as the starting genotype dataset.

Rows represent SNP markers. The first columns contain marker metadata, including chromosome and physical position. The remaining columns contain genotype dosage values for individual plants.

### geno_mat_imp_accessionLevel.rds

Processed accession-level mean allele-dosage matrix generated from `gt_s3011_mask6.csv` using the genotype quality-control and accession-aggregation scripts in the `scripts/` directory.

Rows represent accessions and columns represent SNP markers.

## Processing workflow

The plant-level genotype matrix is processed using:

1. `01_genotype_QC.R`
2. `02_aggregate_accession_genotypes.R`

These scripts perform genotype quality filtering, map individual plants to accessions, calculate mean allele dosage per accession, impute residual missing values, and remove zero-variance markers.

## Notes

The genotype files correspond to the datasets used in the analyses reported in the accompanying manuscript.
