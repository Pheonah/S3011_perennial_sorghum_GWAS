# S3011 Perennial Sorghum GWAS and Local Ancestry Analysis

## Overview

This repository contains the data and reproducible analysis pipeline
supporting the manuscript:

**Genetic Architecture of Rhizome Development in Diploid Perennial
Sorghum Reveals Candidate Regions Supported by Local Ancestry**

The study investigates the genetic architecture of
perenniality-associated traits in the diploid S3011 sorghum population
developed from an interspecific cross between *Sorghum bicolor* and
*Sorghum halepense*. Analyses include accession-level mixed-model
phenotyping, exploratory GWAS, complementary marker-based genome scans,
local ancestry analysis, linkage disequilibrium characterization,
multilocus dosage-profile analysis, candidate gene annotation, and
comparison with historical perenniality QTL.

## Repository Structure

-   scripts/
-   data/raw/
-   data/processed/
-   metadata/
-   results/figures/
-   results/tables/
-   results/intermediate/
-   results/manuscript_figures/
-   results/reproducibility/

## Experimental Population

The S3011 population was developed from a cross between KS105A (*Sorghum
bicolor*) and Gypsum 9N (*Sorghum halepense*). Diploid F4 families
derived through triploid intermediates were evaluated for
perenniality-associated traits.

## Phenotypic Traits

-   Rhizome number
-   Overwinter regrowth
-   Plant height
-   Lateral branching

Accession-level BLUPs estimated using REML-based linear mixed models
were used in downstream analyses.

## Analysis Workflow

1.  Genotype quality control
2.  Accession genotype aggregation
3.  BLUP estimation
4.  GWAS
5.  Marker-based genome scan
6.  GWAS/marker overlap
7.  Candidate gene annotation
8.  Local ancestry analysis
9.  LD analysis
10. Haplotype analysis
11. Historical QTL comparison
12. Manuscript figure generation
13. Reproducibility metadata

## Reproducibility

Running the scripts sequentially reproduces the analyses, figures, and
tables. Session information, package versions, checksums, and Git
metadata are generated automatically.

## Citation

Nabukalu P. et al. *Genetic Architecture of Rhizome Development in
Diploid Perennial Sorghum Reveals Candidate Regions Supported by Local
Ancestry*. BMC Plant Biology (under review).

## Contact

Dr. Pheonah Nabukalu\
The Land Institute\
nabukalu@landinstitute.org
