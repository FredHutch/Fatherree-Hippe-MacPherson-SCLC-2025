# Fatherree-Hippe-MacPherson-SCLC-2025

Analysis code used in Fatherree, Hippe, et al., 
_Transcriptional analyses of patient derived xenograft models reveal distinct subsets of small cell lung cancer_

## Scripts

1. `R/01-sclc-nmf-rf-classifier-manuscript.R`
  - Defines the classifier function and saves it with the required global variables
  - Reads in the global variables from the file `Data/29a-rf-model-manuscript-data-v04.RData`
  - Writes the file `Output/sclc-nmf-rf-classifier-manuscript.RData` with the 
  function and global variables
2. `R/02-sclc-nmf-rf-classifier-manuscript-example.R`
  - Example script that uses the original training data (from supplemental Table S3 
  in the manuscript, `Data/table-s3.csv`) to illustrate all of the steps for 
  classifying new data and the corresponding output
  - Reads in `Output/sclc-nmf-rf-classifier-manuscript.RData` as created by the 
  first script (`R/01-sclc-nmf-rf-classifier-manuscript.R`)
  - Reads in `Data/table-s3.csv` for the test data
  