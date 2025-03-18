################################################################
# Example script for using SCLC NMF Cluster Classifier
################################################################

rm(list = ls());

# packages required
#   tidyverse packages
#     magrittr
#     tibble
#     dplyr
#     rlang
#     stringr
#     readr (for example data)
#  parsnip

# loads 1 function and 3 tibbles that are used by that function
#     pred_nmf_cluster_manucript - function to predict cluster from RNAseq data matrix
#     
#     pred_ensembl_manuscript - tibble of Ensembl IDs used by the RF model (used by function above)
#     cluster_name_manuscript - tibble of cluster numbers and names (used by function above)
#     rf_model_manuscript - RF model object (no need to access directly, is used by function above)
# change path to RData file as needed
load(file = "Output/sclc-nmf-rf-classifier-manuscript.RData");

# list of Ensembl IDs needed for RF model
pred_ensembl_manuscript

# script assumes an RNAseq data matrix or data.frame called X
# rows are samples and columns are genes
# rownames are unique sample IDs
# colnames are unique Ensembl IDs
# elements are TPM values
#     not log transformed
#     zero values are ok

# example data matrix, the NMF cluster training data
library(magrittr);
X = readr::read_csv("Data/table-s3.csv", skip = 1) %>%
  dplyr::filter(used_in_nmf == "Y") %>%
  dplyr::select(-c(gene, used_in_nmf, rank_NMF)) %>%
  tibble::column_to_rownames("ensembl_id") %>%
  t();


# for testing on subsets of rows and columns
#X = X[1:5, , drop = F]; # this is ok
#X = X[ , -2, drop = F]; # this will cause an error below due to missing gene

# first call checks X for any missing Ensembl IDs
# Ensembl ID version numbers are ignored, unless there are multiple versions of the same ID
# In that case, unneeded versions should be dropped from X first
missing_ensembl_ids = pred_nmf_cluster_manuscript(X, return_missing_ensembl_id = T, ignore_ensembl_version = T);
missing_ensembl_ids

# if any Ensembl IDs are missing, they should be added to X.

# after confirming no missing Ensembl IDs, the subtypes can be predicted from X
stopifnot(length(missing_ensembl_ids) == 0);
predicted_subtypes = pred_nmf_cluster_manuscript(X, return_missing_ensembl_id = F, ignore_ensembl_version = T);

# the subtype/cluster number, name, and probabilities are returned
predicted_subtypes

# compare cluster names and numbers
predicted_subtypes %>% dplyr::distinct(cluster, cluster_name) %>% dplyr::arrange(cluster)

# check that newly classifications match the original classifications
predicted_subtypes %>%
  dplyr::mutate(
    cluster_name_test = sample_id %>% stringr::str_split("_") %>% purrr::map_chr(dplyr::nth, 1)
  ) %>%
  dplyr::count(cluster_name, cluster_name_test) %>%
  dplyr::filter(cluster_name != cluster_name_test) %>%
  nrow() %>%
  magrittr::equals(0) %>%
  stopifnot()

# check that the cluster numbers match the column names
predicted_subtypes %>%
  dplyr::mutate(
    cluster_test = dplyr::across(matches("^cluster_prob_"), identity) %>% apply(1, which.max)
  ) %>%
  dplyr::count(cluster, cluster_test) %>%
  dplyr::filter(cluster != cluster_test) %>%
  nrow() %>%
  magrittr::equals(0) %>%
  stopifnot()


