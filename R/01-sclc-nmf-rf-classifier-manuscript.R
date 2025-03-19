rm(list=ls());

###
# Load global variables used by classifier function
###

load(file = "Data/29a-rf-model-manuscript-data-v04.RData");

###
# SCLC NMF Cluster Classifier
#
# X = TPM matrix or data.frame, to be classified into one of 8 NMF clusters
#     rows are samples and columns are genes
#     rownames are unique sample IDs
#     colnames are unique Ensembl IDs
#     elements are TPM values
#         not log transformed
#         zero values are ok
#     all 8437 genes listed in global variable pred_ensembl_manuscript need to be present in X
#
# return_missing_ensembl_id = logical, if TRUE then return a data.frame of Ensembl IDs that are missing from X 
#                             rather than classifying X. Otherwise, return the cluster assignments for each sample in X.
#
# ignore_ensembl_version = logical, if TRUE then the version number in Ensembl IDs are not considered when matching genes in X.
#                          If this creates ambiguities, e.g. two Ensembl IDs that match except for the version number, then an error is thrown.
#
#
# Returns a tibble with sample IDs, cluster assignments, and cluster probabilities, one row per row of X
###

pred_nmf_cluster_manuscript = function(X, return_missing_ensembl_id = F, ignore_ensembl_version = T) {
  require(magrittr);
  
  missing_rf = !exists("rf_model_manuscript");
  missing_ensembl = !exists("pred_ensembl_manuscript");
  missing_name = !exists("cluster_name_manuscript")
  if (missing_rf || missing_ensembl || missing_name) {
    rlang::abort("The random forest model is missing, please reload 29a-rf-model-manuscript-v04.RData")
  }
  
  stopifnot(nrow(rf_model_manuscript) == 1);
  
  stopifnot(is.matrix(X) || is.data.frame(X));
  if (is.matrix(X)) {
    X %<>% as.data.frame();
  }
  stopifnot(tibble::has_rownames(X));
  
  cols = colnames(X);
  stopifnot(length(cols) > 0);
  
  x = cols %>% grepl("^ENSG[0-9]{11}(\\.[0-9]+)?$", .);
  nx = sum(!x);
  if (nx == 1) {
    rlang::abort(stringr::str_glue("1 column name is not in Ensembl ID format: {cols[!x]}"));
  } else if (nx > 1) {
    rlang::abort(stringr::str_glue("{nx} column names are not in Ensembl ID format"));
  }
  
  if (ignore_ensembl_version) {
    cols2 = cols %>% stringr::str_replace("\\.[0-9]+$", "");
    tmp = tibble::tibble(cols, cols2) %>%
      dplyr::group_by(cols2) %>%
      dplyr::filter(dplyr::n() > 1) %>%
      dplyr::ungroup();
    nx = tmp %>% nrow();
    if (nx == 1) {
      rlang::abort(stringr::str_glue("ignore_ensembl_version = T cannot be used became two versions of the same Ensembl ID are present: {tmp$cols %>% toString()}"));
    } else if (nx > 1) {
      rlang::abort(stringr::str_glue("ignore_ensembl_version = T cannot be used became {tmp %>% dplyr::distinct(cols2) %>% nrow()} Ensembl IDs have multiple versions present"));
    }
    
    cols = cols2;
    rm(cols2, tmp);
    colnames(X) = cols;
  }
  
  missing_ensembl_id = pred_ensembl_manuscript %>%
    dplyr::filter(
      !ensembl_id %in% cols,
      !ensembl_id2 %in% cols
    );
  nx = missing_ensembl_id %>% nrow();
  if (nx > 0) {
    if (return_missing_ensembl_id) {
      return(missing_ensembl_id);
    } else if (nx == 1) {
      rlang::abort(stringr::str_glue_data(missing_ensembl_id, "1 Ensembl ID needed by the model is missing: {ensembl_id} or {ensembl_id2}"));
    } else {
      rlang::abort(stringr::str_glue("{nx} Ensembl IDs needed by the model are missing; use return_missing_ensembl_id = T to get a list"));
    }
  } else if (return_missing_ensembl_id) {
    message("No Ensembl IDs are missing; use return_missing_ensembl_id = F to get the predicted clusters");
    return(c());
  }
  
  x = cols %>% grepl("^ENSG[0-9]{11}\\.[0-9]+$", .);
  nx = sum(!x);
  if (nx > 0) {
    # missing version numbers, so add them
    x_ensembl = tibble::tibble(ensembl_id = dplyr::if_else(grepl("\\.", cols), cols, NA_character_)) %>%
      dplyr::mutate(
        ensembl_id2 = dplyr::if_else(is.na(ensembl_id), cols, NA_character_)
      );
    
    x_ensembl %<>%
      dplyr::left_join(
        pred_ensembl_manuscript,
        by = c("ensembl_id2")
      ) %>%
      dplyr::mutate(
        ensembl_id.x = dplyr::if_else(is.na(ensembl_id.x), ensembl_id.y, ensembl_id.x)
      ) %>%
      dplyr::rename(ensembl_id = ensembl_id.x) %>%
      dplyr::select(-ensembl_id.y)
    
    # rename to Ensembl ID with the version number
    colnames(X) = x_ensembl$ensembl_id;
  }
  
  # RF model expects format tpm_ENSEMBLID.VERSION
  colnames(X) = X %>% colnames() %>% paste("tpm", ., sep = "_");
  
  cl = X %>%
    rownames() %>%
    tibble::tibble(sample_id = .) %>%
    dplyr::mutate(
      cluster = rf_model_manuscript$rf_mod[[1]] %>% predict(new_data = X, type = "class") %>% dplyr::pull(.pred_class)
    ) %>%
    dplyr::bind_cols(
      rf_model_manuscript$rf_mod[[1]] %>% predict(new_data = X, type = "prob")
    ) %>%
    dplyr::rename_with(~ .x %>% stringr::str_replace("^\\.pred_", "cluster_prob_")) %>%
    dplyr::left_join(
      cluster_name_manuscript,
      by = "cluster"
    ) %>%
    dplyr::relocate(cluster_name, .after = "cluster");
  
  cl %<>%
    tidyr::pivot_longer(matches("^cluster_prob_[0-9]+$"), names_prefix = "cluster_prob_") %>%
    dplyr::left_join(
      cluster_name_manuscript %>%
        dplyr::arrange(cluster_ord) %>%
        dplyr::transmute(
          name = cluster %>% as.character(),
          name2 = cluster_ord %>% factor() %>% forcats::fct_inorder()
        ),
      by = "name"
    ) %>%
    dplyr::select(-name) %>%
    dplyr::rename(name = name2) %>%
    tidyr::pivot_wider(names_sort = T, names_prefix = "cluster_prob_") %>%
    dplyr::mutate(
      cluster = cluster_ord
    ) %>%
    dplyr::select(-cluster_ord);
  
  return(cl);
}


###
# Export the classifier function and associated global variables
###

save(pred_nmf_cluster_manuscript, rf_model_manuscript, pred_ensembl_manuscript, 
     cluster_name_manuscript, file = "Output/sclc-nmf-rf-classifier-manuscript.RData")


