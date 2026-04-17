############################
# 0. Load packages
############################
library(tidyverse)

############################
# 1. Input loader
############################
# Supports .rds, .csv, and .tsv input formats
# Assumes CpGs in rows and samples in columns
load_input <- function(file) {
  tolower(tools::file_ext(file)) -> ext
  
  if (ext == "rds") {
    readRDS(file) -> data
  } else if (ext == "csv") {
    read.csv(file, row.names = 1, check.names = FALSE) -> data
  } else if (ext %in% c("tsv", "txt")) {
    read.delim(file, row.names = 1, check.names = FALSE) -> data
  } else {
    stop("Unsupported file format")
  }
  
  as.data.frame(data) -> data
  
  return(data)
}

############################
# 2. cAge
############################
run_cAge <- function(input_file) {
  
  ############################
  # Load data
  ############################
  load_input(input_file) -> data
  
  ############################
  # Load model data
  ############################
  read.delim("data/elnet_coefficients_linear.tsv") -> coefficients
  read.delim("data/elnet_coefficients_log.tsv") -> coefficients_log
  read.delim("data/cpg_meanbeta_gs20k.tsv") -> means
  
  coefficients$CpG_Site -> rownames(coefficients)
  coefficients_log$CpG_Site -> rownames(coefficients_log)
  means$cpg -> rownames(means)
  
  30.7873138968084 -> intercept
  3.87249881475691 -> intercept_log
  
  ############################
  # Preprocessing
  ############################
  
  # Split linear and quadratic CpGs
  coefficients[grep("_2", rownames(coefficients)), , drop = FALSE] -> coef_2
  gsub("_2", "", rownames(coef_2)) -> coef_2_simp
  coefficients[!rownames(coefficients) %in% rownames(coef_2), , drop = FALSE] -> coef
  
  coefficients_log[grep("_2", rownames(coefficients_log)), , drop = FALSE] -> coef_log_2
  gsub("_2", "", rownames(coef_log_2)) -> coef_log_2_simp
  coefficients_log[!rownames(coefficients_log) %in% rownames(coef_log_2), , drop = FALSE] -> coef_log
  
  union(rownames(coef), rownames(coef_log)) -> cpgs_linear
  union(coef_2_simp, coef_log_2_simp) -> cpgs_squared
  union(cpgs_linear, cpgs_squared) -> all_cpgs
  
  # Ensure CpGs are rows
  if (ncol(data) > nrow(data)) {
    t(data) -> data
  }
  
  # Subset to required CpGs
  data[intersect(rownames(data), all_cpgs), ] -> coef_data
  
  # Convert M-values to beta values if needed
  if (max(coef_data, na.rm = TRUE) > 1) {
    (2^coef_data / (2^coef_data + 1)) -> coef_data
  }
  
  ############################
  # Imputation
  ############################
  
  # Add missing CpGs using mean beta values
  setdiff(all_cpgs, rownames(coef_data)) -> missing_cpgs
  
  if (length(missing_cpgs) > 0) {
    matrix(nrow = length(missing_cpgs), ncol = ncol(coef_data)) -> mat
    rownames(mat) <- missing_cpgs
    colnames(mat) <- colnames(coef_data)
    
    means[missing_cpgs, "mean"] -> mvals
    mat[,] <- mvals
    
    rbind(coef_data, mat) -> coef_data
  }
  
  # Replace remaining NA values with probe-wise mean
  apply(coef_data, 1, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  }) %>% t() -> coef_data
  
  ############################
  # Prediction
  ############################
  
  # Linear model
  coef_data[rownames(coef), ] -> scores
  coef_data[coef_2_simp, ]^2 -> scores_quad
  paste0(rownames(scores_quad), "_2") -> rownames(scores_quad)
  rbind(scores, scores_quad) -> scores_linear
  
  coefficients[rownames(scores_linear), "Coefficient"] -> coefs
  colSums(scores_linear * coefs) + intercept -> pred_linear
  
  # Log model (used for younger samples)
  coef_data[rownames(coef_log), ] -> scores_log
  coef_data[coef_log_2_simp, ]^2 -> scores_quad_log
  paste0(rownames(scores_quad_log), "_2") -> rownames(scores_quad_log)
  rbind(scores_log, scores_quad_log) -> scores_log
  
  coefficients_log[rownames(scores_log), "Coefficient"] -> coefs_log
  colSums(scores_log * coefs_log) + intercept_log -> pred_log
  exp(pred_log) -> pred_log
  
  # Combine models based on age threshold
  names(pred_linear[pred_linear > 20]) -> over20
  names(pred_linear[pred_linear <= 20]) -> under20
  
  c(pred_log[under20], pred_linear[over20]) -> predictions
  
  ############################
  # Output
  ############################
  
  data.frame(Sample = names(predictions), cAge = predictions)
}

############################
# 3. Garma
############################
run_garma <- function(input_file) {
  
  ############################
  # Load data
  ############################
  load_input(input_file) -> beta
  rownames(beta) -> beta$site
  
  ############################
  # Load model data
  ############################
  read.csv("data/Garma_linear.csv") -> Linear
  read.csv("data/Garma_quadratic.csv") -> Quadratic
  read.csv("data/Garma_mean.csv") -> mean_beta
  
  16.24957834871524 -> intercept
  
  ############################
  # Prediction
  ############################
  
  # Linear component
  Linear %>%
    left_join(mean_beta, by = "site") %>%
    left_join(beta, by = "site") %>%
    pivot_longer(-c(1:5),
                 names_to = "Sample", values_to = "beta") %>%
    mutate(beta = if_else(is.na(beta), Mean.beta, beta),
           val = beta * Coefficient) %>%
    group_by(Sample) %>%
    summarise(Linear = sum(val), .groups = "drop") -> lin
  
  # Quadratic component
  Quadratic %>%
    left_join(mean_beta, by = "site") %>%
    left_join(beta, by = "site") %>%
    pivot_longer(-c(1:5),
                 names_to = "Sample", values_to = "beta") %>%
    mutate(beta = if_else(is.na(beta), Mean.beta, beta),
           val = (beta^2) * Coefficient) %>%
    group_by(Sample) %>%
    summarise(Quadratic = sum(val), .groups = "drop") -> quad
  
  ############################
  # Output
  ############################
  
  left_join(lin, quad, by = "Sample") %>%
    mutate(garma = Linear + Quadratic + intercept) %>%
    select(Sample, garma)
}

############################
# 4. PAYA
############################
run_paya <- function(input_file) {
  
  ############################
  # Load data
  ############################
  load_input(input_file) -> beta
  
  ############################
  # Load model data
  ############################
  read.csv("data/PAYA_coeff.csv") -> coeff
  coeff$Site -> rownames(coeff)
  
  read.csv("data/site_values_mean.csv") -> imputation
  imputation$cpg -> rownames(imputation)
  imputation[, "mean", drop = FALSE] -> imputation
  
  ############################
  # Imputation
  ############################
  
  # Add missing CpGs using mean beta values
  setdiff(rownames(coeff), rownames(beta)) -> missing
  
  if (length(missing) > 0) {
    imputation[missing, , drop = FALSE] -> imp
    imp[, rep(1, ncol(beta)), drop = FALSE] -> imp_expanded
    colnames(imp_expanded) <- colnames(beta)
    rbind(beta, imp_expanded) -> beta
  }
  
  # Align CpGs to model order
  beta[rownames(coeff), , drop = FALSE] -> sub
  
  ############################
  # Prediction
  ############################
  
  -0.05872492 -> intercept
  
  as.matrix(sub) -> sub
  storage.mode(sub) <- "numeric"
  sub[is.na(sub)] <- 0
  
  colSums(sub * coeff$Coefficient) + intercept -> pred
  
  ############################
  # Inverse-transformation
  ############################
  
  invFage <- function(x, adult_age = 20) {
    ifelse(x < 0,
           (1 + adult_age) * exp(x) - 1,
           (1 + adult_age) * x + adult_age)
  }
  
  invFage(pred) -> pred
  
  ############################
  # Output
  ############################
  
  data.frame(Sample = names(pred), paya = pred)
}

############################
# 5. Wrapper
############################
# Runs all clocks and returns a combined results table
run_clocks <- function(input_file, out_file) {
  
  list(
    cAge = run_cAge(input_file),
    garma = run_garma(input_file),
    paya = run_paya(input_file)
  ) -> res_list
  
  Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE), res_list) -> final
  
  rowMeans(final[, -1], na.rm = TRUE) -> final$mean
  
  final[, c("Sample", "cAge", "garma", "paya", "mean")] -> final
  colnames(final) <- c("Sample", "cAge", "Garma", "PAYA", "Ensaya")
  
  write.csv(final, out_file, row.names = FALSE)
  
  return(final)
}