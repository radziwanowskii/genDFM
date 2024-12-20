# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

singular_boot_for_irf <- function(Y, VAR, residuals) {
    .Call(`_genDFM_singular_boot_for_irf`, Y, VAR, residuals)
}

compute_irf_with_constant_vertical <- function(VAR, shock_var, horizon, shock_size = 1) {
    .Call(`_genDFM_compute_irf_with_constant_vertical`, VAR, shock_var, horizon, shock_size)
}

compute_all_irf <- function(VAR, loadings, shock_var, horizon, shock_size = 1) {
    .Call(`_genDFM_compute_all_irf`, VAR, loadings, shock_var, horizon, shock_size)
}

pca_cpp_arma <- function(X, n_components = 0L, center = TRUE, scale = FALSE) {
    .Call(`_genDFM_pca_cpp_arma`, X, n_components, center, scale)
}

compute_residuals <- function(Y, B, p) {
    .Call(`_genDFM_compute_residuals`, Y, B, p)
}

multivariate_multiple_regression_arma <- function(Y, X1, X2) {
    .Call(`_genDFM_multivariate_multiple_regression_arma`, Y, X1, X2)
}

favar_prep_arma <- function(data, n_components, n_lags, data_slow, data_direct) {
    .Call(`_genDFM_favar_prep_arma`, data, n_components, n_lags, data_slow, data_direct)
}

bootstrap_irf_arma_for_fyff <- function(boot_data, VAR, loadings, shock_var, horizon, no_boot, n_lags, shock_size = 1) {
    .Call(`_genDFM_bootstrap_irf_arma_for_fyff`, boot_data, VAR, loadings, shock_var, horizon, no_boot, n_lags, shock_size)
}

bootstrap_irf_arma_optimal <- function(boot_data, VAR, loadings, shock_var, horizon, no_boot, n_lags, shock_size = 1, direct = FALSE) {
    .Call(`_genDFM_bootstrap_irf_arma_optimal`, boot_data, VAR, loadings, shock_var, horizon, no_boot, n_lags, shock_size, direct)
}

row_quantiles <- function(data, probs) {
    .Call(`_genDFM_row_quantiles`, data, probs)
}

gdfm_forni_2000 <- function(matrix, M, q) {
    .Call(`_genDFM_gdfm_forni_2000`, matrix, M, q)
}

estimate_projection <- function(T_x0, Z, x_t) {
    .Call(`_genDFM_estimate_projection`, T_x0, Z, x_t)
}

estimate_projection_singular <- function(T_x0, T_xh, Z, x_T) {
    .Call(`_genDFM_estimate_projection_singular`, T_x0, T_xh, Z, x_T)
}

forecast_projection <- function(T_x0, T_xh, Z, x_T) {
    .Call(`_genDFM_forecast_projection`, T_x0, T_xh, Z, x_T)
}

restricted_gdfm <- function(matrix, M, q, r) {
    .Call(`_genDFM_restricted_gdfm`, matrix, M, q, r)
}

compute_sorted_eigenvalues <- function(matrix, sequence) {
    .Call(`_genDFM_compute_sorted_eigenvalues`, matrix, sequence)
}

Onatski_2009_test <- function(sorted_eigenvalues, CVs, k0, k1) {
    .Call(`_genDFM_Onatski_2009_test`, sorted_eigenvalues, CVs, k0, k1)
}

Onatski_2010_test <- function(data, rmax) {
    .Call(`_genDFM_Onatski_2010_test`, data, rmax)
}

