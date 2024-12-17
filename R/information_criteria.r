#' Test k=k_0 vs k_0<k<=k_1
#'
#' This function performs the Onatski (2009) test to compare k=k_0 vs k_0<k<=k_1.
#'
#' @param data A matrix of data to be tested.
#' @param omega0 Optional parameter omega0. Default is NULL.
#' @param sequence Optional parameter sequence. Default is NULL.
#' @param k0 Parameter k0. Default is 1.
#' @param k1 Parameter k1. Default is 5.
#' @param test_size Size of the test. Default is 5.
#' @return The result of the Onatski (2009) test.
#' @references
#' Onatski, A. (2009). Testing hypotheses about the number of factors in large factor models. Econometrica, 77(5), 1447–1479.
#' @name onatski_2009_test_r
#' @export
onatski_2009_test_r <- function(data, omega0 = NULL, sequence = NULL, k0 = 1, k1 = 5, test_size = 5) {
  smoothed_periodogram <- NULL
  if (is.null(omega0) && is.null(sequence)) {
    omega0 = 0
    sequence = seq(4:min(floor(3.3 * sqrt(length(data[, 1]))), floor(length(data[, 1]) / 2 - 0.1)))
  } else if (is.null(sequence)) {
    m0 = floor(length(data[, 1]) * omega0 / (2 * pi))
    m1 = min(m0, floor(3.3 * sqrt(length(data[, 1]))))
    sequence = seq(m0 - m1 + 4, min(m0 - m1 + 3 + floor(3.3 * sqrt(length(data[, 1])), floor(length(data[, 1]) / 2 - 0.1))))
  }
  CV_matrix <- matrix(c(
    2.75, 3.62, 4.15, 4.54, 4.89, 5.20, 5.45, 5.70,
    3.33, 4.31, 4.91, 5.40, 5.77, 6.13, 6.42, 6.66,
    3.50, 4.49, 5.13, 5.62, 6.03, 6.39, 6.67, 6.92,
    3.69, 4.72, 5.37, 5.91, 6.31, 6.68, 6.95, 7.25,
    3.92, 4.99, 5.66, 6.24, 6.62, 7.00, 7.32, 7.59,
    4.20, 5.31, 6.03, 6.57, 7.00, 7.41, 7.74, 8.04,
    4.52, 5.73, 6.46, 7.01, 7.50, 7.95, 8.29, 8.59,
    5.02, 6.26, 6.97, 7.63, 8.16, 8.61, 9.06, 9.36,
    5.62, 6.91, 7.79, 8.48, 9.06, 9.64, 10.11, 10.44,
    6.55, 8.15, 9.06, 9.93, 10.47, 11.27, 11.75, 12.13,
    8.74, 10.52, 11.67, 12.56, 13.42, 14.26, 14.88, 15.25
  ), nrow = 11, byrow = TRUE)

  cvs <- CV_matrix[length(CV_matrix[, 1]) - test_size + 1, ]
  sorted_eigenvalues <- compute_sorted_eigenvalues(data, sequence)
  return(Onatski_2009_test(sorted_eigenvalues, cvs, 1, 7))
}

#' Determine Number of Dynamic Factors
#'
#' This function determines the number of dynamic factors in the data using specified methods.
#'
#' @param data A matrix of data to be analyzed.
#' @param method Method to determine factors. Must be one of "on", "av", or "hl". Default is "av".
#' @param max_factors Maximum number of factors. Default is 10.
#' @return The number of dynamic factors.
#' @details
#' This function uses the following methods to determine the number of dynamic factors:
#' \itemize{
#'   \item \code{"on"}: Onatski (2009) test.
#'   \item \code{"av"}: Avarucci et al. (2022) test - relies on the `fnets` package.
#'   \item \code{"hl"}: Hallin & Liska (2007) test - relies on the `fnets` package.
#' }
#' @references
#' Onatski, A. (2009). Testing hypotheses about the number of factors in large factor models. Econometrica, 77(5), 1447–1479.
#' Avarucci, M., Cavicchioli, M., Mario, F., & Zaffaroni, P. (2021, 01). The main business cycle shock(s). frequency-band estimation of the number of dynamic factors.
#' Hallin, M., & Liska, R. (2007, 02). Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102,603-617.
#' @import fnets
#' @name determine_number_of_dynamic_factors
#' @export
determine_number_of_dynamic_factors <- function(data, method = "av", max_factors = 10) {
  
  if(sum(is.na(data)) > 0){
    stop("Data contains missing values")
  }
  if(!is.matrix(data)){
    data <- as.matrix(data)
  }
  if (!is.numeric(data)) {
    stop("Data must be numeric")
  }
    
    switch(method,
                     "on" = {
                         cat("Number of factors determined by the Onatski (2009) test\n")
                         return(onatski_2009_test_r(data, k0 = 1, k1 = max_factors, test_size = 5))
                     },
                     "av" = {
                         cat("Number of factors determined by the Avarucci_et_al_2022 test\n")
                         return(fnets_ic_wrapper(fnets::factor.number(data, fm.restricted = FALSE, method = "er", q.max = max_factors)))
                     },
                     "hl" = {
                         cat("Number of factors determined by the Hallin & Liska (2007) test\n")
                         return(fnets_ic_wrapper(fnets::factor.number(data, fm.restricted = FALSE, method = "ic", q.max = max_factors)))
                     },
                      stop("Invalid method")
    )
}

#' Determine Number of Static Factors
#'
#' This function determines the number of static factors in the data using specified methods.
#'
#' @param data A matrix of data to be analyzed.
#' @param method Method to determine factors. Must be one of "on", "al", "ah", or "bn". Default is "on".
#' @param max_factors Maximum number of factors. Default is 10.
#' @return The number of static factors.
#' @details
#' This function uses the following methods to determine the number of static factors:
#' \itemize{
#'   \item \code{"on"}: Onatski (2010) test.
#'   \item \code{"al"}: Alessi, Barigozzi and Capasso (2010) test - relies on the `fnets` package.
#'   \item \code{"ah"}: Ahn and Horenstein (2013) test - relies on the `fnets` package.
#'   \item \code{"bn"}: Bai and Ng (2002) test - relies on the `dfms` package.
#' }
#' @references
#' Onatski, A. (2010). Determining the number of factors from empirical distribution of eigenvalues. The Review of Economics and Statistics, 92(4), 1004–1016.
#' Alessi, L., Barigozzi, M., & Capasso, M. (2010, December 1). Improved penalization for determining the number of factors in approximate factor models. Statistics Probability Letters, 80(23-24), 1806–1813
#' Ahn, S., & Horenstein, A. (2013, May 1). Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203–1227
#' Bai, J., & Ng, S. (2003). Determining the number of factors in approximate factor models.
#' @import fnets
#' @import dfms
#' @name determine_number_of_static_factors
#' @export
determine_number_of_static_factors <- function(data, method = "on", max_factors = 10) {
  if(sum(is.na(data)) > 0){
    stop("Data contains missing values")
  }
  if(!is.matrix(data)){
    data <- as.matrix(data)
  }
  if (!is.numeric(data)) {
    stop("Data must be numeric")
  }
    
    switch(method,
                     "on" = {
                         cat("Number of factors determined by the Onatski (2010) test\n")
                         return(Onatski_2010_test(data, max_factors))
                     },
                     "al" = {
                         cat("Number of factors determined by the modified Alessi, Barigozzi and Capasso (2010)\n")
                         return(fnets_ic_wrapper(fnets::factor.number(data, fm.restricted = TRUE, method = "ic", q.max = max_factors)))
                     },
                     "ah" = {
                         cat("Number of factors determined by the modified Ahn and Horenstein (2013)\n")
                         return(fnets_ic_wrapper(fnets::factor.number(data, fm.restricted = TRUE, method = "er", q.max = max_factors)))
                     },
                     "bn" = {
                         cat("Number of factors determined by the modified Bai and Ng (2002)\n")
                         return(dfms_ic_wrapper(dfms::ICr(data, max_factors)))
                     },
                     stop("Invalid method"))
}

#' Wrap FNets IC Results
#'
#' This function wraps the results from the FNets IC method into a data frame.
#'
#' @param fnets_ic A list or matrix containing the FNets IC results.
#' @return A data frame with the number of factors determined by the FNets IC method.
#' @references
#' Avarucci, M., Cavicchioli, M., Mario, F., & Zaffaroni, P. (2021, 01). The main business cycle shock(s). frequency-band estimation of the number of dynamic factors.
#' @name fnets_ic_wrapper
#' @export
fnets_ic_wrapper <- function(fnets_ic){
  df <- as.data.frame(fnets_ic)
  colnames(df) <- "factors"
  rownames(df) <- paste0("IC", seq_len(nrow(df)))
  return(df)
}

#' Wrap DFMS IC Results
#'
#' This function wraps the results from the DFMS IC method into a data frame.
#'
#' @param dfms_ic A list containing the DFMS IC results, specifically the "r.star" element.
#' @return A data frame with the number of factors determined by the DFMS IC method.
#' @references
#' Bai, J., & Ng, S. (2002). Determining the Number of Factors in Approximate Factor Models. \emph{Econometrica}, 70(1), 191-221.
#' @name dfms_ic_wrapper
#' @export
dfms_ic_wrapper <- function(dfms_ic){
  df <- as.data.frame(dfms_ic[["r.star"]])
  colnames(df) <- "r"
  return(df)
}