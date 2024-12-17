#' Estimate Dynamic Factor Model (DFM)
#'
#' This function estimates a Dynamic Factor Model (DFM) using various methods.
#'
#' @param data a T x N dimensional matrix of data to be analyzed. Data can be in the form of a matrix or data frame.
#' @param r Number of static factors in the model. If NULL, it will be determined automatically. 
#' @param q Number of dynamic factors in the model. If NULL, it will be determined automatically. 
#' @param q_max Maximum number of dynamic factors to consider. Default is 10.
#' @param M Number of lags for the spectral density matrix. If NULL, it will be calculated as (2/3) * length(data[,1])^(1/3). Default is NULL.
#' @param model The model to use for estimation. Must be one of "fhlr_2000", "fhlr_2005", "sw_2002", "doz_2011", "doz_2012", or "BM_2014". Default is "fhlr_2000".
#' \itemize{
#'   \item \code{"fhlr_2000"}: {Generalized Dynamic-Factor Model (Forni, Hallin, Lippi, and Reichlin, 2000).}
#'   \item \code{"fhlr_2005"}: {Restricted Generalized Dynamic Factor Model (Forni, Hallin, Lippi, and Reichlin 2005).}
#'   \item \code{"sw_2002"}: {Principal Components method (Stock and Watson, 2002).}
#'   \item \code{"doz_2011"}: Doz, Giannone, and Reichlin (2011) - relies on the `dfms` package.
#'   \item \code{"doz_2012"}: Doz, Giannone, and Reichlin (2012) - relies on the `dfms` package.
#'   \item \code{"BM_2014"}: Banbura and Modugno (2014) - relies on the `dfms` package.
#' }
#' @return A list containing the results of the DFM estimation. The structure of the list depends on the chosen model.
#' @references
#' Forni, M., Hallin, M., Lippi, M., & Reichlin, L. (2000). The generalized dynamic-factormodel: Identification and estimation. The Review of Economics and Statistics.
#' Forni, M., Hallin, M., Lippi, M., & Reichlin, L. (2005). The generalized dynamicfactor model. Journal of the American Statistical Association, 100(471), 830–840.
#' Stock, J. H., & Watson, M. W. (2002a). Forecasting using principal componentsfrom a large number of predictors. Journal of the American Statistical Association, 97(460), 1167–1179.
#' Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for large ap-proximate dynamic factor models based on kalman filtering. Journal of Econometrics, 164(1), 188-205. 
#' Doz, C., Giannone, D., & Reichlin, L. (2012, 11). A Quasi–Maximum LikelihoodApproach for Large, Approximate Dynamic Factor Models. The Review of Economics and Statistics, 94(4), 1014-1024.
#'Banbura, M., & Modugno, M. (2013). Maximum likelihood estimation of large factor model on datasets with arbitrary pattern of missing data. Journal of Applied Econometrics.
#' @import dfms
#' @name estimate_DFM
#' @export
estimate_DFM <- function(data, r = NULL, q = NULL, q_max = 10, M = NULL, model = "fhlr_2000"){
  if(sum(is.na(data)) > 0){
    stop("Data contains missing values")
  }
  if(is.null(M)){
    M <- round((2/3)*length(data[,1])^(1/3))
  }
  if(!is.matrix(data)){
    data <- as.matrix(data)
  }
  if (!is.numeric(data)) {
    stop("Data must be numeric")
  }
  scaled_data <- scale(data)
  switch(model,
         "fhlr_2000" = {
           if(is.null(q)){
             q <- as.numeric(determine_number_of_dynamic_factors(scaled_data,"av",q_max))
           }
           gdfm <- gdfm_forni_2000(data,M,q)
           generalized_dfm.out <- list(filter_estimate = gdfm[[1]], common_component = gdfm[[2]])
           class(generalized_dfm.out) <- "generalized_dfm"
           return(generalized_dfm.out)
         },
         "fhlr_2005" = {
           if(is.null(r)){
             r <- determine_number_of_static_factors(scaled_data,"on",q_max)
           }
           if(is.null(q)){
             q <- as.numeric(determine_number_of_dynamic_factors(scaled_data,"av",q_max))
           }
           restridcted_gdfm <- restricted_gdfm(data,M,q,r)
           generalized_dfm.out <- list(covariance_matrix = restridcted_gdfm[[1]], 
                                       Z = restridcted_gdfm[[2]], Factors = restridcted_gdfm[[3]], Common_component = restridcted_gdfm[[4]],
                                       Max_lag = restridcted_gdfm[[5]], Original_data = restridcted_gdfm[[6]])
           class(generalized_dfm.out) <- "restricted_gdfm"
           return(generalized_dfm.out)
         },
         "sw_2002" = {
           if(is.null(r)){
             r <- determine_number_of_static_factors(scaled_data,"on",q_max)
           }
           factors <- pca_cpp_arma(data, n_components = r)
           pca.out <- list(factors = factors, data = data)
           class(pca.out) <- "sw_2002"
           return(pca.out)
         },
         "doz_2011" = {
           if(is.null(r)){
             r <- determine_number_of_static_factors(scaled_data,"on",q_max)
           }
           return(dfms::DFM(data, r, em.method = "none"))
         },
         "doz_2012" = {
           if(is.null(r)){
             r <- determine_number_of_static_factors(scaled_data,"on",q_max)
           }
           return(dfms::DFM(data, r, em.method = "none"))
         },
         "BM_2014" = {
           if(is.null(r)){
             r <- determine_number_of_static_factors(scaled_data,"on",q_max)
           }
           return(dfms::DFM(data, r, em.method = "BM"))
         },
         stop("Invalid model specified")
  )
}
#' Plot Restricted Generalized Dynamic Factor Model (restricted_gdfm)
#'
#' This function plots the components of a Restricted Generalized Dynamic Factor Model (restricted_gdfm).
#'
#' @param x An object of class "restricted_gdfm".
#' @param main Title of the plot. Default is NULL.
#' @param xlab Label for the x-axis. Default is "Time".
#' @param ylab Label for the y-axis. Default is "Factors".
#' @return A plot of the DFM components.
#' @name plot.restricted_gdfm
#' @export
plot.restricted_gdfm <- function(x, main = NULL, xlab = "Time", ylab = "Factors") {
  plot_series_matrix(x[[3]], main = main, xlab = xlab, ylab = ylab)
}
#' Plot Generalized Dynamic Factor Model (generalized_dfm)
#'
#' This function plots the components of a Generalized Dynamic Factor Model (generalized_dfm).
#'
#' @param x An object of class "generalized_dfm".
#' @param main Title of the plot. Default is NULL.
#' @param xlab Label for the x-axis. Default is "Time".
#' @param ylab Label for the y-axis. Default is "Common Component".
#' @return A plot of the DFM components.
#' @name plot.generalized_dfm
#' @export
plot.generalized_dfm <- function(x, main = NULL, xlab = "Time", ylab = "Common Component") {
  plot_series_matrix(x[[2]], main = main, xlab = xlab, ylab = ylab)
}




#' Predict Restricted Generalized Dynamic Factor Model (restricted_gdfm)
#'
#' This function makes predictions based on a Restricted Generalized Dynamic Factor Model (restricted_gdfm) object.
#'
#' @param object An object of class "restricted_gdfm".
#' @param h The forecast horizon.
#' @param ... Additional arguments (not used).
#' @return The forecasted values or a message indicating that prediction is not supported for the given class type.
#' @name predict.restricted_gdfm
#' @export
predict.restricted_gdfm <- function(object, h, ...) {
  M <- object[[5]]
  if (h > M) {
    print("Only forecasts up to Max lag M are allowed")
    return(NULL)
  } else {
    cov_cube <- object[[1]]
    Z <- object[[2]]
    original_data <- object[[6]]
    x_T <- original_data[nrow(original_data), ]
    return(forecast_projection(cov_cube[,,M+1], cov_cube[,,M+1+h], Z, x_T))
  }
}

plot_series_matrix <- function(mat, time = NULL, main_title = "Time Series Plot", ylab = "Value", xlab = "Time") {
  if (!is.matrix(mat)) {
    stop("Input must be a matrix.")
  }
  
  # Use default time if not provided
  if (is.null(time)) {
    time <- 1:nrow(mat)
  }
  
  # Check if time matches rows
  if (length(time) != nrow(mat)) {
    stop("Length of 'time' must match the number of rows in the matrix.")
  }
  
  # Generate Factor labels
  factor_labels <- paste("Factor", seq_len(ncol(mat)))
  
  # Create the plot
  matplot(time, mat, type = "l", lwd = 2, lty = 1, col = rainbow(ncol(mat)),
          xlab = xlab, ylab = ylab, main = main_title)
  
  # Add a legend with Factor labels
  legend("topright", legend = factor_labels,
         col = rainbow(ncol(mat)), lty = 1, lwd = 2)
}
