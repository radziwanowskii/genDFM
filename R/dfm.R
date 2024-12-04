#' Estimate Dynamic Factor Model (DFM)
#'
#' This function estimates a Dynamic Factor Model (DFM) using various methods.
#'
#' @param data A matrix of data to be analyzed.
#' @param r Number of static factors. If NULL, it will be determined automatically. Default is NULL.
#' @param q Number of dynamic factors. If NULL, it will be determined automatically. Default is NULL.
#' @param q_max Maximum number of dynamic factors to consider. Default is 10.
#' @param M Number of lags for the spectral density matrix. If NULL, it will be calculated as (2/3) * length(data[,1])^(1/3). Default is NULL.
#' @param model The model to use for estimation. Must be one of "fhlr_2000", "fhlr_2005", "sw_2002", "doz_2011", "doz_2012", or "BM_2014". Default is "fhlr_2000".
#' \itemize{
#'   \item \code{"fhlr_2000"}: Forni, Hallin, Lippi, and Reichlin (2000).
#'   \item \code{"fhlr_2005"}: Forni, Hallin, Lippi, and Reichlin (2005).
#'   \item \code{"sw_2002"}: Stock and Watson (2002).
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
  if(is.null(M)){
    M <- round((2/3)*length(data[,1])^(1/3))
  }
  if(!is.matrix(data)){
    data <- as.matrix(data)
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
#' Plot Dynamic Factor Model (DFM)
#'
#' This function plots the components of a Dynamic Factor Model (DFM) based on its class type.
#'
#' @param x An object of class "generalized_dfm", "restricted_gdfm", "dfm", or "sw_2002".
#' @param ... Additional arguments passed to the plot function.
#' @return A plot of the DFM components.
#' @name plot.DFM
#' @export
plot.DFM <- function(x, ...) {
  if(class(x) == "generalized_dfm"){
    plot_series_matrix(x[[2]], "Common Component", ...)
  }
  else if(class(x) == "restricted_gdfm"){
    plot_series_matrix(x[[3]], "Factors", ...)
  }
  else if(class(x) == "dfm"){
    plot(x, ...)
  }
  else if(class(x) == "sw_2002"){
    plot_series_matrix(x[[1]], "Factors", ...)
  }
  else{
    stop("Unknown class type")
  }
}

#' Predict Dynamic Factor Model (DFM)
#'
#' This function makes predictions based on a Dynamic Factor Model (DFM) object.
#'
#' @param dfm An object of class "generalized_dfm", "restricted_gdfm", "dfm", or "sw_2002".
#' @param h The forecast horizon.
#' @return The forecasted values or a message indicating that prediction is not supported for the given class type.
#' @name predict.DFM
#' @export
predict.DFM <- function(dfm, h) {
  class_type <- class(dfm)
  
  result <- switch(class_type,
                   "generalized_dfm" = {
                     print("No prediction for generalized_dfm")
                     NULL
                   },
                   "restricted_gdfm" = {
                     M <- dfm[[5]]
                     if(h > M){
                       print("Only forecasts up to Max lag M are allowed")
                       NULL
                     } else {
                       cov_cube <- dfm[[1]]
                       Z <- dfm[[2]]
                       original_data <- dfm[[6]]
                       x_T <- original_data[nrow(original_data),]
                       forecast_projection(cov_cube[,,M+1], cov_cube[,,M+1+h], Z, x_T)
                     }
                   },
                   "dfm" = {
                     predict(dfm, n.ahead = h)
                   },
                   "sw_2002" = {
                     print("Currently we support no prediction for sw_2002")
                     NULL
                   },
                   {
                     stop("Unknown class type")
                   }
  )
  return(result)
}