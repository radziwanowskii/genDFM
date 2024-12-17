#' @import Rcpp
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @useDynLib genDFM, .registration = TRUE
NULL
#' Estimate Factor-Augmented Vector Autoregression (FAVAR)
#'
#' This function estimates a Factor-Augmented Vector Autoregression (FAVAR) model.
#'
#' @param data A data frame or matrix of data to be analyzed.
#' @param r Number of static factors. Default is 3.
#' @param n_lags Number of lags in the VAR model. Default is 13.
#' @param slow_indices Indices of slow-moving variables. 
#' @param direct_indices Indices of directly observed variables.
#' @return A list containing the FAVAR model components:
#' \itemize{
#'   \item \code{data_for_var}: Data prepared for VAR.
#'   \item \code{var}: VAR model.
#'   \item \code{factors}: Estimated factors.
#'   \item \code{loadings}: Factor loadings.
#' }
#' @references
#' Bernanke, B. S., Boivin, J., & Eliasz, P. (2005, 02). Measuring the Effects of Monetary Policy: A Factor-Augmented Vector Autoregressive (FAVAR) Approach*. The Quarterly Journal of Economics, 120(1), 387-422
#' @name FAVAR
#' @export
FAVAR <- function(data, r = 3, n_lags = 13, slow_indices = NULL, direct_indices = NULL){
  # Check if the data contains missing values
  if(sum(is.na(data)) > 0){
    stop("Data contains missing values")
  }
  # Convert the data to a matrix
  if(!is.matrix(data)){
    data <- as.matrix(data)
  }
  if (!is.numeric(data)) {
    stop("Data must be numeric")
  }
  # Check if slow indices are provided
  if(is.null(slow_indices) || length(slow_indices) < 2){
    stop("Slow indices must be provided")
  }
   # Check if direct indices are provided
  if(is.null(direct_indices) || length(direct_indices) < 2){
    stop("Direct indices must be provided")
  }

  # Create matrices of direct and slow variables
  data_direct <- data[,direct_indices]
  data_slow <- data[,slow_indices]
  
  # Run the C++ code
  favar <- favar_prep_arma(data,r,n_lags,data_slow,data_direct)
  
  # Name the columns of the data
  # Assign names to the data_for_var matrix
   
  factor_names <- paste0("factor ", 1:r)
  if(is.numeric(direct_indices)){
  direct_var_names <- colnames(data)[direct_indices]
  }
  else{
  direct_var_names <- direct_indices
  }
  colnames(favar[[1]]) <- c(factor_names, direct_var_names)
  colnames(favar[[2]]) <- c(factor_names, direct_var_names)
  colnames(favar[[3]]) <- factor_names
  colnames(favar[[4]]) <- colnames(data)
  favar.out <- list(data_for_var = favar[[1]], var = favar[[2]], factors = favar[[3]], loadings = favar[[4]])
  class(favar.out) <- "FAVAR"
  return(favar.out)
}

#' Compute Impulse Response Functions (IRFs) for FAVAR
#'
#' This function computes the Impulse Response Functions (IRFs) for a Factor-Augmented Vector Autoregression (FAVAR) model.
#'
#' @param favar An object of class "FAVAR".
#' @param shock_variable The variable to shock.
#' @param response_variable The variable to respond. Default is NULL.
#' @param n_ahead Number of periods ahead to compute the IRF. Default is 49.
#' @param n_boot Number of bootstrap samples. Default is 500.
#' @param cumulative Whether to compute cumulative IRFs. Default is FALSE.
#' @param shock The type of shock used in the model. "unit" implements a one-unit shock to the schock variable, "sd" implements a shock of one standard deviation. Default is "unit".
#' @param plot Whether to plot the IRF. Default is TRUE.
#' @param alpha Significance level for the confidence intervals in the plot. Default is 0.05.
#' @return A list containing the IRF and bootstrapped IRF:
#' \itemize{
#'   \item \code{irf}: The impulse response function.
#'   \item \code{boot}: The bootstrapped values of the impulse response function.
#' }
#' @references
#' Bernanke, B. S., Boivin, J., & Eliasz, P. (2005, 02). Measuring the Effects of Monetary Policy: A Factor-Augmented Vector Autoregressive (FAVAR) Approach*. The Quarterly Journal of Economics, 120(1), 387-422
#' @name favar_IRF
#' @export
favar_IRF <- function(favar, shock_variable, response_variable = NULL,
n_ahead = 49, n_boot = 500, cumulative = FALSE, shock = "unit", plot = TRUE, alpha = 0.05){
  
  if(class(favar) != "FAVAR"){
    stop("Input must be of class FAVAR")
  }

  boot_data <- favar[[1]]
  VAR <- favar[[2]]
  loadings <- favar[[4]]
  n_lags <- (length(VAR[,1])-1)/length(VAR[1,])

   # Number of directly observed variables
  num_direct_vars <- length(VAR[1,]) - length(loadings[,1])
   # Create the loadings matrix with additional rows of zeros for each direct variable
  loadings_for_irf <- rbind(loadings, matrix(0, nrow = num_direct_vars, ncol = ncol(loadings)))

  #Check if the shock variable is a string or a number
  if(is.numeric(shock_variable)){
    shock_index <- shock_variable
  }
  else{
    shock_index <- match(shock_variable,colnames(VAR))
  }
   # Check if the shock is a unit or a sd, if sd then compute
  if(shock == "unit"){
    shock_size <- 1
  }
  else if(shock == "sd"){
    var_residuals <- compute_residuals(boot_data,VAR,n_lags)
    shock_size <- sd(var_residuals[,shock_index])
  }
  # Compute IRFs and bootstrapped IRFs depending if the variable is directly observed
  if(!is.na(match(response_variable,colnames(VAR)))){
  all_irf <- compute_irf_with_constant_vertical(VAR,shock_index-1,n_ahead,shock_size)
  boot_irf <- bootstrap_irf_arma_optimal(boot_data,VAR,loadings_for_irf, shock_index-1,n_ahead,n_boot,n_lags,shock_size,direct = TRUE)
  }
  else{
  all_irf <- compute_all_irf(VAR,loadings_for_irf,shock_index-1,n_ahead,shock_size)
  boot_irf <- bootstrap_irf_arma_optimal(boot_data,VAR,loadings_for_irf,shock_index-1,n_ahead,n_boot,n_lags,shock_size,direct = FALSE)
  }
  
  # Cumulate the IRF if necessary
  if(cumulative == TRUE){
    all_irf <- apply(all_irf, 2, cumsum)
    boot_irf <- apply(boot_irf, c(2, 3), cumsum)
  }

  # Return the IRF of the chosen variables
  if(is.null(response_variable)){
    irf_result <- all_irf
    boot_result <- boot_irf
  }
  else if(!is.na(match(response_variable,colnames(VAR)))){
  
    response_variable_index <- match(response_variable,colnames(VAR))
    irf_result <- all_irf[,response_variable_index]
    boot_result <- boot_irf[,response_variable_index,]
  }
  else{
    response_variable_index <- match(response_variable,colnames(loadings))
    irf_result <- all_irf[,response_variable_index]
    boot_result <- boot_irf[,response_variable_index,]
  }
  # Return the IRF and the bootstrapped IRF
irf.out <- list(irf = irf_result, boot = boot_result)
class(irf.out) <- "irf"

if(plot == TRUE){
    #plot.irf(irf_result,boot_result,alpha)
    plot.irf(irf.out,alpha)
  }

return(irf.out)
}

#' Plot Impulse Response Functions (IRFs)
#'
#' This function plots the Impulse Response Functions (IRFs) along with their confidence intervals.
#'
#' @param x An object of class "irf".
#' @param alpha Significance level for the confidence intervals in the plot. Default is 0.05.
#' @return A plot of the IRFs with confidence intervals.
#' @name plot.irf
#' @export
plot.irf <- function(x, alpha = 0.05) {
  irf <- x[[1]]
  boot_irf <- x[[2]]
  
  n_ahead <- length(irf)
  quantiles <- row_quantiles(boot_irf, c(alpha, 1 - alpha))
  
  # Create the plot
  plot(1:n_ahead, irf, type = "l", col = "#0b6306", lwd = 1.5, ylim = range(c(quantiles, irf)),
       xlab = "Horizon", ylab = "Response", main = "Impulse Response Functions")
  
  # Add the confidence intervals
  polygon(c(1:n_ahead, rev(1:n_ahead)), c(quantiles[, 1], rev(quantiles[, 2])),
          col = rgb(100, 182, 117, maxColorValue = 255, alpha = 127), border = NA)
  
  # Add the IRF line again to ensure it is on top
  lines(1:n_ahead, irf, col = "#0b6306", lwd = 1.5)
  
  # Add a horizontal line at y = 0
  abline(h = 0, lty = "twodash", col = "black", lwd = 1)
}