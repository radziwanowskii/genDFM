#' Estimate Factor-Augmented Vector Autoregression (FAVAR)
#'
#' This function estimates a Factor-Augmented Vector Autoregression (FAVAR) model.
#'
#' @param data A data frame or matrix of data to be analyzed.
#' @param r Number of static factors. Default is 3.
#' @param n_lags Number of lags in the VAR model. Default is 13.
#' @param slow_indices Indices of slow-moving variables. Default is NULL.
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
FAVAR <- function(data, r = 3, n_lags = 13, slow_indices = NULL, direct_indices){
  # Convert the data to a matrix
  data_matrix <- as.matrix(data)
  # Create matrix of direct variables
  data_direct <- as.matrix(data[,direct_indices])

  # Check if slow indices are provided
  if(is.null(slow_indices)){
  }
  else{
    data_slow <- as.matrix(data[,slow_indices])
    favar <- favar_prep_arma(data_matrix,r,n_lags,data_slow,data_direct)
  }
  # Name the columns of the data
  # Assign names to the data_for_var matrix
   
  factor_names <- paste0("factor ", 1:r)
  if(is.numeric(direct_indices)){
  direct_var_names <- colnames(data_matrix)[direct_indices]
  }
  else{
  direct_var_names <- direct_indices
  }
  colnames(favar[[1]]) <- c(factor_names, direct_var_names)
  colnames(favar[[2]]) <- c(factor_names, direct_var_names)
  colnames(favar[[3]]) <- factor_names
  colnames(favar[[4]]) <- colnames(data_matrix)
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
#' @param shock The type of shock ("unit" or "sd"). Default is "unit".
#' @param plot Whether to plot the IRF. Default is TRUE.
#' @param alpha Significance level for the confidence intervals. Default is 0.05.
#' @return A list containing the IRF and bootstrapped IRF:
#' \itemize{
#'   \item \code{irf}: The impulse response function.
#'   \item \code{boot}: The bootstrapped impulse response function.
#' }
#' @references
#' Bernanke, B. S., Boivin, J., & Eliasz, P. (2005, 02). Measuring the Effects of Monetary Policy: A Factor-Augmented Vector Autoregressive (FAVAR) Approach*. The Quarterly Journal of Economics, 120(1), 387-422
#' @name favar_IRF
#' @export
favar_IRF <- function(favar, shock_variable, response_variable = NULL,
n_ahead = 49, n_boot = 500, cumulative = FALSE, shock = "unit", plot = TRUE, alpha = 0.05){
  
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
  if(plot == TRUE){
    plot_IRF(irf_result,boot_result,alpha)
  }
 
  # Return the IRF and the bootstrapped IRF
irf.out <- list(irf = irf_result, boot = boot_result)
class(irf.out) <- "IRF"

return(irf.out)
}

#' Plot Impulse Response Functions (IRFs)
#'
#' This function plots the Impulse Response Functions (IRFs) along with their confidence intervals.
#'
#' @param irf The impulse response function.
#' @param boot_irf The bootstrapped impulse response function.
#' @param alpha Significance level for the confidence intervals. Default is 0.05.
#' @return A plot of the IRFs with confidence intervals.
#' @name plot_IRF
#' @export
plot_IRF <- function(irf, boot_irf, alpha = 0.05) {
  n_ahead <- length(irf)
  quantiles <- row_quantiles(boot_irf, c(alpha, 1-alpha))
    
  p <- ggplot(data = data.frame(x = 1:n_ahead, lower = quantiles[, 1], upper = quantiles[, 2], irf = irf)) +
    geom_ribbon(aes(x = x, ymin = lower, ymax = upper), fill = "#64b675", alpha = 0.5) +
    geom_line(aes(x = x, y = irf), color = "#0b6306", linewidth = 1.5) +
    geom_hline(yintercept = 0, linetype = "twodash", color = "black", size = 1) +
    theme(
      #panel.background = element_rect(fill = "#d3d3d3"),  # Darker background for the plot panel
     # plot.background = element_rect(fill = "#a9a9a9"),   # Darker background for the entire plot
      axis.title = element_text(size = 14, face = "bold", color = "black"),  # Customize axis titles
      axis.text = element_text(size = 12, color = "black")  # Customize axis text
    ) +
    labs(x = "Horizon", y = "Y Axis Label")  # Add axis labels

    print(p)  # Print the plot
}