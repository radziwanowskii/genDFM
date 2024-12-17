#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h> 
#include <complex>
#include <algorithm>
#include <utility>
#include <random>

using namespace arma;


// [[Rcpp::export]]
arma::mat singular_boot_for_irf(const arma::mat& Y, const arma::mat& VAR, const arma::mat& residuals) {
  
  int k = Y.n_cols;            // Number of variables
  int p = (VAR.n_rows - 1) / k; // Number of lags (subtract 1 for the constant row)
  int n_obs = Y.n_rows;         // Number of observations
  
  // Initialize the bootstrapped time series data with the first p observations from the original data
  arma::mat Y_bootstrap(n_obs, k, arma::fill::zeros);
  Y_bootstrap.rows(0, p - 1) = Y.rows(0, p - 1);  // Copy the first p rows from original data
  
  // Resample residuals with replacement using RcppArmadillo's sample function
  arma::mat resampled_residuals(n_obs - p, k);
  
  // Generate a vector of random indices for resampling (using RcppArmadillo::sample)
  arma::uvec indices = Rcpp::RcppArmadillo::sample(arma::regspace<arma::uvec>(0, n_obs - p - 1), n_obs - p, true);
  
  for (int i = 0; i < n_obs - p; ++i) {
    resampled_residuals.row(i) = residuals.row(indices[i]);
  }
  
  // Generate the bootstrapped series for each time step starting from p
  for (int t = p; t < n_obs; ++t) {
    
    // Initialize the new observation with the constant term
    arma::vec Y_new = VAR.row(0).t();  // The first row contains the constants
    
    // Add the contribution from the lagged terms (VAR coefficients)
    for (int lag = 1; lag <= p; ++lag) {
      // Get the coefficients for the corresponding lag (stacked format, skip the constant row)
      arma::mat VAR_block = VAR.rows(1 + (lag - 1) * k, lag * k);  // Rows for the given lag
      
      // Add the lagged effect to the new observation
      Y_new += VAR_block.t() * Y_bootstrap.row(t - lag).t();  // Transpose the block to match dimensions
    }
    
    // Add the resampled residual
    Y_new += resampled_residuals.row(t - p).t();
    
    // Store the new bootstrapped observation
    Y_bootstrap.row(t) = Y_new.t();
  }
  
  return Y_bootstrap;
}

// [[Rcpp::export]]
arma::mat compute_irf_with_constant_vertical(const arma::mat& VAR, int shock_var, int horizon, double shock_size = 1) {
  int k = VAR.n_cols;  // Number of variables in the VAR system
  int p = (VAR.n_rows - 1) / k;  // Number of lags (first row is the constant)
  
  // Extract the constant vector (first row of VAR)
  arma::vec constant = VAR.row(0).t();
  
  // Initialize the matrix to store the IRFs for all variables at each time step
  arma::mat irf_matrix(k, horizon, arma::fill::zeros);
  
  // Initial impulse: only shock the specified variable at time 0
  irf_matrix(shock_var, 0) = shock_size;  // Set shock to 1 for the chosen variable at t=0
  
  // Loop over the horizon (starting at time 1 because the initial shock is already stored at time 0)
  for (int h = 1; h < horizon; ++h) {
    
    // Initialize the response at time h
    arma::vec response = constant;  // Start with the constant term
    
    // Loop over lags
    for (int lag = 1; lag <= std::min(h, p); ++lag) {
      // Get the VAR coefficients for lag 'lag' (skip the first row since it's the constant)
      arma::mat var_coeff_block = VAR.rows(1 + (lag - 1) * k, lag * k);  // Coefficients for lag 'lag'
      
      // Add the contribution of the response from lag 'lag' to the current response at time 'h'
      response += var_coeff_block.t() * irf_matrix.col(h - lag);
    }
    
    // Store the computed response at time h (including constant for level shift)
    irf_matrix.col(h) = response;
  }
  
  // Return the IRF matrix containing the impulse res
  return irf_matrix.t();
}

// [[Rcpp::export]]
arma::mat compute_all_irf(const arma::mat& VAR, const arma::mat& loadings, int shock_var, int horizon, double shock_size = 1) {
  arma::mat irf = compute_irf_with_constant_vertical(VAR, shock_var, horizon, shock_size);
  return irf * loadings;
}


// [[Rcpp::export]]
arma::mat pca_cpp_arma(const arma::mat& X, int n_components = 0, bool center = true, bool scale = false) {
  arma::mat X_centered = X;
  
  // Center the matrix
  if (center) {
    arma::rowvec mean = arma::mean(X, 0);
    X_centered.each_row() -= mean;
  }
  
  // Scale the matrix
  if (scale) {
    arma::rowvec std_dev = arma::stddev(X_centered, 0, 0);
    X_centered.each_row() /= std_dev;
  }
  
  // Compute the covariance matrix
  arma::mat cov_matrix = arma::cov(X_centered);
  
  // Eigen decomposition
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  arma::eig_sym(eigenvalues, eigenvectors, cov_matrix);
  
  // Reverse the order of eigenvalues and eigenvectors
  eigenvalues = arma::reverse(eigenvalues);
  eigenvectors = arma::fliplr(eigenvectors);
  
  // Determine the number of components to return
  int max_components = eigenvalues.n_elem;
  if (n_components <= 0 || n_components > max_components) {
    n_components = max_components;
  }
  
  // Select the top n_components
  arma::vec selected_eigenvalues = eigenvalues.head(n_components);
  arma::mat selected_eigenvectors = eigenvectors.cols(0, n_components - 1);
  arma::mat scores = X_centered * selected_eigenvectors;
  
  return scores;
}

arma::mat VAR_estimate_arma(const arma::mat& Y, int p) {
  int n_obs = Y.n_rows;
  int n_vars = Y.n_cols;
  
  // Create the lagged matrix (X)
  arma::mat X(n_obs - p, n_vars * p + 1, arma::fill::ones);  // Initialize with ones for intercept
  for (int i = 0; i < p; i++) {
    X.cols(i * n_vars + 1, (i + 1) * n_vars) = Y.rows(p - i - 1, n_obs - i - 2);
  }
  
  // Response matrix
  arma::mat Y_lagged = Y.rows(p, n_obs - 1);
  
  // Estimate coefficients using OLS: B = (X'X)^-1 * X'Y
  arma::mat B = arma::inv(X.t() * X) * X.t() * Y_lagged;
  
  return B;
}
// [[Rcpp::export]]
arma::mat compute_residuals(const arma::mat& Y, const arma::mat& B, int p) {
  int n_obs = Y.n_rows;
  int n_vars = Y.n_cols;
  
  // Create the lagged matrix (X)
  arma::mat X(n_obs - p, n_vars * p + 1, arma::fill::ones);  // Initialize with ones for intercept
  for (int i = 0; i < p; i++) {
    X.cols(i * n_vars + 1, (i + 1) * n_vars) = Y.rows(p - i - 1, n_obs - i - 2);
  }
  
  // Response matrix
  arma::mat Y_lagged = Y.rows(p, n_obs - 1);
  
  // Calculate residuals
  arma::mat residuals = Y_lagged - X * B;
  
  return residuals;
}
arma::mat multivariate_multiple_regression_arma_2(const arma::mat& X, const arma::mat& Y) {
  // Estimate coefficients using OLS: B = (X'X)^-1 * X'Y
  arma::mat B = arma::inv(X.t() * X) * X.t() * Y;
  return B;
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat multivariate_multiple_regression_arma(const arma::mat& Y, 
                                                const arma::mat& X1, 
                                                const arma::mat& X2) {
  // Combine X1 and X2 into a single matrix X
  arma::mat X = arma::join_horiz(X1, X2);
  
  // Add a column of ones to X to include an intercept
  arma::mat X_with_intercept = arma::join_horiz(arma::ones<arma::vec>(X.n_rows), X);
  
  // Calculate the coefficients matrix B using the normal equation: B = (X'X)^-1 * X'Y
  arma::mat XtX = X_with_intercept.t() * X_with_intercept;
  arma::mat XtY = X_with_intercept.t() * Y;
  arma::mat B = arma::solve(XtX, XtY); // Using solve for stability
  
  return B;
}
// [[Rcpp::export]]
Rcpp::List favar_prep_arma(const arma::mat& data, int n_components,
                           int n_lags, const arma::mat& data_slow, const arma::mat& data_direct){
  
  // Subset the data matrix by the direct variable names
  arma::mat F_direct = data_direct;
  // Compute the principal components of all variables
  arma::mat C = pca_cpp_arma(data, n_components, false, false);
  
  // Compute the principal components of the slow variables
  arma::mat F_slow = pca_cpp_arma(data_slow, n_components, false, false);
  
  // Combine F_slow and F_direct into a single matrix
  //arma::mat F_slow_direct = join_horiz(F_slow, F_direct);
  
  arma::mat coefficients = multivariate_multiple_regression_arma(C, F_slow, F_direct);
  // Extract the coefficients of the F_direct variables
  int num_slow_vars = F_slow.n_cols;
  int num_direct_vars = F_direct.n_cols;
  
  arma::mat direct_coefficients = coefficients.submat(num_slow_vars + 1, 0, num_slow_vars + num_direct_vars, coefficients.n_cols - 1);
  
  // Compute the factors F
  arma::mat F = C - F_direct * direct_coefficients;
  
  // Combine F and F_direct into a single matrix
  arma::mat data_combined = join_horiz(F, F_direct);
  
  // Create the VAR
  arma::mat FAVAR = VAR_estimate_arma(data_combined, n_lags);
  
  // Compute Loadings 
  arma::mat loadings = multivariate_multiple_regression_arma_2(F, data);
  
  return Rcpp::List::create(
    Rcpp::Named("data_for_var") = data_combined,
    Rcpp::Named("FAVAR") = FAVAR,
    Rcpp::Named("Factors") = F,
    Rcpp::Named("Loadings") = loadings
  );
}
// [[Rcpp::export]]
arma::cube bootstrap_irf_arma_for_fyff(const arma::mat& boot_data,  const arma::mat& VAR,
                                       const arma::mat& loadings, int shock_var, int horizon, int no_boot, int n_lags, double shock_size = 1){
  // Initialize the cube to store the impulse response functions (IRFs) for each bootstrap sample
  arma::cube irf_cube(horizon, boot_data.n_cols, no_boot, arma::fill::zeros);
  // Compute the residuals for the original data
  arma::mat residuals = compute_residuals(boot_data, VAR, n_lags);
  for(int i = 0; i<no_boot; i++){
    
    // Resample the data using bootstrap
    arma::mat boot_sample = singular_boot_for_irf(boot_data, VAR, residuals);
    
    // Estimate the VAR coefficients using the bootstrapped sample
    arma::mat boot_VAR = VAR_estimate_arma(boot_sample, n_lags);
    
    // Compute the impulse response function (IRF) for the bootstrapped sample
    arma::mat irf_matrix = compute_irf_with_constant_vertical(boot_VAR, shock_var, horizon, shock_size);
    
    // Store the computed IRF in the cube
    irf_cube.slice(i) = irf_matrix;
  }
  return irf_cube;
}

// [[Rcpp::export]]
arma::cube bootstrap_irf_arma_optimal(const arma::mat& boot_data,  const arma::mat& VAR,
                                      const arma::mat& loadings, int shock_var, int horizon, int no_boot, int n_lags, double shock_size = 1,
                                      bool direct = false){
  // Initialize the cube to store the impulse response functions (IRFs) of the directly observed variables
  arma::cube irf_cube_direct(horizon, boot_data.n_cols, no_boot, arma::fill::zeros);
  // Compute the residuals for the original data
  arma::mat residuals = compute_residuals(boot_data, VAR, n_lags);
  for(int i = 0; i<no_boot; i++){
    
    // Resample the data using bootstrap
    arma::mat boot_sample = singular_boot_for_irf(boot_data, VAR, residuals);
    
    // Estimate the VAR coefficients using the bootstrapped sample
    arma::mat boot_VAR = VAR_estimate_arma(boot_sample, n_lags);
    
    // Compute the impulse response function (IRF) for the bootstrapped sample, arma::stddev(residuals.col(3))
    arma::mat irf_matrix = compute_irf_with_constant_vertical(boot_VAR, shock_var, horizon, shock_size);
    
    // Store the computed IRF in the cube
    //irf_cube.slice(i) = irf_matrix*loadings;
    // Store the directly observed variables in the cube
    irf_cube_direct.slice(i) = irf_matrix;
  }
  // Multiply each slice by loadings if needed
  if(direct == false) {
    // Initialize the cube to store the impulse response functions (IRFs) for each bootstrap sample
    arma::cube irf_cube_multiplied(horizon, loadings.n_cols, no_boot, arma::fill::zeros);
    for(int i = 0; i < irf_cube_direct.n_slices; i++) {
      irf_cube_multiplied.slice(i) = irf_cube_direct.slice(i) * loadings;
    }
    return irf_cube_multiplied;
  } else{
    return irf_cube_direct;
  }
}

// Function to calculate quantiles for each row of a matrix
// probs is a vector of quantiles (e.g., c(0.05, 0.95))
// [[Rcpp::export]]
arma::mat row_quantiles(const arma::mat& data, const arma::vec& probs) {
  int n_rows = data.n_rows;
  int n_cols = data.n_cols;
  int n_probs = probs.n_elem;
  
  // Initialize the matrix to store quantiles
  arma::mat quantiles(n_rows, n_probs, arma::fill::zeros);
  
  for (int i = 0; i < n_rows; ++i) {
    arma::rowvec row = data.row(i);
    arma::rowvec sorted_row = arma::sort(row);
    
    for (int j = 0; j < n_probs; ++j) {
      double p = probs(j);
      double pos = p * (n_cols - 1);
      int lower = std::floor(pos);
      int upper = std::ceil(pos);
      double weight = pos - lower;
      
      if (lower == upper) {
        quantiles(i, j) = sorted_row(lower);
      } else {
        quantiles(i, j) = (1 - weight) * sorted_row(lower) + weight * sorted_row(upper);
      }
    }
  }
  
  return quantiles;
}
arma::cube compute_covariance_with_lags_and_store_in_cube(const arma::mat& data, int max_lag) {
  // Number of observations and variables
  int n_obs = data.n_rows;
  int n_vars = data.n_cols;
  
  // Initialize a cube to store covariance matrices for each lag
  arma::cube cov_cube(n_vars, n_vars, max_lag + 1, arma::fill::zeros);
  
  for (int lag = 0; lag <= max_lag; ++lag) {
    if (lag == 0) {
      // Covariance matrix for lag 0 (original data)
      cov_cube.slice(lag) = arma::cov(data);
    } else {
      // Create lagged version of the data
      arma::mat lagged_data = data.rows(lag, n_obs - 1);
      arma::mat original_data = data.rows(0, n_obs - lag - 1);
      
      // Compute the covariance matrix between original data and lagged data
      arma::mat cov_matrix = arma::cov(original_data, lagged_data);
      cov_cube.slice(lag) = cov_matrix;
    }
  }
  
  return cov_cube;
}
arma::cx_cube compute_dft_for_frequencies(const arma::cube& cov_cube, const arma::vec& frequencies,
                                          arma::vec weights, arma::vec k) {
  int n_slices = cov_cube.n_slices;
  int n_rows = cov_cube.n_rows;
  int n_cols = cov_cube.n_cols;
  int n_freqs = frequencies.n_elem;
  
  // Initialize a complex cube to store the DFT results for each frequency
  arma::cx_cube dft_cube(n_rows, n_cols, n_freqs, arma::fill::zeros);
  
  for (int f = 0; f < n_freqs; ++f) {
    double freq = frequencies(f);
    
    for (int lag = 0; lag < n_slices; ++lag) {
      // Compute the DFT manually for the specified frequency
      std::complex<double> exp_term = std::exp(std::complex<double>(0, -freq * k(lag)));
      dft_cube.slice(f) += cov_cube.slice(lag) * weights(lag) * exp_term;
    }
  }
  
  return dft_cube;
}
arma::cx_mat compute_first_q_eigenvectors(const arma::cx_mat& matrix, int q) {
  // Ensure the matrix is square
  if (matrix.n_rows != matrix.n_cols) {
    Rcpp::stop("The input matrix must be square.");
  }
  
  // Ensure q is less than or equal to the number of rows/columns
  if (q > matrix.n_rows) {
    Rcpp::stop("q must be less than or equal to the number of rows/columns of the matrix.");
  }
  
  // Initialize matrices to store eigenvalues and eigenvectors
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  
  // Compute eigenvalues and eigenvectors for a general (possibly complex) matrix
  bool status = arma::eig_gen(eigval, eigvec, matrix);
  
  // Check if the computation was successful
  if (!status) {
    Rcpp::stop("Eigen decomposition failed.");
  }
  
  // Return the first q eigenvectors
  return eigvec.cols(0, q - 1);
}
arma::cx_cube apply_eigenvectors_to_cube(const arma::cx_cube& cube, int q) {
  // Get the dimensions of the input cube
  int n_rows = cube.n_rows;
  int n_slices = cube.n_slices;
  
  // Initialize a cube to store the eigenvectors
  arma::cx_cube eigenvectors_cube(n_rows, q, n_slices, arma::fill::zeros);
  
  // Loop over each slice and apply the compute_first_q_eigenvectors function
  for (int i = 0; i < n_slices; ++i) {
    eigenvectors_cube.slice(i) = compute_first_q_eigenvectors(cube.slice(i), q);
  }
  
  return eigenvectors_cube;
}
arma::cx_mat compute_transposed_conjugate(const arma::cx_mat& matrix) {
  // Compute the complex conjugate
  arma::cx_mat conj_matrix = arma::conj(matrix);
  
  // Compute the transpose of the complex conjugate
  arma::cx_mat trans_conj_matrix = arma::trans(conj_matrix);
  
  return trans_conj_matrix;
}
arma::cx_cube apply_transposed_conjugate_to_cube(const arma::cx_cube& cube) {
  // Get the dimensions of the input cube
  int n_rows = cube.n_rows;
  int n_cols = cube.n_cols;
  int n_slices = cube.n_slices;
  
  // Initialize a cube to store the transposed conjugates
  arma::cx_cube trans_conj_cube(n_cols, n_rows, n_slices, arma::fill::zeros);
  
  // Loop over each slice and apply the compute_transposed_conjugate function
  for (int i = 0; i < n_slices; ++i) {
    trans_conj_cube.slice(i) = compute_transposed_conjugate(cube.slice(i));
  }
  
  return trans_conj_cube;
}

//alternative filter estimation giving nxn output
arma::cx_cube filter_estimate_alternative(const arma::cx_cube& transposed_conjugate_cube, const arma::cx_cube& eigenvectors_cube) {
  // Get the dimensions of the input cube
  int n_rows = eigenvectors_cube.n_rows;
  int n_slices = eigenvectors_cube.n_slices;
  
  // Initialize a cube to store the filtered estimates
  arma::cx_cube filtered_cube(n_rows, n_rows, n_slices, arma::fill::zeros);
  
  // Loop over each slice and apply the filtering operation
  for (int i = 0; i < n_slices; ++i) {
    filtered_cube.slice(i) =  eigenvectors_cube.slice(i) * transposed_conjugate_cube.slice(i);
  }
  
  return filtered_cube;
}
arma::cx_cube compute_inverse_dft_for_frequencies(const arma::cx_cube& dft_cube, const arma::vec& frequencies,
                                                  const arma::vec& k) {
  int n_rows = dft_cube.n_rows;
  int n_cols = dft_cube.n_cols;
  int n_freqs = frequencies.n_elem;
  
  // Initialize a complex cube to store the inverse DFT results for each frequency
  arma::cx_cube ifft_cube(n_rows, n_cols, k.n_elem, arma::fill::zeros);
  
  for (int i = 0; i < k.n_elem; i++){
    //int k_aux <- k[i]
    
    for (int f = 0; f < n_freqs; ++f) {
      double freq = frequencies(f);
      
      // Compute the inverse DFT manually for the specified frequency
      std::complex<double> exp_term = std::exp(std::complex<double>(0, freq * k(i)));
      ifft_cube.slice(i) += dft_cube.slice(f) * exp_term;
      
    }
  }
  // Normalize the result by the total number of elements
  ifft_cube /= static_cast<double>(n_freqs);
  
  return ifft_cube;
}

// Application of Forni et al. (2000) method to estimate dfms by dynamic principal components
arma::cube dynamic_pca(const arma::mat& matrix, int M, int q){
  
  //###1 Compute covariance matrices from -M to M
  
  // Compute the covariance matrix
  arma::cube cov_matrices = compute_covariance_with_lags_and_store_in_cube(matrix, M);
  arma::cube cov_matrices_transposed(cov_matrices.n_cols, cov_matrices.n_rows, cov_matrices.n_slices-1, arma::fill::zeros);
  
  // Transpose each matrix in the cube apart from the 0th lag
  for (int i = 0; i < cov_matrices.n_slices-1; ++i) {
    cov_matrices_transposed.slice(i) = arma::trans(cov_matrices.slice(cov_matrices.n_slices-i-1));
  }
  // create the two_sided sequence of covariance matrices
  arma::cube cov_sequence = join_slices(cov_matrices_transposed, cov_matrices);
  
  // compute the barlett weights
  arma::vec k = arma::abs(arma::linspace<arma::vec>(-M, M, 2*M+1));
  arma::vec weights = 1 - k/(M+1);
  
  // compute frequencies
  arma::vec h = arma::abs(arma::linspace<arma::vec>(0, 2*M, 2*M+1));
  arma::vec frequencies = (2*M_PI * h)/(2*M+1);
  
  //###2 Compute the sigma matrix by the Fourier transform
  // compute fourier transforms
  arma::cx_cube dft_cube = compute_dft_for_frequencies(cov_sequence, weights, frequencies, k);
  //###3 Compute eigenvectors and their transposed conjugates and estimate the filter K
  // compute the eigenvectors
  arma::cx_cube eigenvectors_cube = apply_eigenvectors_to_cube(dft_cube, q);
  
  // compute the transposed conjugates
  arma::cx_cube transposed_conjugate_cube = apply_transposed_conjugate_to_cube(eigenvectors_cube);
  
  // compute the filtered estimates
  //arma::cx_cube filtered_cube = filter_estimate(transposed_conjugate_cube, eigenvectors_cube);
  arma::cx_cube filtered_cube = filter_estimate_alternative(transposed_conjugate_cube, eigenvectors_cube);
  //###4 Compute the inverse Fourier transform to estimate the Filter K_
  // compute the inverse fourier transform
  arma::cx_cube ifft_cube = compute_inverse_dft_for_frequencies(filtered_cube, frequencies, k);
  
  // keep only the real part
  arma::cube ifft_cube_real = arma::real(ifft_cube);
  
  return ifft_cube_real;      
}
// [[Rcpp::export]]
Rcpp::List gdfm_forni_2000(const arma::mat& matrix, int M, int q){
  int t = matrix.n_rows - M;
  arma::mat subset_matrix = matrix.rows(t-M-1,t+M-1);
  arma::cube dfm = dynamic_pca(matrix, M, q);
  arma::vec filter(matrix.n_cols, arma::fill::zeros);
  for(int i = 0;i<(2*M+1);i++){
    filter += arma::trans(subset_matrix.row(i) * dfm.slice(i));
  }
  arma::mat common_component = matrix * filter;
  return Rcpp::List::create(
    Rcpp::Named("filter") = filter,
    Rcpp::Named("common_component") = common_component);
}
std::pair<arma::cx_mat, arma::cx_vec> estimate_row_eigenvectors(const arma::cx_mat& matrix) {
  // Ensure the matrix is square
  if (matrix.n_rows != matrix.n_cols) {
    Rcpp::stop("The input matrix must be square.");
  }
  
  // Initialize matrices to store eigenvalues and eigenvectors
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  
  // Compute eigenvalues and eigenvectors for a general (possibly complex) matrix
  bool status = arma::eig_gen(eigval, eigvec, matrix);
  
  // Check if the computation was successful
  if (!status) {
    Rcpp::stop("Eigen decomposition failed.");
  }
  
  // Transpose the eigenvector matrix to get row eigenvectors
  arma::cx_mat row_eigenvectors = eigvec.t();
  
  // Return a pair containing row eigenvectors and eigenvalues
  return std::make_pair(row_eigenvectors, eigval);
}
std::pair<arma::cx_cube, arma::cx_mat> apply_row_eigenvectors_to_cube(const arma::cx_cube& cube) {
  // Get the dimensions of the input cube
  int n_rows = cube.n_rows;
  int n_slices = cube.n_slices;
  
  // Initialize a cube to store the eigenvectors
  arma::cx_cube row_eigenvectors_cube(n_rows, n_rows, n_slices, arma::fill::zeros);
  arma::cx_mat eigenvalues_matrix(n_rows, n_slices, arma::fill::zeros);
  // Loop over each slice and apply the compute_first_q_eigenvectors function
  for (int i = 0; i < n_slices; ++i) {
    // Estimate row eigenvectors and eigenvalues for the current slice
    std::pair<arma::cx_mat, arma::cx_vec> result = estimate_row_eigenvectors(cube.slice(i));
    
    // Extract row eigenvectors and eigenvalues from the pair
    arma::cx_mat row_eigenvectors = result.first;
    arma::cx_vec eigenvalues = result.second;
    
    // Assign the extracted values to the corresponding slices
    row_eigenvectors_cube.slice(i) = row_eigenvectors;
    eigenvalues_matrix.col(i) = eigenvalues; // Assuming eigenvalues are stored in a column vector
  }
  
  return std::make_pair(row_eigenvectors_cube, eigenvalues_matrix);
}
// Estimate spectral density maastricht as in Forni et al. (2000)
arma::cx_cube spec_density_estimate(arma::cx_mat eigenvalues_matrix, arma::cx_cube row_eigenvectors_cube, arma::cx_cube transposed_conjugate_row_eigen){
  int n = row_eigenvectors_cube.slice(0).n_cols;
  int p = row_eigenvectors_cube.n_slices;
  int q = eigenvalues_matrix.n_rows;
  
  arma::cx_cube density_estimate(n, n, p, arma::fill::zeros);
  for(int i = 0; i < p; i++){
    for(int j=0;j<q;j++){
      density_estimate.slice(i) += eigenvalues_matrix(j,i) * transposed_conjugate_row_eigen.slice(i).col(j) * row_eigenvectors_cube.slice(i).row(j);
    }
  }
  return density_estimate;
}

// Function to compute inverse Fourier transform of the spectral density estimates
arma::cube compute_inverse_fourier(const arma::cx_cube& spectral_density) {
  int nrow = spectral_density.n_rows;
  int ncol = spectral_density.n_cols;
  int n_frequencies = spectral_density.n_slices;
  
  // Output cube to store the covariance matrices for each lag
  arma::cube covariance_matrices(nrow, ncol, n_frequencies, arma::fill::zeros);
  
  // Perform the inverse Fourier transform over the frequency dimension (3rd dimension)
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      // Extract the (i, j) element across all frequencies (this gives a vector of size n_frequencies)
      arma::cx_vec freq_vec = spectral_density.tube(i, j);
      
      // Compute the inverse FFT of the frequency vector
      arma::cx_vec ifft_result = arma::ifft(freq_vec);
      
      // Store the real part of the inverse FFT result (covariance values) into the covariance matrix cube
      covariance_matrices.tube(i, j) = arma::real(ifft_result);
    }
  }
  return covariance_matrices;
}

// [[Rcpp::export]]
arma::mat estimate_projection(const arma::mat& T_x0, const arma::mat& Z, const arma::mat& x_t) {
  
  // Step 1: Compute Z^T * T_x0 * Z
  arma::mat Z_T_Tx0_Z = Z.t() * T_x0 * Z;
  
  // Step 2: Invert Z^T * T_x0 * Z
  arma::mat Z_T_Tx0_Z_inv = arma::inv(Z_T_Tx0_Z);
  
  arma::mat Zx = x_t * Z;
  
  // Step 3: Compute the projection: [T_x0 Z (Z^T T_x0 Z)^(-1)] Z^T x_t
  arma::mat projection = Zx * (T_x0 * Z * Z_T_Tx0_Z_inv).t();
  
  return projection;
}
// [[Rcpp::export]]
arma::mat estimate_projection_singular(const arma::mat& T_x0,const arma::mat& T_xh, const arma::mat& Z, const arma::vec& x_T) {
  
  // Step 1: Compute Z^T * T_x0 * Z
  arma::mat Z_T_Tx0_Z = Z.t() * T_x0 * Z;
  
  //Create Z_x_T
  arma::vec Z_x_T = Z.t() * x_T;
  
  // Step 2: Invert Z^T * T_x0 * Z
  arma::mat Z_T_Tx0_Z_inv = arma::inv(Z_T_Tx0_Z);
  
  // Step 3: Compute the projection: [T_x0 Z (Z^T T_x0 Z)^(-1)] Z^T x_t
  arma::mat projection = (T_xh * Z * Z_T_Tx0_Z_inv) * Z_x_T;
  
  return projection;
}
// [[Rcpp::export]]
arma::mat forecast_projection(const arma::mat& T_x0,const arma::mat& T_xh, const arma::mat& Z, const arma::vec& x_T) {
  
  // Step 1: Compute Z^T * T_x0 * Z
  arma::mat Z_T_Tx0_Z = Z.t() * T_x0 * Z;
  
  //Create Z_x_T
  arma::vec Z_x_T = Z.t() * x_T;
  
  // Step 2: Invert Z^T * T_x0 * Z
  arma::mat Z_T_Tx0_Z_inv = arma::inv(Z_T_Tx0_Z);
  
  // Step 3: Compute the projection: [T_x0 Z (Z^T T_x0 Z)^(-1)] Z^T x_t
  arma::mat projection = (T_xh * Z * Z_T_Tx0_Z_inv) * Z_x_T;
  
  return projection;
}
// Application of Forni et al. (2005) method to estimate restricted gdfm
// M - max lag, q - number of common factors
// [[Rcpp::export]]
Rcpp::List restricted_gdfm(const arma::mat& matrix, int M, int q, int r){
  int n = matrix.n_cols;
  // Compute the covariance matrix
  arma::cube cov_matrices = compute_covariance_with_lags_and_store_in_cube(matrix, M);
  arma::cube cov_matrices_transposed(cov_matrices.n_cols, cov_matrices.n_rows, cov_matrices.n_slices-1, arma::fill::zeros);
  
  // Transpose each matrix in the cube apart from the 0th lag
  for (int i = 0; i < cov_matrices.n_slices-1; ++i) {
    cov_matrices_transposed.slice(i) = arma::trans(cov_matrices.slice(cov_matrices.n_slices-i-1));
  }
  
  // create the two_sided sequence of covariance matrices
  arma::cube cov_sequence = join_slices(cov_matrices_transposed, cov_matrices);
  
  // compute the barlett weights
  arma::vec k = arma::abs(arma::linspace<arma::vec>(-M, M, 2*M+1));
  arma::vec weights = 1 - k/(M+1);
  
  // compute frequencies
  arma::vec h = arma::abs(arma::linspace<arma::vec>(0, 2*M, 2*M+1));
  arma::vec frequencies = (2*M_PI * h)/(2*M+1);
  
  // compute fourier transforms
  arma::cx_cube dft_cube = compute_dft_for_frequencies(cov_sequence, weights, frequencies, k);
  
  // compute the row eigenvectors and eigenvalues
  std::pair<arma::cx_cube, arma::cx_mat> row_eigenvectors_list = apply_row_eigenvectors_to_cube(dft_cube);
  arma::cx_cube row_eigenvectors_cube = row_eigenvectors_list.first;
  arma::cx_mat eigenvalues_matrix = row_eigenvectors_list.second;
  
  // compute the transposed conjugates of the row eigenvectors
  arma::cx_cube transposed_conjugate_row_eigen = apply_transposed_conjugate_to_cube(row_eigenvectors_cube);
  
  // Subset to get only the first q eigenvalues, row eigenvectors, and conjugate transposes
  arma::cx_mat first_q_eigenvalues_matrix = eigenvalues_matrix(arma::span(0, q - 1), arma::span::all);    
  arma::cx_cube first_q_row_eigenvectors_cube = row_eigenvectors_cube(arma::span(0, q - 1), arma::span::all, arma::span::all);
  arma::cx_cube first_q_transposed_conjugate_row_eigen = transposed_conjugate_row_eigen(arma::span::all,arma::span(0, q - 1), arma::span::all);
  
  // Subset to get only the last n-q eigenvalues, row eigenvectors, and conjugate transposes
  arma::cx_mat rest_eigenvalues_matrix = eigenvalues_matrix(arma::span(q, n - 1), arma::span::all);
  arma::cx_cube rest_row_eigenvectors_cube = row_eigenvectors_cube(arma::span(q, n - 1), arma::span::all, arma::span::all);
  arma::cx_cube rest_transposed_conjugate_row_eigen = transposed_conjugate_row_eigen(arma::span::all, arma::span(q, n - 1), arma::span::all);
  
  // compute estimates of spectral density matrices
  // both cubes are of dimension nvar-nvar-2M+1
  arma::cx_cube common_density_estimate = spec_density_estimate(first_q_eigenvalues_matrix, first_q_row_eigenvectors_cube, first_q_transposed_conjugate_row_eigen);
  arma::cx_cube idiosyncratic_density_estimate = spec_density_estimate(rest_eigenvalues_matrix, rest_row_eigenvectors_cube, rest_transposed_conjugate_row_eigen);
  
  // testing: another function to compute the ift
  //covariance matix should be real (and they seem to be real)
  arma::cube common_covariance_2 = compute_inverse_fourier(common_density_estimate); 
  arma::cube idiosyncratic_covariance_2 = compute_inverse_fourier(idiosyncratic_density_estimate);
  
  // I am following the genealogy dfm paper which suggest: dynamic pca to recover the covariance matrix + static pca on it
  arma::mat static_pca = pca_cpp_arma(common_covariance_2.slice(M), r, false, false);
  
  // Estimate Zx_t vector as in Forni et al. (2005)
  arma::mat Zx_t = matrix * static_pca;
  
  // Estimate the common component
  arma::mat common_component = estimate_projection(common_covariance_2.slice(M), static_pca, matrix);
  
  
  return Rcpp::List::create(
    Rcpp::Named("common_covariance_matrx") = common_covariance_2,
    Rcpp::Named("Z") = static_pca,
    Rcpp::Named("Factors") = Zx_t,
    Rcpp::Named("common_component") = common_component,
    Rcpp::Named("Max_lag") = M,
    Rcpp::Named("Original_data") = matrix);
  
}

// [[Rcpp::export]]
arma::vec compute_sorted_eigenvalues(const arma::mat& matrix,arma::vec& sequence){
  
  int n_rows = matrix.n_rows;
  int n_cols = matrix.n_cols;
  int n_freqs = sequence.n_elem;
  // Compute the frequencies (as in the paper)
  arma::vec frequencies = 2 * M_PI * sequence / n_rows;
  arma::cx_mat dft_matrix(n_cols, n_freqs, arma::fill::zeros);
  for (int f = 0; f < n_freqs; ++f) {
    double freq = frequencies(f);
    // Compute the DFT manually for the specified frequency
   // std::complex<double> exp_term = std::exp(std::complex<double>(0, -2 * M_PI * freq / n_rows));
    for(int i=0; i<n_rows;i++){
      std::complex<double> exp_term = std::exp(std::complex<double>(0, -2 * M_PI * freq * (i+1) / n_rows));
      dft_matrix.col(f) += matrix.row(i).t() * exp_term;
    }
  }
  // Compute the transposed conjugate of the DFT matrix
  arma::cx_mat transposed_conjugate = compute_transposed_conjugate(dft_matrix);
  
  // Compute the periodogram
  arma::cx_mat periodogram(n_cols, n_cols, arma::fill::zeros);
  for(int i=0;i<n_freqs;i++){
    periodogram +=  dft_matrix.col(i) * transposed_conjugate.row(i) / (2 * M_PI * n_freqs);
  }
  //compute the eigenvalues
  arma::cx_vec eigval;
  arma::eig_gen(eigval, periodogram);
  arma::vec eigval_real = arma::real(eigval);
  arma::vec sorted_eigenvalues = arma::sort(arma::abs(eigval_real),"descend");
  
  return sorted_eigenvalues;
}

// [[Rcpp::export]]
double Onatski_2009_test(const arma::vec& sorted_eigenvalues, const arma::vec& CVs,
                         int k0, int k1){
  
  for(int k = k0-1; k < k1-2; k++){
    
    //change 10 exact value later
    arma::vec eigenvalue_ratios = arma::vec(k1-k-1, arma::fill::zeros);
    for(int i = k;i < k1-1; i++){
      
      eigenvalue_ratios(i-k) = (sorted_eigenvalues(i) - sorted_eigenvalues(i+1))/(sorted_eigenvalues(i+1) - sorted_eigenvalues(i+2));
    }
    
    double test_statistic = arma::max(eigenvalue_ratios);
    if(test_statistic <= CVs(k1-k-2)){
      cout<<"Hypothesis "<<(k+1)<<" is not rejected"<<endl;
      return (k+1);
    }
  }
  cout<<"All hypotheses rejected"<<endl;
  return k1;
}

// [[Rcpp::export]]
int Onatski_2010_test(const arma::mat& data,int rmax) {
  
  arma::mat covariance_matrix = arma::cov(data);
  
  // Compute the eigenvalues of the covariance matrix
  arma::vec eigenvalues;
  arma::eig_sym(eigenvalues, covariance_matrix);
  arma::vec sorted_eigenvalues = arma::sort(eigenvalues,"descend");
  
  //compute the covariance matrix
  int iter = 0;
  int max_iter = 10;
  int j = rmax+1;
  while(iter < max_iter){
    
    // Subset the eigenvalues vector
    arma::vec subset_eigenvalues = sorted_eigenvalues(arma::span(j-1, j+3));
    
    // Create the sequence (j-1)^(2/3), ..., (j+3)^(2/3)
    arma::vec sequence = arma::linspace<arma::vec>(j - 1, j + 3, 5);
    sequence = arma::pow(sequence, 2.0 / 3.0);
    
    // Create the design matrix X with a column of ones and the sequence
    arma::mat X(sequence.n_elem, 2);
    X.col(0).ones();
    X.col(1) = sequence;
    
    // Perform OLS regression to obtain the coefficients
    arma::vec coefficients = arma::solve(X, subset_eigenvalues);
    
    //Compute the betahat
    double delta = 2 * std::abs(coefficients(1));
    
    // Subset rmax biggest eigenvalues from 0 to rmax+1
    arma::vec subset_sorted_eigenvalues = sorted_eigenvalues(arma::span(0, rmax));
    // Create vector of differences
    arma::vec diff_eigenvalues = subset_sorted_eigenvalues(arma::span(0, rmax-1)) - subset_sorted_eigenvalues(arma::span(1, rmax));

    //subset eigenvalues greater than delta
    //arma::vec subset_eigenvalues_delta = subset_eigenvalues.elem(find(subset_eigenvalues > delta));
    arma::uvec indices = find(diff_eigenvalues > delta);
    int r_hat = 0;
    if (indices.is_empty()) {
      r_hat = 0;
    }
    else{
     // Subset r_hat 
      r_hat = arma::max(indices) + 1;
    }
    //set j
    int new_j = r_hat + 1;
    
    if(new_j == j){
      cout<<"Convergence reached after "<<iter<<" iterations"<<endl;
      break;
    }
    else
    {
      j = new_j;
      iter = iter+1;
    }
  }    
  return j;
}


