
test_that("wrong inputs",{
  data(MacroPL)
  data_with_NA <- MacroPL[,-1]
  data_with_NA[1,1] <- NA

  expect_error(FAVAR(data_with_NA, r = 3, n_lags = 13, slow_indices = c(1:22), direct_indices = 1:3), "missing values")
  expect_error(FAVAR(MacroPL, r = 3, n_lags = 13, slow_indices = c(1:22), direct_indices = 1:3), "numeric")
  expect_error(FAVAR(MacroPL[,-1], r = 3, n_lags = 13, direct_indices = 1:3), "Slow indices")
  expect_error(FAVAR(MacroPL[,-1], r = 3, n_lags = 13, slow_indices = 1:3), "Direct indices")
})

test_that("wrong dimensions", {
  data(MacroPL)
  expect_error(FAVAR(MacroPL[,-1], r = c(1,2), n_lags = 13, slow_indices = c(1:22), direct_indices = 1:3))
  expect_error(FAVAR(MacroPL[,-1], r = 3, n_lags = c(1,3), slow_indices = c(1:22), direct_indices = 1:3))
  expect_error(FAVAR(MacroPL[,2], r = 3, n_lags = 13, slow_indices = c(1:22), direct_indices = 1:3))
  
  expect_error(FAVAR(MacroPL[,-1], r = 3, n_lags = 13, slow_indices = c(1:22), direct_indices = 1), "Direct indices")
  expect_error(FAVAR(MacroPL[,-1], r = 3, n_lags = 13, slow_indices = 1, direct_indices = 1:3), "Slow indices")
})

test_that("irf function", {
  data(MacroPL)
  favar <- FAVAR(scale(MacroPL[,-1]), r = 3, n_lags = 13, slow_indices = c(1:22), direct_indices = c("Ref","Lom"))

  expect_no_error(favar_IRF(favar, shock_variable = "Ref", response_variable = "Ref", n_ahead = 49, shock = "unit", plot = FALSE))
  expect_no_error(favar_IRF(favar, shock_variable = "Ref", response_variable = "Ref", n_ahead = 49, shock = "sd", plot = FALSE))
  expect_error(favar_IRF(favar, shock_variable = "Ref", response_variable = c("Ref","CPI"), n_ahead = 49, plot = FALSE))
   expect_error(favar_IRF(favar, shock_variable = c("Ref","CPI"), response_variable = c("Ref"), n_ahead = 49, plot = FALSE))

  expect_no_error(favar_IRF(favar, shock_variable = "Ref", response_variable = "Ref", n_ahead = 49, shock = "unit", plot = TRUE))
})




