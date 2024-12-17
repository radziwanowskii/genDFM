test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("basic operations work", {
  expect_no_error(determine_number_of_dynamic_factors(data = scale(MacroPL[,-1]), method = "on", max_factors = 10))
  expect_no_error(determine_number_of_dynamic_factors(data = scale(MacroPL[,-1]), method = "av", max_factors = 10))
  expect_no_error(determine_number_of_dynamic_factors(data = scale(MacroPL[,-1]), method = "hl", max_factors = 10))

  expect_no_error(determine_number_of_static_factors(data = scale(MacroPL[,-1]), method = "on", max_factors = 10))
  expect_no_error(determine_number_of_static_factors(data = scale(MacroPL[,-1]), method = "al", max_factors = 10))
  expect_no_error(determine_number_of_static_factors(data = scale(MacroPL[,-1]), method = "ah", max_factors = 10))
  expect_no_error(determine_number_of_static_factors(data = scale(MacroPL[,-1]), method = "bn", max_factors = 10))

})

test_that("wrong inputs", {
  df_with_na <- scale(MacroPL[,-1])
  df_with_na[1,1] <- NA
  expect_error(determine_number_of_dynamic_factors(data = scale(MacroPL[,-1]), method = "ol", max_factors = 10),"Invalid method")
  expect_error(determine_number_of_static_factors(data = scale(MacroPL[,-1]), method = "ol", max_factors = 10),"Invalid method")

  expect_error(determine_number_of_dynamic_factors(data = scale(MacroPL), method = "on", max_factors = 10),"numeric")
  expect_error(determine_number_of_static_factors(data = scale(MacroPL), method = "on", max_factors = 10),"numeric")

  expect_error(determine_number_of_dynamic_factors(data = df_with_na, method = "on", max_factors = 10),"missing values")
  expect_error(determine_number_of_dynamic_factors(data = df_with_na, method = "on", max_factors = 10),"missing values")
})
