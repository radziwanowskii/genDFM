test_that("corect classes are returned", {
  data(MacroPL)
  dfm1 <- estimate_DFM(scale(MacroPL[,-1]), model = "fhlr_2000")
  expect_s3_class(dfm1, "generalized_dfm")
  dfm2 <- estimate_DFM(scale(MacroPL[,-1]), r = 5, model = "fhlr_2005")
  expect_s3_class(dfm2, "restricted_gdfm")
  dfm3 <- estimate_DFM(scale(MacroPL[,-1]), model = "sw_2002")
  expect_s3_class(dfm3, "sw_2002")
  dfm4 <- estimate_DFM(scale(MacroPL[,-1]), model = "doz_2011")
  expect_s3_class(dfm4, "dfm")
  dfm5 <- estimate_DFM(scale(MacroPL[,-1]), model = "doz_2012")
  expect_s3_class(dfm5, "dfm")
  dfm6 <- estimate_DFM(scale(MacroPL[,-1]), model = "BM_2014")
  expect_s3_class(dfm6, "dfm")
})

test_that("misspecified models",  {
  data(MacroPL)
  expect_error(estimate_DFM(scale(MacroPL[,-1]), model = "xxx"),"Invalid model specified")
  expect_error(determine_number_of_static_factors(scale(MacroPL[,-1]),"xxx"), "Invalid method")
  expect_error(determine_number_of_dynamic_factors(scale(MacroPL[,-1]),"xxx"), "Invalid method")
})

test_that("correct inputs", {
  data(MacroPL)
  df <- data.frame(MacroPL[,-1])
  expect_no_error(estimate_DFM(df, model = "fhlr_2000"))
  expect_no_error(estimate_DFM(df, model = "fhlr_2005"))
  expect_no_error(estimate_DFM(df, model = "sw_2002"))
  expect_no_error(estimate_DFM(df, model = "doz_2011"))
  expect_no_error(estimate_DFM(df, model = "doz_2012"))
  expect_no_error(estimate_DFM(df, model = "BM_2014"))

  expect_equal(estimate_DFM(df, model = "fhlr_2000"), estimate_DFM(MacroPL[,-1], model = "fhlr_2000"))
  expect_equal(estimate_DFM(df, model = "fhlr_2005"), estimate_DFM(MacroPL[,-1], model = "fhlr_2005"))
  expect_equal(estimate_DFM(df, model = "sw_2002"), estimate_DFM(MacroPL[,-1], model = "sw_2002"))
  expect_equal(estimate_DFM(df, model = "doz_2011"), estimate_DFM(MacroPL[,-1], model = "doz_2011"))
  expect_equal(estimate_DFM(df, model = "doz_2012"), estimate_DFM(MacroPL[,-1], model = "doz_2012"))
  expect_equal(estimate_DFM(df, model = "BM_2014"), estimate_DFM(MacroPL[,-1], model = "BM_2014"))
})

test_that("incorrect inputs", {
  data(MacroPL)
  data_with_NA <- data.frame(MacroPL[,-1])
  data_with_NA[1,1] <- NA

  expect_error(estimate_DFM(MacroPL,model = "fhlr_2000"),"numeric")
  expect_error(estimate_DFM(MacroPL,model = "fhlr_2005"),"numeric")
  expect_error(estimate_DFM(MacroPL,model = "sw_2002"),"numeric")
  expect_error(estimate_DFM(MacroPL,model = "doz_2011"),"numeric")
  expect_error(estimate_DFM(MacroPL,model = "doz_2012"),"numeric")
  expect_error(estimate_DFM(MacroPL,model = "BM_2014"),"numeric")

  expect_error(estimate_DFM(data_with_NA, model = "fhlr_2000"), "missing values")
  expect_error(estimate_DFM(data_with_NA, model = "fhlr_2005"), "missing values")
  expect_error(estimate_DFM(data_with_NA, model = "sw_2002"), "missing values")
  expect_error(estimate_DFM(data_with_NA, model = "doz_2011"), "missing values")
  expect_error(estimate_DFM(data_with_NA, model = "doz_2012"), "missing values")
  expect_error(estimate_DFM(data_with_NA, model = "BM_2014"), "missing values")
})

test_that("dimensions work",{
  data(MacroPL)
  expect_error(estimate_DFM(MacroPL[,-1],model = "fhlr_2000",q = c(1,2)))
  expect_error(estimate_DFM(MacroPL[,-1],model = "fhlr_2005",r = c(1,2)))
  expect_error(estimate_DFM(MacroPL[,-1],model = "sw_2002",r = c(1,2)))
  expect_error(estimate_DFM(MacroPL[,-1],model = "doz_2011",r = c(1,2)))
  expect_error(estimate_DFM(MacroPL[,-1],model = "doz_2012",r = c(1,2)))
  expect_error(estimate_DFM(MacroPL[,-1],model = "BM_2014",r = c(1,2)))

  expect_error(estimate_DFM(MacroPL[,2],model = "fhlr_2000"))
  expect_error(estimate_DFM(MacroPL[,2],model = "fhlr_2005"))
  expect_error(estimate_DFM(MacroPL[,2],model = "sw_2002"))
  expect_error(estimate_DFM(MacroPL[,2],model = "doz_2011"))
  expect_error(estimate_DFM(MacroPL[,2],model = "doz_2012"))
  expect_error(estimate_DFM(MacroPL[,2],model = "BM_2014"))
})

test_that("predict_dfm",{
  dfm <- estimate_DFM(scale(MacroPL[,-1]),model = 'fhlr_2005')
  expect_no_error(predict.DFM(dfm,3))
})

