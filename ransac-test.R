# please run after running the scripts 'topdown-ransac-utils.R' and
# 'topdown-ransac-example.R'

library(testthat)

test_that("ransaclm deals with problematic inputs", {
  data_with_NA <- data_simple
  data_with_NA[101:110, ] <- NA
  expect_warning(ransaclm(y ~ . - inlier,
                          data = data_with_NA, error_threshold = 2,
                          inlier_threshold = 50, seed = 20171111
  ))
  expect_error(ransaclm(y ~ . - inlier,
                        data = 1:5, error_threshold = 2,
                        inlier_threshold = 50, seed = 20171111
  ))
  data_n_smaller_p <- data_simple
  data_n_smaller_p["a"] <- 1:100
  expect_error(ransaclm(y ~ . - inlier,
                        data = data_n_smaller_p, error_threshold = 2,
                        inlier_threshold = 50, seed = 20171111,
                        minimum_observations = 1
  ))
  data_factor <- data_simple
  data_factor["fac"] <- as.factor(1:100)
  expect_error(ransaclm(y ~ . - inlier,
                        data = data_factor, error_threshold = 2,
                        inlier_threshold = 50, seed = 20171111
  ))
  data_only_NA <- data.frame(y = rep(NA, 100),
                             x = rep(NA, 100),
                             inlier = rep(NA, 100))
  expect_error(ransaclm(y ~ . - inlier,
                          data = data_only_NA, error_threshold = 2,
                          inlier_threshold = 50, seed = 20171111
  ))
})


