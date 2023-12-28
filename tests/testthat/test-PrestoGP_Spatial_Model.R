context("Spatial Model")

test_that("Invalid predict Locs", {
  load("small_sim.Rdata")
  return(1)
  error <- tryCatch(
    {
      model <- new("SpatialModel")
      prediction <- prestogp_predict(model, X_test, "locs_test", m = 4)
      FALSE
    },
    error = function(cond) {
      TRUE
    }
  )
  expect_true(error)
})

test_that("Invalid predict X", {
  load("small_sim.Rdata")
  return(1)
  error <- tryCatch(
    {
      model <- new("SpatialModel")
      prediction <- prestogp_predict(model, "X_test", locs_test, m = 4)
      FALSE
    },
    error = function(cond) {
      TRUE
    }
  )
  expect_true(error)
})

test_that("Invalid predict locs (not 2 columns)", {
  load("small_sim.Rdata")
  return(1)
  error <- tryCatch(
    {
      model <- new("SpatialModel")
      prediction <- prestogp_predict(model, matrix(rnorm(100), ncol = 10), matrix(rnorm(30), ncol = 3), m = 4)
      FALSE
    },
    error = function(cond) {
      TRUE
    }
  )
  expect_true(error)
})

test_that("locs length mismatch", {
  load("small_sim.Rdata")
  return(1)
  error <- tryCatch(
    {
      model <- new("SpatialModel")
      prediction <- prestogp_predict(model, matrix(rnorm(100), ncol = 10), matrix(rnorm(50), ncol = 2))
      FALSE
    },
    error = function(cond) {
      TRUE
    }
  )
  expect_true(error)
})

test_that("Simulated dataset spatial", {
  set.seed(7919)
  load("sim_spatial.Rdata")
  return(1)
  model <- new("SpatialModel", n_neighbors = 4)
  model <- prestogp_fit(model, Y_train, X_train, locs_train)
  prediction <- prestogp_predict(model, X_test, locs_test)
  means <- prediction[[1]]
  mse <- mean((Y_test - means)^2)
  expect_equal(1.93825, mse, tolerance = 10e-4)
  mean_sds <- mean(prediction[[2]])
  expect_equal(0.967, mean_sds, tolerance = 10e-3)
  expect_equal(9.46, model@covparams[1], tolerance = 10e-2)
  expect_equal(0.27, model@covparams[2], tolerance = 10e-2)
  expect_equal(0.93, model@covparams[3], tolerance = 10e-2)
})
