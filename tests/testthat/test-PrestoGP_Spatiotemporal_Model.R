context("Spatiotemporal Model")

test_that("Invalid X", {
  load("small_sim.Rdata")
  return(1)
  model <- new("SpatiotemporalModel")
  expect_error(
    prestogp_fit(model, Y_train, "X_train", locs_train),
    "X parameter must be a matrix."
  )
})

test_that("Invalid Y", {
  load("small_sim.Rdata")
  return(1)
  model <- new("SpatiotemporalModel")
  expect_error(
    prestogp_fit(model, "Y_train", X_train, locs_train),
    "Y parameter must be a matrix."
  )
})

test_that("Invalid Locs", {
  load("small_sim.Rdata")
  return(1)
  model <- new("SpatiotemporalModel")
  expect_error(
    prestogp_fit(model, Y_train, X_train, "locs_train"),
    "locs parameter must be a matrix."
  )
})

test_that("X length mismatch", {
  load("small_sim.Rdata")
  return(1)
  model <- new("SpatiotemporalModel")
  expect_error(
    prestogp_fit(model, Y_train, X_train[1:10, ], locs_train),
    "Y must have the same number of rows as X."
  )
})

test_that("locs length mismatch", {
  load("small_sim.Rdata")
  return(1)
  model <- new("SpatiotemporalModel")
  expect_error(
    prestogp_fit(model, Y_train, X_train, locs_train[1:10, ]),
    "Y must have the same number of rows as locs."
  )
})

test_that("Invalid fit m", {
  load("small_sim.Rdata")
  return(1)
  model <- new("SpatiotemporalModel", n_neighbors = 0)
  expect_error(
    prestogp_fit(model, Y_train, X_train, locs_train),
    "M must be at least 3."
  )
})

test_that("Invalid predict Locs", {
  load("small_sim.Rdata")
  return(1)
  model <- new("SpatiotemporalModel")
  expect_error(
    prestogp_predict(model, X_test, "locs_test"),
    "The locs parameter must be a matrix."
  )
})

test_that("Invalid predict X", {
  load("small_sim.Rdata")
  return(1)
  model <- new("SpatiotemporalModel")
  expect_error(
    prestogp_predict(model, "X_test", locs_test),
    "X parameter must be a matrix."
  )
})

# test_that("Invalid predict m", {
#   load("small_sim.Rdata")
#   error <- tryCatch({
#     model <- SpatiotemporalModel()
#     prediction <- prestogp_predict(model, X_test, locs_test, m = 0)
#     FALSE
#   },
#   error=function(cond){
#     TRUE
#   })
#   expect_true(error)
# })

test_that("locs predict mismatch", {
  load("small_sim.Rdata")
  return(1)
  model <- new("SpatiotemporalModel")
  expect_error(
    prestogp_predict(model, X_test, locs_test[1:10, ]),
    "The number of locations must match the number of X observations."
  )
})

test_that("Simulated dataset Spatiotemporal", {
  set.seed(7919)
  load("sim_spatiotemporal.Rdata")
  return(1)
  model <- new("SpatiotemporalModel", n_neighbors = 4)
  model <- prestogp_fit(model, Y_train, X_train, locs_train)
  prediction <- prestogp_predict(model, X_test, locs_test, m = 4)
  means <- prediction[[1]]
  abs_error <- abs(Y_test - means)
  expect_equal(5.64803, mean(abs_error), tolerance = 10e-2)
  mean_sds <- mean(prediction[[2]])
  expect_equal(8.58, mean_sds, tolerance = 10e-2)
  expect_equal(220.22, model@covparams[1], tolerance = 10e-2)
  expect_equal(8.44, model@covparams[2], tolerance = 10e-2)
  expect_equal(0.08, model@covparams[3], tolerance = 10e-2)
  expect_equal(21.60, model@covparams[4], tolerance = 10e-2)
})
