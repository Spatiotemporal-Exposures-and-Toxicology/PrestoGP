context("Spatial Full Model")

test_that("Simulated dataset spatial full", {
  set.seed(7919)
  load("sim_spatial.Rdata")
  return(1)
  model <- FullSpatialModel()
  model <- prestogp_fit(model, Y_train, X_train, locs_train)
  prediction <- prestogp_predict(model, X_test, locs_test)
  means <- prediction[[1]]
  mse <- mean((Y_test - means)^2)
  expect_equal(1.94, mse, tolerance = 10e-2)
  expect_equal(9.46, model@covparams[1], tolerance = 10e-2)
  expect_equal(0.27, model@covparams[2], tolerance = 10e-2)
  expect_equal(0.93, model@covparams[3], tolerance = 10e-2)
})
