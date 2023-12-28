context("Spatiotemporal Full Model")

test_that("Simulated dataset spatiotemporal full", {
  set.seed(7919)
  load("sim_spatiotemporal.Rdata")
  return(1)
  model <- SpatiotemporalFullModel()
  model <- prestogp_fit(model, Y_train, X_train, locs_train)
  prediction <- prestogp_predict(model, X_test, locs_test)
  means <- prediction[[1]]
  abs_error <- abs(Y_test - means)
  expect_equal(5.68, mean(abs_error), tolerance = 10e-2)
  expect_equal(236.25, model@covparams[1], tolerance = 10e-2)
  expect_equal(7.51, model@covparams[2], tolerance = 10e-2)
  expect_equal(0.08, model@covparams[3], tolerance = 10e-2)
  expect_equal(20.29, model@covparams[4], tolerance = 10e-2)
})
