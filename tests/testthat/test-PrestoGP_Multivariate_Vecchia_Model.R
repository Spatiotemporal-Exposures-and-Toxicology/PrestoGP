context("Multivariate Vecchia Model")

test_that("Invalid locs input", {
    model <- new("MultivariateVecchiaModel")
    expect_error(
        prestogp_fit(model, as.matrix(1:3), as.matrix(1:3), 1:3),
        "locs parameter must be a matrix or a list."
    )
})

test_that("Simulated dataset multivariate spatial", {
    source("sim_multivariate_big.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 25)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list, X.st, locs.list,
        scaling = c(1, 1), apanasovich = TRUE, verbose = FALSE,
        optim.control = list(
            trace = 0, maxit = 5000,
            reltol = 1e-3
        )
    )
    beta.out <- as.vector(pgp.mmodel1@beta)
    params.out <- pgp.mmodel1@covparams

    expect_length(beta.out, 31)
    expect_length(params.out, 15)
    expect_equal(beta.out, c(
        0, 0.95, 0.93, 0.92, 0.88, rep(0, 6), 0.56,
        0.7, 1.12, 1, rep(0, 6), 0.94, 0.85,
        1.03, 0.94, rep(0, 6)
    ), tolerance = 0.06)
    expect_equal(params.out[1], 1.6, tolerance = 2.3)
    expect_equal(params.out[2], 2.9, tolerance = 4)
    expect_equal(params.out[3], 3.9, tolerance = 3.7)
    expect_equal(params.out[4], 0.41, tolerance = 2.4)
    expect_equal(params.out[5], 0.35, tolerance = 1.9)
    expect_equal(params.out[6], 0.41, tolerance = 1)
    expect_equal(params.out[7], 0.61, tolerance = 1.2)
    expect_equal(params.out[8], 0.43, tolerance = 0.9)
    expect_equal(params.out[9], 0.87, tolerance = 1.5)
    expect_equal(params.out[10], 1.8, tolerance = 1.6)
    expect_equal(params.out[11], 2.4, tolerance = 2.2)
    expect_equal(params.out[12], 1.4, tolerance = 1.2)
    expect_equal(params.out[13], 0.17, tolerance = 0.6)
    expect_equal(params.out[14], 0.33, tolerance = 0.4)
    expect_equal(params.out[15], 0.14, tolerance = 0.4)
})

test_that("Simulated dataset multivariate spatiotemporal", {
    source("sim_multivariate_big_st.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 25)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list, X.st, locs.list,
        scaling = c(1, 1, 2), verbose = FALSE,
        optim.control = list(
            trace = 0, maxit = 5000,
            reltol = 1e-3
        )
    )
    beta.out <- as.vector(pgp.mmodel1@beta)
    params.out <- pgp.mmodel1@covparams

    expect_length(beta.out, 31)
    expect_length(params.out, 18)
    expect_equal(beta.out, c(
        0, 0.91, 0.86, 0.82, 0.97, rep(0, 6), 0.95,
        0.97, 0.92, 0.78, rep(0, 6), 0.8, 0.97,
        1.04, 0.81, rep(0, 6)
    ), tolerance = 1.1)
})

test_that("Simulated dataset multivariate spatial prediction", {
    source("sim_multivariate_big_pred.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 25)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
        locs.list.otr,
        scaling = c(1, 1),
        apanasovich = TRUE, verbose = FALSE,
        optim.control = list(
            trace = 0, maxit = 5000,
            reltol = 1e-3
        )
    )

    pgp.mmodel1.pred <- prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst)

    mse <- mean((pgp.mmodel1.pred$means - unlist(y.list.otst))^2)

    expect_equal(mse, 1.99, tolerance = 0.1)
})
