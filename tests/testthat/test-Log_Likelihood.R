context("Log Likelihood Functions")

test_that("negloglik_vecchia_ST", {
    set.seed(7919)
    load("small_sim.Rdata")
    return(1) # this test isn't useful at the moment
    params <- c(0.5, 0.5, 0.5, 0.5)
    locs <- locs_train
    X <- X_train
    Y <- Y_train
    locs.scaled <- locs / c(params[2], params[2], params[3])
    vecchia_approx <- vecchia_specify(locs.scaled, 25)
    beta.hat <- rep(0, ncol(X))
    Y.hat <- as.matrix(X) %*% beta.hat
    res <- as.double(Y - Y.hat)
    ord_locs <- locs[vecchia_approx$ord, ]
    locs_mat <- cbind(ord_locs[, 1], ord_locs[, 2], ord_locs[, 3])
    result <- negloglik_vecchia_ST(log(params), locs_mat, res, vecchia_approx)
    expect_equal(1597.9808688, result, tolerance = 10e-5)
    vecchia.result <- optim(
        par = log(params), fn = negloglik_vecchia_ST,
        locs = locs_mat, res = res, vecchia.approx = vecchia_approx, method = "Nelder-Mead",
        control = list(trace = 0)
    )
    params_optim <- exp(vecchia.result$par)
    result_optim <- negloglik_vecchia_ST(log(params_optim), locs_mat, res, vecchia_approx)
    expect_equal(-6861.65048095, result_optim, tolerance = 1e-5)
})

test_that("negloglik.full", {
    set.seed(1234)

    y.orig <- rnorm(100)

    locs.dim <- seq(1, 10, length = sqrt(length(y.orig)))
    locs <- as.matrix(expand.grid(locs.dim, locs.dim))

    covmat.true <- Matern(rdist(locs), range = 1, smoothness = 0.5) + diag(100)

    y <- y.orig %*% chol(covmat.true)
    y <- as.vector(y)

    params.init <- rep(NA, 4)
    params.init[1] <- 0.9 * var(y)
    params.init[2] <- 1
    params.init[3] <- 0.5
    params.init[4] <- 0.1 * var(y)

    params.init <- create.initial.values.flex(
        params.init[1],
        params.init[2],
        params.init[3],
        params.init[4],
        1, 1
    )

    d <- rdist(locs)

    res.optim.NM <- optim(
        par = params.init, fn = negloglik.full, locs = locs, y = y,
        control = list(maxit = 5000)
    )

    LL.full <- negloglik.full(res.optim.NM$par, locs, y)

    params.final <- c(
        exp(res.optim.NM$par[1:2]),
        gtools::inv.logit(res.optim.NM$par[3], 0, 2.5),
        exp(res.optim.NM$par[4])
    )

    pgp.params <- create.initial.values.flex(
        params.final[1],
        params.final[2],
        params.final[3],
        params.final[4],
        1, 1
    )
    pseq <- create.param.sequence(1)

    LL.full.pgp <- mvnegloglik.full(pgp.params, list(locs), y, pseq)

    vec.approx <- vecchia_specify(locs, 99)

    LL.vecchia <- -1 * vecchia_likelihood(
        y, vec.approx, params.final[1:3],
        params.final[4]
    )

    vec.approx.pgp <- vecchia_Mspecify(list(locs), 99)
    vec.U.pgp <- createUMultivariate(vec.approx.pgp, c(params.final, 1))

    LL.pgp <- -1 * GPvecchia:::vecchia_likelihood_U(y, vec.U.pgp)

    expect_equal(173.315, LL.full, tolerance = 1e-3)
    # Univariate likelihood should equal the multivariate likelihood
    expect_equal(LL.full, LL.full.pgp, tolerance = 1e-3)
    # Full likelihood should equal both the univariate and multivariate
    # Vecchia approximations
    expect_equal(LL.full, LL.vecchia, tolerance = 1e-3)
    expect_equal(LL.full, LL.pgp, tolerance = 1e-3)
})

test_that("mvnegloglik", {
    source("sim_multivariate_big.R")
    P <- 3
    Y <- cbind(runif(10), runif(10), runif(10))
    cor.matrix <- cor(Y)
    cov_mat <- c(cor.matrix[upper.tri(cor.matrix)])
    logparams <- create.initial.values.flex(
        rep(0.9, P), # marginal variance
        rep(0.5, P), # range
        rep(0.5, P), # smoothness
        rep(0.1, P), # nuggets
        cov_mat,
        P
    )
    pseq <- create.param.sequence(P)
    vec.approx <- vecchia_Mspecify(locs.list, 25)
    neg_likelihood <- mvnegloglik(
        logparams, vec.approx,
        unlist(y.list), pseq, P
    )
    expect_equal(34870.57, neg_likelihood, tolerance = 1e-2)
})

test_that("mvnegloglik_ST", {
    source("sim_multivariate_big_st.R")
    P <- 3
    Y <- cbind(runif(10), runif(10), runif(10))
    cor.matrix <- cor(Y)
    cov_mat <- c(cor.matrix[upper.tri(cor.matrix)])
    logparams <- create.initial.values.flex(
        rep(0.9, P), # marginal variance
        rep(c(2, 3), P), # range
        rep(0.5, P), # smoothness
        rep(0.1, P), # nuggets
        cov_mat,
        P
    )
    pseq <- create.param.sequence(P, 2)
    vec.approx <- vecchia_Mspecify(locs.list, 25)
    neg_likelihood <- mvnegloglik_ST(
        logparams, vec.approx,
        unlist(y.list), pseq, P, c(1, 1, 2), 2
    )
    expect_equal(35106.73, neg_likelihood, tolerance = 1e-2)

    vec.approx2 <- vec.approx
    for (i in 1:P) {
        vec.approx2$locsord[, 1:2] <- vec.approx$locsord[, 1:2] / 2
        vec.approx2$locsord[, 3] <- vec.approx$locsord[, 3] / 3
    }
    logparams2 <- logparams
    logparams2[pseq[2, 1]:pseq[2, 2]] <- 0
    neg_likelihood2 <- mvnegloglik_ST(
        logparams2, vec.approx2,
        unlist(y.list), pseq, P, c(1, 1, 2), 2
    )
    expect_equal(neg_likelihood, neg_likelihood2, tolerance = 1e-3)

    logparams <- create.initial.values.flex(
        rep(0.9, P), # marginal variance
        2:7, # range
        rep(0.5, P), # smoothness
        rep(0.1, P), # nuggets
        cov_mat,
        P
    )
    neg_likelihood <- mvnegloglik_ST(
        logparams, vec.approx,
        unlist(y.list), pseq, P, c(1, 1, 2), 2
    )
    expect_equal(36354.9, neg_likelihood, tolerance = 1e-2)

    vec.approx2 <- vec.approx
    vec.approx2$locsord[vec.approx$ondx == 1, 1:2] <-
        vec.approx$locsord[vec.approx$ondx == 1, 1:2] / 2
    vec.approx2$locsord[vec.approx$ondx == 1, 3] <-
        vec.approx$locsord[vec.approx$ondx == 1, 3] / 3
    vec.approx2$locsord[vec.approx$ondx == 2, 1:2] <-
        vec.approx$locsord[vec.approx$ondx == 2, 1:2] / 4
    vec.approx2$locsord[vec.approx$ondx == 2, 3] <-
        vec.approx$locsord[vec.approx$ondx == 2, 3] / 5
    vec.approx2$locsord[vec.approx$ondx == 3, 1:2] <-
        vec.approx$locsord[vec.approx$ondx == 3, 1:2] / 6
    vec.approx2$locsord[vec.approx$ondx == 3, 3] <-
        vec.approx$locsord[vec.approx$ondx == 3, 3] / 7

    neg_likelihood2 <- mvnegloglik_ST(
        logparams2, vec.approx2,
        unlist(y.list), pseq, P, c(1, 1, 2), 2
    )
    expect_equal(neg_likelihood, neg_likelihood2, tolerance = 1e-3)
})

test_that("mvnegloglik.full", {
    source("sim_multivariate_lik.R")

    pseq <- create.param.sequence(3)

    param.marg.var <- 0.9 * unlist(lapply(y.list, var))
    param.marg.scale <- rep(1, 3)
    param.marg.smooth <- rep(0.5, 3)
    param.marg.nugget <- 0.1 * unlist(lapply(y.list, var))
    param.rho <- rep(0.5, 3)
    params.init <- create.initial.values.flex(
        param.marg.var,
        param.marg.scale,
        param.marg.smooth,
        param.marg.nugget,
        param.rho, 3
    )

    res.optim.NM <- optim(
        par = params.init, fn = mvnegloglik.full,
        locs = locs.list, y = unlist(y.list),
        param.seq = pseq,
        method = "Nelder-Mead",
        control = list(trace = 0, maxit = 10000, reltol = 1e-4)
    )

    LL.full.mv <- mvnegloglik.full(
        res.optim.NM$par, locs.list,
        unlist(y.list), pseq
    )

    param.seq.begin <- pseq[, 1]
    param.seq.end <- pseq[, 2]
    params.init.final <- res.optim.NM$par
    params.init.final.t <- c(
        exp(params.init.final[param.seq.begin[1]:
        param.seq.end[2]]),
        gtools::inv.logit(params.init.final[
            param.seq.begin[3]:
            param.seq.end[3]
        ], 0, 2.5),
        exp(params.init.final[param.seq.begin[4]:
        param.seq.end[4]]),
        tanh(params.init.final[
            param.seq.begin[5]:param.seq.end[5]
        ])
    )

    params.final.flex.test <- create.initial.values.flex(
        params.init.final.t[param.seq.begin[1]:param.seq.end[1]],
        params.init.final.t[param.seq.begin[2]:param.seq.end[2]],
        params.init.final.t[param.seq.begin[3]:param.seq.end[3]],
        params.init.final.t[param.seq.begin[4]:param.seq.end[4]],
        params.init.final.t[param.seq.begin[5]:param.seq.end[5]], 3
    )

    cov.list <- create.cov.upper.flex(
        3,
        params.init.final.t[param.seq.begin[1]:param.seq.end[1]],
        params.init.final.t[param.seq.begin[2]:param.seq.end[2]],
        params.init.final.t[param.seq.begin[3]:param.seq.end[3]],
        params.init.final.t[param.seq.begin[4]:param.seq.end[4]],
        params.init.final.t[param.seq.begin[5]:param.seq.end[5]]
    )

    cov.mat <- cat.covariances(
        locs.list, cov.list$variance,
        cov.list$range,
        cov.list$smoothness, cov.list$nugget
    )

    LL.full.calc <- -1 * mvtnorm::dmvnorm(unlist(y.list),
        rep(0, length(unlist(y.list))),
        cov.mat,
        log = TRUE
    )

    vec.mapprox <- vecchia_Mspecify(locs.list, length(unlist(y.list)) - 1)
    U.mobj <- createUMultivariate(vec.mapprox, params.init.final.t)

    LL.vecchia.mv <- -1 * GPvecchia:::vecchia_likelihood_U(unlist(y.list), U.mobj)

    expect_equal(541.31, LL.full.mv, tolerance = 1e-3)
    expect_equal(LL.full.calc, LL.full.mv, tolerance = 1e-3)
    # Full likelihood should equal the Vecchia likelihood
    expect_equal(LL.full.mv, LL.vecchia.mv, tolerance = 1e-3)
    # Check create.initial.values.flex:
    expect_equal(params.init.final, params.final.flex.test, tolerance = 1e-3)
    # Check create.cov.upper.flex:
    cov.list.true <- readRDS("covlist.rds")
    expect_equal(cov.list$variance, cov.list.true$variance, tolerance = 1e-3)
    expect_equal(cov.list$range, cov.list.true$range, tolerance = 1e-3)
    expect_equal(cov.list$smoothness, cov.list.true$smoothness, tolerance = 1e-3)
    expect_equal(cov.list$nugget, cov.list.true$nugget, tolerance = 1e-3)
    # Check cat.covariances:
    cov.mat.true <- readRDS("covmat.rds")
    expect_equal(cov.mat, cov.mat.true, tolerance = 1e-3)
})
