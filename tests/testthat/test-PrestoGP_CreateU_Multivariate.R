context("createUMultivariate")

test_that("create.param.sequence", {
    seq <- create.param.sequence(1)
    colnames(seq) <- NULL
    expect_equal(2, ncol(seq))
    expect_equal(5, nrow(seq))
    expect_equal(c(1, 1), seq[1, ])
    expect_equal(c(2, 2), seq[2, ])
    expect_equal(c(3, 3), seq[3, ])
    expect_equal(c(4, 4), seq[4, ])
    expect_equal(c(5, 5), seq[5, ])

    seq <- create.param.sequence(3)
    colnames(seq) <- NULL
    expect_equal(2, ncol(seq))
    expect_equal(5, nrow(seq))
    expect_equal(c(1, 3), seq[1, ])
    expect_equal(c(4, 6), seq[2, ])
    expect_equal(c(7, 9), seq[3, ])
    expect_equal(c(10, 12), seq[4, ])
    expect_equal(c(13, 15), seq[5, ])

    seq <- create.param.sequence(3, 2)
    colnames(seq) <- NULL
    expect_equal(2, ncol(seq))
    expect_equal(5, nrow(seq))
    expect_equal(c(1, 3), seq[1, ])
    expect_equal(c(4, 9), seq[2, ])
    expect_equal(c(10, 12), seq[3, ])
    expect_equal(c(13, 15), seq[4, ])
    expect_equal(c(16, 18), seq[5, ])
})

test_that("max_min_ordering", {
    set.seed(7919)
    load("multivariate_sim_spatial3.Rdata")
    order <- max_min_ordering(locs_train, fields::rdist)
    exp_order <- c(
        63, 21, 50, 30, 76, 78, 40, 36, 23, 69, 9, 67, 32, 2, 20, 62,
        22, 31, 74, 39, 35, 58, 68, 54, 41, 3, 80, 46, 6, 88, 12, 47,
        72, 42, 13, 83, 25, 52, 11, 60, 24, 1, 28, 84, 29, 64, 66, 81,
        82, 55, 61, 87, 17, 33, 43, 45, 10, 79, 53, 75, 89, 51, 73, 27,
        26, 77, 44, 38, 65, 16, 19, 37, 57, 70, 15, 4, 5, 86, 14, 49,
        85, 34, 48, 59, 18, 8, 71, 7, 56, 90
    )
    order.gpv <- order_maxmin_exact(locs_train)
    expect_equal(exp_order, order)
    expect_equal(order, order.gpv)
})

test_that("sparseNN", {
    set.seed(1212)

    locs <- matrix(nrow = 100, ncol = 2)
    locs[1, ] <- rep(0, 2)
    for (i in 2:nrow(locs)) {
        cur.r <- rnorm(1, 5)
        cur.t <- runif(1, 0, 2 * pi)
        locs[i, ] <- locs[i - 1, ] + c(cur.r * cos(cur.t), cur.r * sin(cur.t))
    }

    mm.order <- order_maxmin_exact(locs)
    olocs <- locs[mm.order, ]
    pgp.nn <- sparseNN(olocs, 5, fields::rdist, "rdist")
    gpv.nn <- GpGp:::find_ordered_nn(olocs, 5)

    indices <- matrix(nrow = nrow(olocs), ncol = 5)
    distances <- indices
    for (i in 1:nrow(olocs)) {
        if (i <= 5) {
            cur.dist <- fields::rdist(
                olocs[(1:(5 + 1)), ][-i, ],
                olocs[i, , drop = FALSE]
            )
            indices[i, ] <- order(cur.dist)
        } else {
            cur.dist <- fields::rdist(olocs[(1:(i - 1)), ], olocs[i, , drop = FALSE])
            indices[i, ] <- order(cur.dist)[1:5]
        }
        distances[i, ] <- cur.dist[indices[i, ]]
    }

    # Should produce the same nearest neighbors as GPvecchia
    expect_equal(pgp.nn$indices[-(1:5), ], gpv.nn[-(1:5), -1])
    # Should obtain the same nearest neighbors and distances when we calculate
    # the neighbors by brute force.
    expect_equal(pgp.nn$indices[-(1:5), ], indices[-(1:5), ])
    expect_equal(pgp.nn$distances[-(1:5), ], distances[-(1:5), ], tolerance = 1e-2)
})
