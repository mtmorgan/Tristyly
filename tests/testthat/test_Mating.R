test_that("mate() works", {
    .mate <- Tristyly:::.mate
    G <- G(0)
    M <- M()

    g <- g0 <- isoplethy()
    for (i in 1:2000)
        g <- .mate(g, G, M)
    expect_lt(sum(abs(g0 - g)), 1e-4)
})
