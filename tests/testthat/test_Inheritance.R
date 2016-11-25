test_that(".gamete_frequency()  works", {
    gamete_frequency <- function(g, G)
        Tristyly:::.gamete_frequency(as_genotype(g), G)

    obs <- colSums(gamete_frequency(c(`sm/sm`=1), G(0)))
    exp <- structure(c(1, 0, 0, 0), .Names = c("sm", "sM", "Sm", "SM"))
    expect_identical(obs, exp)
    
    obs <- colSums(gamete_frequency(c(`sM/sm`=1), G(0)))
    exp <- structure(c(0.5, 0.5, 0, 0), .Names = c("sm", "sM", "Sm", "SM"))
    expect_identical(obs, exp)

    obs <- colSums(gamete_frequency(c(`sM/sM`=1), G(0)))
    exp <- structure(c(0, 1, 0, 0), .Names = c("sm", "sM", "Sm", "SM"))
    expect_identical(obs, exp)

    obs <- colSums(gamete_frequency(c(`Sm/sm`=1), G(0)))
    exp <- structure(c(0.5, 0, 0.5, 0), .Names = c("sm", "sM", "Sm", "SM"))
    expect_identical(obs, exp)

    obs <- colSums(gamete_frequency(c(`SM/sm`=1), G(0)))
    exp <- structure(c(0.5, 0, 0, 0.5), .Names = c("sm", "sM", "Sm", "SM"))
    expect_identical(obs, exp)

    obs <- colSums(gamete_frequency(c(`SM/sm`=1), G(1)))
    exp <- structure(rep(1, 4) / 4, .Names = c("sm", "sM", "Sm", "SM"))
    expect_identical(obs, exp)

    obs <- colSums(gamete_frequency(c(`Sm/sM`=1), G(0)))
    exp <- structure(c(0, 0.5, 0.5, 0), .Names = c("sm", "sM", "Sm", "SM"))
    expect_identical(obs, exp)

    obs <- colSums(gamete_frequency(c(`Sm/sM`=1), G(1)))
    exp <- structure(rep(1, 4) / 4, .Names = c("sm", "sM", "Sm", "SM"))
    expect_identical(obs, exp)

    obs <- colSums(gamete_frequency(c(`SM/sM`=1), G(0)))
    exp <- structure(c(0, 0.5, 0, 0.5), .Names = c("sm", "sM", "Sm", "SM")) 
    expect_identical(obs, exp)
    
    obs <- colSums(gamete_frequency(c(`Sm/Sm`=1), G(0)))
    exp <- structure(c(0, 0, 1, 0), .Names = c("sm", "sM", "Sm", "SM")) 
    expect_identical(obs, exp)

    obs <- colSums(gamete_frequency(c(`SM/Sm`=1), G(0)))
    exp <- structure(c(0, 0, 0.5, 0.5), .Names = c("sm", "sM", "Sm", "SM")) 
    expect_identical(obs, exp)
    
    obs <- colSums(gamete_frequency(c(`SM/SM`=1), G(0)))
    exp <- structure(c(0, 0, 0, 1), .Names = c("sm", "sM", "Sm", "SM")) 
    expect_identical(obs, exp)
})
