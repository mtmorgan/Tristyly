#' Mating
#'
#' @name Mating
NULL

#' @describeIn Mating Create a matrix describing dissasortative mating
#'     and morph fertility. The i,jth entry in the matrix is the
#'     contribution of pollen from morph j to the pollen pool of morph
#'     i. \code{rowSums(M())} are morph-specific female
#'     fertilities. Currently only pure dissassortative mating is
#'     supported
#' @return A 3 x 3 numeric matrix of fertilities
#' @examples
#' M()
#' @export
M <- function() {
    m <- (1 - diag(3)) / 2
    dimnames(m) <- list(
        Female=levels(genetics$Morph),
        Male=levels(genetics$Morph))
    m
}

#' @describeIn Mating Summarize population morph frequencies.
#' @param population A numeric() vector of genotype counts, in the
#'     same form as returned by \code{isoplethy()}, but with
#'     \code{sum(population)} the size of the population.
#' @return named numeric(3) vector of morph frequencies
#' @export
morph_frequency <- function(population) {
    gtype <- as_genotype(population)
    .morph_frequency(gtype)
}

.morph_frequency <- function(gtype) {
    vapply(split(gtype, genetics$Morph), sum, numeric(1))
}

#' @describeIn Mating Mate individuals in \code{population} using the
#'     inheritance model defined by \code{G} and mating patterns
#'     summarized by \code{M}. Only pure disassortative mating and
#'     equal morph-specific female fertility are currently supported.
#' @param G See \code{G()}.
#' @param M See \code{M()}.
#' @return numeric() vector of sampled genotypes.
#' @examples
#' n <- isoplethy(1000)
#' morph_frequency(n)
#' n1 <- mate(n)
#' n1
#' morph_frequency(n1)
#' @importFrom stats rmultinom setNames
#' @export
mate <- function(population, G = Tristyly::G(), M = Tristyly::M()) {
    stopifnot(
        identical(dimnames(G), dimnames(Tristyly::G())),
        identical(dimnames(M), dimnames(Tristyly::M())))
    .mate_population(population, G, M)
}

.mate_population <- function(population, G, M) {
    gtype <- .mate(.as_genotype(population), G, M)
    .sample(gtype, sum(population))
}

.mate <- function(gtype, G, M) {
    ## FIXME: incorporate female fertility differences
    gamete_freq <- .gamete_frequency_by_morph(gtype, G)
    morph_freq <- rowSums(gamete_freq)

    exp <- numeric(length(cannonical_gtype))
    for (morph in names(morph_freq)[morph_freq != 0]) {
        female <- gamete_freq[morph,]
        female <- female / sum(female)

        male <- colSums(gamete_freq * M[morph,])
        male <- male / sum(male)
        if (!all(is.finite(male)))
            next

        offspring <- morph_freq[morph] * outer(female, male)
        if (is.null(exp))
            exp <- offspring
        else exp <- exp + offspring
    }

    ## collapse to cannonical form, ignoring sex of parent
    vapply(split(exp, cannonical_gtype), sum, numeric(1))
}
