#' Mating
#'
#' @name Mating
NULL

.cannonical_gtype <- local({
    ## map between all gentoypes and cannonical (independent of
    ## parent-of-origin) form.
    ##
    ## use of local() and function to allow reference to lazy object
    ## G0, which is not available at package load time.
    value <- NULL
    function() {
        if (is.null(value)) {
            ## first time through
            gametes <- colnames(G(0))
            all <- outer(gametes, gametes, paste, sep="/")
            cannonical <- outer(gametes, gametes, function(x, y) {
                paste(pmax(x, y), pmin(x, y), sep="/")
            })
            value <<- data.frame(
                All=as.vector(all),
                Cannonical=factor(
                    as.vector(cannonical),
                    levels=genetics$Genotype),
                stringsAsFactors=FALSE)
        }
        value
    }
})

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
    .stopifnot_is_gtype(population)
    stopifnot(
        identical(dimnames(G), dimnames(G0)),
        identical(dimnames(M), dimnames(Tristyly::M())))
    gtype <- .mate(as_genotype(population), G, M)
    .sample(gtype, sum(population))
}

.mate <- function(gtype, G, M) {
    ## FIXME: incorporate female fertility differences
    gamete_freq <- .gamete_frequency_by_morph(gtype, G)
    morph_freq <- rowSums(gamete_freq)

    exp <- NULL
    for (morph in names(morph_freq)[morph_freq != 0]) {
        female <- gamete_freq[morph,]
        female <- female / sum(female)

        male_freq <- morph_freq * M[morph,]
        male <- colSums(gamete_freq * male_freq)
        male <- male / sum(male)
        if (!all(is.finite(male)))
            next

        offspring <- morph_freq[morph] * outer(female, male)
        if (is.null(exp))
            exp <- offspring
        else exp <- exp + offspring
    }

    ## collapse to cannonical form, ignoring sex of parent
    vapply(split(exp, .cannonical_gtype()$Cannonical), sum, numeric(1))
}
