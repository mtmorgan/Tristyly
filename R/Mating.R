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
                xstringsAsFactors=FALSE)
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
#'     same form as returned by \code{gtype_init()}, but with
#'     \code{sum(population)} the size of the population.
#' @return named numeric(3) vector of morph frequencies
#' @export
morph_frequency <- function(population) {
    gtype <- gtype_init(population)
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
#' gtype <- gtype_init()
#' n <- setNames(as.integer(gtype * 1000), names(gtype))
#' morph_frequency(n)
#' n1 <- mate(n)
#' n1
#' morph_frequency(n1)
#' @export
mate <- function(population, G = Tristyly::G(), M = Tristyly::M()) {
    stopifnot(
        ## population is checked in call to .mate
        identical(dimnames(G), dimnames(G0)),
        identical(dimnames(M), dimnames(Tristyly::M())))
    expected_frequency <- .mate(population, G, M)
    obs <- rmultinom(1L, sum(population), expected_frequency)
    setNames(as.vector(obs), rownames(obs))
}

.mate <- function(gtype, G, M) {
    ## FIXME: incorporate female fertility differences
    gamete_freq <- gamete_frequency_by_morph(gtype, G)
    morph_freq <- rowSums(gamete_freq)

    result <- NULL
    for (morph in names(morph_freq)) {
        female <- gamete_freq[morph,]
        n <- sum(female)
        if (!n)
            next
        female <- female / n

        male <- colSums(gamete_freq * M[morph,])
        n <- sum(male)
        if (!n)
            next
        male <- male / n

        offspring <- morph_freq[morph] * outer(female, male)
        if (is.null(result))
            result <- offspring
        else result <- result + offspring
    }

    ## collapse to cannonical form, ignoring sex of parent
    vapply(split(result, .cannonical_gtype()$Cannonical), sum, numeric(1))
}
