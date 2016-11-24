#' Inheritance
#'
#' @name Inheritance
NULL

#' @describeIn Inheritance Create a matrix summaring gametes produced
#'     by each genotype, allowing for arbitrary recombination between
#'     the S and M loci.
#' @param r numeric(1) recombination rate, between 0 (complete
#'     linkage) and 1 (free recombination).
#'
#' @examples
#' G(0)
#' G(1)
#' G(.1)
#' 
#' @export
G <- function(r = 0) {
    stopifnot(is.numeric(r), length(r) == 1L, !is.na(r), r >= 0, r <= 1)
    (1 - r) * G0 + r * G1
}

#' @describeIn Inheritance Create initial, approximately isoplethic,
#'     genotype frequencies in standard form.
#' @param n missing or integer(1). When missing, return genotype
#'     frequencies in populations at approximate isoplethic
#'     equilibirum. Otherwise, return a sample of size \code{n} drawn
#'     from a poualation at approximate isoplethy.
#' @return A vector of genotype frequencies (i.e., non-negative values
#'     summing to 1)
#' @examples
#' isoplethy()      # approximate isoplethy
#' isoplethy(100)  # sample from isoplethic population
#' @export
isoplethy <- function(n) {
    gtype <- setNames(
        c(0.333, 0.309, 0.024, 0.122, 0.122, 0.045, 0.045, 0, 0,
          0),
        c("sm/sm", "sM/sm", "sM/sM", "Sm/sm", "SM/sm", "Sm/sM",
          "SM/sM", "Sm/Sm", "SM/Sm", "SM/SM"))
    gtype <- as_genotype(gtype)
    if (!missing(n))
        gtype <- rmultinom(1, n, gtype)[,1]
    gtype
}

.stopifnot_is_single_integer <- function(n) {
    stopifnot(is.numeric(n), length(n) == 1L, !is.na(n))
}

.stopifnot_is_gtype <- function(gtype) {
    stopifnot(
        is.numeric(gtype), length(gtype) != 0, !anyNA(gtype),
        !is.null(names(gtype)), !anyDuplicated(names(gtype)),
        all(names(gtype) %in% genetics$Genotype),
        all(gtype >= 0),
        sum(gtype) > 0)
}

#' @describeIn Inheritance Create a vector of genotype frequencies in
#'     standard form.
#' @param gtype named numeric() vector of genotype frequencies. The
#'     vector must have at least one positive value; it cannot contain
#'     NAs. The names of the vector correspond to the \code{Genotype}
#'     column of the \code{genetics} data object; names cannot be
#'     duplicated. If missing, genotypes are set to approximate
#'     isoplethy under complete disassortative mating.
#' @export
as_genotype <- function(gtype) {
    .stopifnot_is_gtype(gtype)
    ## place frequencies in standard form
    result <- setNames(numeric(nrow(genetics)), genetics$Genotype)
    result[names(gtype)] <- gtype
    result / sum(result)
}

#' @describeIn Inheritance Calculate gamete freqencies from current
#'     genotype frequencies and mode of inheritance.
#' @param G numeric() matrix of genotype-to-gamete transition
#'     probabilities, from \code{G()}.
#' @return numeric() matrix of genotype x gamete expected frequencies
#' @export
gamete_frequency <- function(gtype, G) {
    gtype <- as_genotype(gtype)
    .gamete_frequency(gtype, G)
}

.gamete_frequency <- function(gtype, G)
    G * gtype

#' @describeIn Inheritance Calculate gamete frequencies produced by
#'     each morph from current genotype frequencies and mode of
#'     inheritance.
#' @return numeric() matrix of morph x gamete expected frequencies.
#' @export
gamete_frequency_by_morph <- function(gtype, G) {
    gtype <- as_genotype(gtype)
    .gamete_frequency_by_morph(gtype, G)
}

.gamete_frequency_by_morph <- function(gtype, G) {
    gametes <- .gamete_frequency(gtype, G)
    as.matrix(rowsum(gametes, genetics$Morph))
}    
