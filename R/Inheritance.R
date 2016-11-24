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

#' @describeIn Inheritance Create initial genotype frequencies in
#'     standard form.
#' @param gtype named numeric() vector of genotype frequencies. The
#'     vector must have at least one positive value; it cannot contain
#'     NAs. The names of the vector correspond to the \code{Genotype}
#'     column of the \code{genetics} data object; names cannot be
#'     duplicated. If missing, genotypes are set to approximate
#'     isoplethy under complete disassortative mating.
#' @return A vector of genotype frequencies (i.e., non-negative values
#'     summing to 1)
#' @examples
#' gtype_init()    # approximate isoplethy
#' gtype_init(sm/sm=100, sM/sm=70, sM/sM=20, Sm/sm=50)
#' @export
gtype_init <- function(gtype) {
    if (missing(gtype)) {
        ## Approximate isoplethy
        gtype <- setNames(
            c(0.333, 0.309, 0.024, 0.122, 0.122, 0.045, 0.045, 0, 0,
              0),
            c("sm/sm", "sM/sm", "sM/sM", "Sm/sm", "SM/sm", "Sm/sM",
              "SM/sM", "Sm/Sm", "SM/Sm", "SM/SM"))
    }
    stopifnot(
        is.numeric(gtype), length(gtype) != 0, !anyNA(gtype),
        !is.null(names(gtype)), !anyDuplicated(names(gtype)),
        all(names(gtype) %in% genetics$Genotype),
        all(gtype >= 0),
        sum(gtype) > 0)
    .gtype_frequency(gtype)
}

.gtype_frequency <- function(gtype) {
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
    gtype <- gtype_init(gtype)
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
    gtype <- gtype_init(gtype)
    .gamete_frequency_by_morph(gtype, G)
}

.gamete_frequency_by_morph <- function(gtype, G) {
    gametes <- .gamete_frequency(gtype, G)
    as.matrix(rowsum(gametes, genetics$Morph))
}    
