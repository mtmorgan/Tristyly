#' Iteration
#'
#' @name Iteration
NULL

#' @describeIn Iteration Iterate a single population for a fixed
#'     number of generations.
#' @param population See \code{\link{Mating}}.
#' @param generations the number of generations to iterate; the initial
#'     population is generation 1.
#' @param G See \code{\link{G}()}.
#' @param M See \code{\link{M}()}.
#' @return a tibble (data.frame) summarizing morph frequencies in each
#'     generation, including the first (generation 1)
#' @examples
#' iterate(isoplethy(10), 1000)
#' @importFrom tibble tibble
#' @export
iterate <-
    function(population, generations, G = Tristyly::G(), M = Tristyly::M())
{
    .stopifnot_is_gtype(population)
    stopifnot(
        is.numeric(generations), length(generations) == 1L, !is.na(generations),
        generations > 0L)

    result <- numeric(3L * generations)
    result[1:3] <- morph_frequency(population)
    for (i in seq_len(generations - 1L)) {
        population <- .mate_population(population, G, M)

        result[i * 3L + 1:3] <- morph_frequency(population)
        if (sum(population != 0) <= 1L) {
            idx <- seq(i * 3L + 1L, generations * 3L)
            result[idx] <- morph_frequency(population)
            break
        }
    }

    tibble(
        Generation=rep(seq_len(generations), each=3),
        Morph=factor(
            rep(levels(genetics$Morph), generations),
            levels=levels(genetics$Morph)),
        Frequency=result)
}       

#' @describeIn Iteration Iterate independent populations to loss of
#'     one morph
#' @param N integer(1) Population size
#' @param verbose missing or numeric(1) When present, report progress
#'     every \code{verbose} iterations.
#' @return A \code{tibble} (data.frame) with columns \code{Generation}
#'     (generation of loss of first morph) and \code{Morph_lost}
#'     (morph lost). If two morphs are lost in the same generation,
#'     the value of \code{Morph} is \code{NA}.
#' @examples
#' di <- replicate_to_dimorphism(10, 100)
#' table(di$Morph_lost)
#' plot(ecdf(di$Generation), xlab="Generation")
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
replicate_to_dimorphism <-
    function(N, times, G = Tristyly::G(), M = Tristyly::M(), verbose=FALSE)
{
    stopifnot(
        is.numeric(N), length(N) == 1, !is.na(N), N > 0,
        is.numeric(times), length(times) == 1, !is.na(times), times > 0)
    if (verbose)
        pb <- txtProgressBar(style=3)

    result <- vector("list", times)
    for (i in seq_len(times)) {
        if (verbose && (i %% verbose) == 0)
            setTxtProgressBar(pb, i / times)
        result[[i]] <- .iterate_to_morphism(isoplethy(N), 1, 2L, G, M, FALSE)
    }

    if (verbose)
        close(pb)

    generation <- vapply(result, "[[", integer(1), 1)
    morph <- vapply(result, function(elt) {
        elt <- elt$Morph[[1]]
        if (sum(elt == 0) > 1L) {
            NA_character_
        } else
            names(which.min(elt))
    }, character(1))
    tibble(
        Generation=generation,
        Morph_lost=factor(morph, levels=levels(genetics$Morph))
    )
}

#' @describeIn Iteration Repeatedly iterate a population until loss of
#'     one morph.
#' @param times numeric(1) times to iterate the population to
#'     monomorphism.
#' @return A \code{tibble} (data.frame) with columns \code{Generation}
#'     (generation of loss of first morph) and \code{Morph_lost}
#'     (morph lost). If two morphs are lost in the same generation,
#'     the value of \code{Morph} is \code{NA}.
#' @examples
#' population <- isoplethy(20)
#' morph_frequency(population)
#' di <- iterate_to_dimorphism(population, 100)
#' table(di$Morph_lost)
#' plot(ecdf(di$Generation), xlab="Generation")
#' @export
iterate_to_dimorphism <-
    function(population, times, G = Tristyly::G(), M = Tristyly::M(),
             verbose=FALSE)
{
    result <- .iterate_to_morphism(population, times, 2L, G, M, verbose)

    result$Morph <- vapply(result$Morph, function(elt) {
        if (sum(elt == 0) > 1L) {
            NA_character_
        } else
            names(which.min(elt))
    }, character(1))
    tibble(
        Generation=result$Generation,
        Morph_lost=factor(result$Morph, levels=levels(genetics$Morph))
    )
}

#' @describeIn Iteration Repeatedly iterate a population until
#'     monomorphism.
#' @return A \code{tibble} (data.frame) with columns \code{Generation}
#'     (generation until loss of two morphs) and \code{Morph_kept}
#'     (the morph remaining in the population).
#' @examples
#' population <- isoplethy(10)
#' morph_frequency(population)
#' mono <- iterate_to_monomorphism(population, 30)
#' table(mono$Morph_kept)
#' plot(ecdf(mono$Generation), xlab="Generation")
#' @export
iterate_to_monomorphism <-
    function(population, times, G = Tristyly::G(), M = Tristyly::M(),
             verbose=FALSE)
{
    result <- .iterate_to_morphism(population, times, 1L, G, M, verbose)

    result$Morph <- vapply(result$Morph, function(elt) {
        names(which.max(elt))
    }, character(1))
    tibble(
        Generation=result$Generation,
        Morph_kept=factor(result$Morph, levels=levels(genetics$Morph))
    )
}

.iterate_to_morphism <-
    function(population, times, n_morphs, G, M, verbose)
{
    stopifnot(
        is.numeric(times), length(times) == 1L, !is.na(times), times > 0L,
        is.logical(verbose) || is.numeric(verbose),
        length(verbose) == 1L, !is.na(verbose))

    generation <- integer(times)
    morph <- vector("list", times)
    for (i in seq_len(times)) {
        if (verbose && (i %% verbose) == 0L)
            message(i)

        result <- .iterate1_to_morphism(population, n_morphs, G, M)

        generation[i] <- result$Iteration
        morph[[i]] <- morph_frequency(result$Population)
    }

    list(Generation=generation, Morph=morph)
}

.iterate1_to_morphism <- function(population, n_morphs, G, M) {
    iteration <- 0L
    repeat {
        iteration <- iteration + 1L
        if (sum(.morph_frequency(population) != 0) <= n_morphs)
            break
        population <- .mate_population(population, G, M)
    }
    list(Iteration=iteration, Population=population)
}
