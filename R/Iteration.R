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
#' population <- rmultinom(1, 10, gtype_init())[,1]
#' iterate(population, 1000)
#' @importFrom tibble tibble
#' @export
iterate <-
    function(population, generations, G = Tristyly::G(), M = Tristyly::M())
{
    stopifnot(
        is.numeric(generations), length(generations) == 1L, !is.na(generations),
        generations > 0L)

    result <- numeric(3L * generations)
    result[1:3] <- morph_frequency(population)
    for (i in seq_len(generations - 1L)) {
        population <- mate(population, G, M)
        result[i * 3L + 1:3] <- morph_frequency(population)
        if (sum(population != 0) == 1L) {
            idx <- seq(i * 3L + 1L, generations * 3L)
            result[idx] <- morph_frequency(population)
            break
        }
    }

    tibble(
        Generation=rep(seq_len(generations), each=3),
        Morph=factor(
            rep(names(morph_frequency(population)), generations),
            levels=levels(genetics$Morph)),
        Frequecy=result)
}       

#' @describeIn Iteration Repeatedly iterate a population until
#'     monomorphism.
#' @param times numeric(1) times to iterate the population to
#'     monomorphism.
#' @param progress.interval numeric(1) report progress every \code{verbose}
#'     iterations. Suppress output with \code{verbose = Inf}.
#' @examples
#' population_size <- 10
#' population <- rmultinom(1, population_size, gtype_init())[,1]
#' morph_frequency(population)
#' mono <- iterate_to_monomorphism(population, 30)
#' table(mono$Genotype)
#' plot(ecdf(mono$Generation), xlab="Generation")
#' @export
iterate_to_monomorphism <-
    function(population, times, G = Tristyly::G(), M = Tristyly::M(),
             progress.interval=10)
{
    result <- .iterate_to_morphism(
        population, times, 1L, G, M, progress.interval)

    result$Morph <-
        vapply(result$Morph, function(elt) names(which.max(elt)), character(1))
    tibble(
        Generation=result$Generation,
        Morph_kept=factor(result$Morph, levels=levels(genetics$Morph))
    )
}

#' @describeIn Iteration Repeatedly iterate a population until loss of
#'     one morph.
#' @export
iterate_to_dimorphism <-
    function(population, times, G = Tristyly::G(), M = Tristyly::M(),
             progress.interval=50)
{
    result <- .iterate_to_morphism(
        population, times, 2L, G, M, progress.interval)

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

.iterate_to_morphism <-
    function(population, times, n_morphs, G, M, progress.interval=10)
{
    stopifnot(
        is.numeric(times), length(times) == 1L, !is.na(times), times > 0L,
        is.numeric(progress.interval), length(progress.interval) == 1L,
        !is.na(progress.interval), progress.interval > 0L)

    generation <- integer(times)
    morph <- vector("list", times)
    for (i in seq_len(times)) {
        if ((i %% progress.interval) == 0L)
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
        if (sum(morph_frequency(population) != 0) <= n_morphs)
            break
        population <- mate(population, G, M)
    }
    list(Iteration=iteration, Population=population)
}
