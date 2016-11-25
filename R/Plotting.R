#' Plotting
#'
#' @name Plotting
NULL

#' @describeIn Plotting Plot simulated changes in morph frequency.
#' @param iterations An object returned by \code{\link{iterate}}.
#' @return 
#' @examples
#' plot_iterations(iterate(isoplethy(100), 1000))
#' @export
plot_iterations <- function(iterations, ..., type="b", add=FALSE) {
    long <- iterations[iterations$Morph == "Long", "Frequency"][[1]]
    mid <- iterations[iterations$Morph == "Mid", "Frequency"][[1]]

    max <- sqrt(1 - .5^2)
    y <- max * long
    x_offset <- y / tan(pi/3)
    x <- x_offset + max * mid / sin(pi/3)

    plot(y ~ x, ..., type=type, xlim=c(0, 1), ylim=c(0, max), axes=FALSE)
    lines(c(0, 1, .5, 0), c(0, 0, max, 0))
}
