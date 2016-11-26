#' Plotting
#'
#' @name Plotting
NULL

#' @describeIn Plotting Plot changes in morph frequency across
#'     generations. The starting population is represented by a solid
#'     dot.
#' @param iterations An object returned by
#'     \code{\link{iterate_stochastic}} or
#'     \code{\link{iterate_deterministic}}.
#' @param ... Additional graphical parameters used by \code{lines()}
#'     to plot population trajectories.
#' @param pch Plotting symbol(s); see \code{points}. If missing, the
#'     first point is a solid dot and subsequent points open circles.
#' @param type character(1) Type of plot to create; default plots
#'     trajectories as points joined by lines.
#' @param add logical(1) Create a new plot (\code{add = FALSE},
#'     default) or add a trajectory to an existing plot.
#' @return None; invoked for side effect.
#' @examples
#' plot_iterations(iterate_stochastic(isoplethy(50), 1000), pch=20)
#' @importFrom graphics lines par plot text
#' @export
plot_iterations <- function(iterations, ..., pch, type="b", add=FALSE) {
    short <- iterations[iterations$Morph == "Short", "Frequency"][[1]]
    mid <- iterations[iterations$Morph == "Mid", "Frequency"][[1]]

    max <- sqrt(1 - .5^2)
    y <- max * short
    x <- y / tan(pi/3) + max * mid / sin(pi/3)

    if (!add)
        .plot_crosby(max)
    if (missing(pch))
        pch <- c(16L, rep(1L, length(x) - 1L))

    lines(y ~ x, ..., pch=pch, type="b")
}

.plot_crosby <- function(max) {
    opar <- par(mar=c(4, 1, 1, 1))
    on.exit(par(opar))
    plot(integer(), type="n", xlim=c(0, 1), ylim=c(0, max), asp=1,
         ann=FALSE, axes=FALSE)
    lines(c(0, 1, .5, 0), c(0, 0, max, 0))
    text(c(0.5, 0.25, 0.75), c(0, max / 2, max / 2), c("Short", "Mid", "Long"),
         pos=c(1, 2, 4), offset=.8, cex=1.8)
}
