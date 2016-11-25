#' Plotting
#'
#' @name Plotting
NULL

#' @describeIn Plotting Plot simulated changes in morph frequency.
#' @param iterations An object returned by \code{\link{iterate}}.
#' @param ... Additional graphical parameters used by \code{lines()}
#'     to plot population trajectories.
#' @param type character(1) Type of plot to create; default plots
#'     trajectories as points joined by lines.
#' @param add logical(1) Create a new plot (\code{add = FALSE},
#'     default) or add a trajectory to an existing plot.
#' @return None; invoked for side effect.
#' @examples
#' plot_iterations(iterate(isoplethy(50), 1000), pch=20)
#' @importFrom graphics lines par plot text
#' @export
plot_iterations <- function(iterations, ..., type="b", add=FALSE) {
    long <- iterations[iterations$Morph == "Long", "Frequency"][[1]]
    mid <- iterations[iterations$Morph == "Mid", "Frequency"][[1]]

    max <- sqrt(1 - .5^2)
    y <- max * long
    x <- y / tan(pi/3) + max * mid / sin(pi/3)

    if (!add)
        .plot_crosby(max)
    lines(y ~ x, ..., type="b")
}

.plot_crosby <- function(max) {
    opar <- par(mar=c(4, 1, 1, 1))
    on.exit(par(opar))
    plot(integer(), type="n", xlim=c(0, 1), ylim=c(0, max), asp=1,
         ann=FALSE, axes=FALSE)
    lines(c(0, 1, .5, 0), c(0, 0, max, 0))
    text(c(0.5, 0.25, 0.75), c(0, max / 2, max / 2), c("Long", "Mid", "Short"),
         pos=c(1, 2, 4), offset=.8, cex=1.8)
}
