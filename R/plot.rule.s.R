#' @title Plot Stopping Rule
#' @description
#' Display a stopping rule graphically as a curve.
#'
#' @param x A rule object calculated by \code{calc.rule.s} function
#' @param ... Other parameters passed to the \code{plot} function.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' pocock.rule <- calc.rule.s(n = 30, tau = 100, p0 = 0.1, type = "Pocock", alpha = 0.05)
#' plot.rule.s(pocock.rule, col = "red")
#' }
plot.rule.s <- function(x,...){
  x <- x$Rule
  NextMethod('plot', type = "l", xlim = c(0, max(x[,1])),
             ylim = c(0, max(x[,2])+1), xlab = "Total follow up time",
             ylab = "# Events",... )
  # plot(x = rule[,1], y = rule[,2], type = "l", xlim = c(0, max(rule[,1])),
  #      ylim = c(0, max(rule[,2])+1), xlab = "Total follow up time",
  #      ylab = "# Events",...)
}
