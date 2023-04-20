#' @title Plot Stopping Rule (Survival Data)
#' @description
#' Display a stopping rule graphically as a curve for time-to-event data
#'
#' @param x A 'rule.s' object calculated by \code{calc.rule.s()} function
#' @param ... Other parameters passed to the \code{plot} function.
#'
#' @return
#' @export
#'
#' @examples
#' pocock.rule <- calc.rule.s(n = 30, tau = 100, p0 = 0.1, type = "Pocock", alpha = 0.05)
#' plot(pocock.rule, col = "red")
#'
plot.rule.s <- function(x,...){
  x = x$Rule
  NextMethod('plot', type = "l", xlim = c(0, max(x[,1])),
             ylim = c(0, max(x[,2])+1), xlab = "Total follow up time",
             ylab = "# Events",... )
}
