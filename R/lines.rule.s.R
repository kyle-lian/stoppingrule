#' @title Add Stopping Rule Curve to Current Plot (Survival Data)
#'
#' @param x A rule object calculated by \code{calc.rule.s} function
#' @param ... Other options to be passed to generic \code{lines} function
#'
#' @return No return value, function solely modifies current plot
#' @export
#'
#' @examples
#' \dontrun{
#' pocock.rule <- calc.rule.s(n = 30, tau = 100, p0 = 0.1, type = "Pocock", alpha = 0.05)
#' OBF.rule <- calc.rule.s(n = 30, tau = 100, p0 = 0.1, type = "OBF", alpha = 0.05)
#' }
#'
lines.rule.s = function(x,...) {
  x = x$Rule
  NextMethod("lines",type='l',...)
}
