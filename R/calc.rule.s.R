#' @title Stopping Rule Calculation (Survival Data)
#' @description
#'  Calculate a stopping rule for safety monitoring
#'
#' @param n Maximum sample size for safety monitoring
#' @param tau Observation period
#' @param p0 The toxicity rate under the null hypothesis
#' @param type The method used for constructing the stopping rule
#' @param param Extra parameter(s) needed for certain stopping rule methods. For Wang-Tsiatis test, this is the Delta parameter. For modified SPRT, this is the targeted alternative toxicity rate p1. For Bayeian Gamma-Poisson model, this is the pair of hyperparameters for the gamma prior on the toxicity rate.
#' @param alpha The desired type I error/false positive rate for the stopping rule
#'
#' @return A list of four: 1. A matrix with two columns: total follow up time and their corresponding rejection boundary. 2. value tau to be stored for later use. 3. The calibration constant value used for calculation. 4. Stopping probability at each stage.
#' @export
#'
#' @examples
#' \dontrun{
#' calc.rule.s(n = 30, tau = 100, p0 = 0.1, type = "Pocock", alpha = 0.05)
#' }
calc.rule.s <- function(n, tau, p0, type, param=NULL, alpha){
  cval <- findconst.s(n = n, tau = tau, p0 = p0, type = type, param = param, alpha = alpha)

  bdry <- calc.bnd.s(n = n, tau = tau, p0 = p0, cval = cval, type = type, param = param)
  stage.stop.prob <- stopping.prob(bnd = bdry, p = p0)$stage.stop.prob
  val <- cbind(floor(bdry$ud), bdry$S)
  colnames(val) <- c("Total follow up time","Reject bdry")
  val2 <- list(Rule = val, tau = tau, cval = cval, stage.stop.prob = stage.stop.prob)
  class(val2) <- "rule.s"
  return(val2)
}
