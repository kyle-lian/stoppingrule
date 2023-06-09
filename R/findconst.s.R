#' @title Search for Calibration Value (Survival Data)
#' @description
#' Internal workhorse function to calculate the calibration constant value that attains level alpha for given method for time-to-event data
#'
#' @param n Maximum sample size for safety monitoring
#' @param tau Observation period
#' @param p0 The toxicity rate under the null hypothesis
#' @param type The method used for constructing the stopping rule
#' @param alpha The desired type I error/false positive rate for the stopping rule
#' @param param Extra parameter(s) needed for certain stopping rule methods. For Wang-Tsiatis test, this is the Delta parameter. For modified SPRT, this is the targeted alternative toxicity rate p1. For Bayeian Gamma-Poisson model, this is the pair of hyperparameters for the gamma prior on the toxicity rate.
#'
#' @return The calibration constant used for subsequent stopping boundary calculation
#'
#'
#'
findconst.s <- function(n, tau, p0, type, alpha, param = NULL){

  inner <- function(n, tau, p0, type, cval, param){
    bnd = calc.bnd.s(n = n, tau = tau, p0 = p0, type = type, cval = cval, param = param)
    return(stopping.prob(bnd, p = p0))
  }

  if (type == 'Pocock'|type == 'OBF'|type == 'WT'){
    l = qnorm(1-alpha)
    u = qnorm(1 - alpha/(n*p0))
  } else if (type == "Bayesian"){
    lambda0 <- -log(1 - p0)/tau
    Umax <- n*tau

    if (length(param) == 1){
      k <- param
      s <- k/lambda0
    } else {
      k <- param[1] # shape
      s <- param[2] # rate
    }
    l = (1 - pgamma(lambda0, shape = (k + 1), rate = (s + Umax))) + 0.01
    u = 0.999
  } else {
    l = exp(0.5*qchisq(1-alpha,1))
    u = exp(0.5*qchisq(1-alpha/(n*p0),1))
  }

  cval <- uniroot(function(cval) {
    inner(n = n, tau = tau, p0 = p0, type = type, param = param, cval = cval)$Stop.prob - alpha
  }, lower = l, upper = u, extendInt = "yes")$root

  return(cval)
}
