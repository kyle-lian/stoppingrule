#' @title Stopping Boundary Calculation (Survival Data)
#' @description
#' Internal workhorse function to calculate stopping boundary for a given method
#' @param n Maximum sample size for safety monitoring
#' @param tau Observation period
#' @param p0 The toxicity rate under the null hypothesis
#' @param cval Critical value for stopping rule method
#' @param type The method used for constructing the stopping rule
#' @param param Extra parameter(s) needed for certain stopping rule methods. For Wang-Tsiatis test, this is the Delta parameter. For modified SPRT, this is the targeted alternative toxicity rate p1. For Bayeian Gamma-Poisson model, this is the pair of hyperparameters for the gamma prior on the toxicity rate.
#'
#' @import pracma
#' @return A list of three: tau, number of events that can trigger a stop, and the corresponding total follow up time.
#' @export
#' @examples
#' \dontrun{
#' calc.bnd.s(n = 30, tau = 100, p0 = 0.1, cval = 4.7, type = "Pocock")
#' calc.bnd.s(n = 30, tau = 100, p0 = 0.1, cval = 3, type = "SPRT", param = 0.3)
#' }
#

calc.bnd.s <- function(n, tau, p0, cval, type, param = NULL){
  # library(pracma)
  lambda0 <- -log(1 - p0)/tau
  W <- cval
  Umax <- n*tau

  if (type == "Pocock"){
    f <- function(U){
      lambda0*U + W*sqrt(lambda0*U)
    }
    f.inverse <- function(D){
      (D + 0.5*(W^2 - sqrt(4*D*W^2 + W^4)))/lambda0
    }
  }

  else if (type == "OBF"){
    f <- function(U){
      lambda0*U + W*sqrt(lambda0*Umax)
    }
    f.inverse <- function(D){
      (D - W*sqrt(lambda0*Umax))/lambda0
    }
  }

  else if (type == "WT"){
    Delta <- param
    f <- function(U){
      lambda0*U + W*sqrt(lambda0)*Umax^(0.5 - Delta)*U^Delta
    }
    f.inverse <- function(d, lower = 0, upper = 2*n*tau){
      uniroot(function(u){f(u) - d}, lower = lower, upper = upper)$root
    }
  }

  else if (type == "SPRT"){
    p1 <- param
    lambda1 <- -log(1 - p1)/tau
    f <- function(U){
      (log(W) + (lambda1 - lambda0)*U)/log(lambda1/lambda0)
    }
    f.inverse <- function(D){
      (D*(log(lambda1) - log(lambda0)) - log(W))/(lambda1 - lambda0)
    }
  }

  else if (type == "MaxSPRT"){
    f <- function(U){
      if (U ==0 ){return(0)}
      else {return((log(W) - lambda0*U)*(lambertWp(1/exp(1)*(log(W)/(lambda0*U) - 1)))^(-1))}
    }
    f.inverse <- function(D){
      -D/lambda0*lambertWp(-exp(-log(W)/D - 1))
    }
  }

  if (type != "Bayesian"){
    S <- seq(from = max(ceiling(f(0)),1), to = ceiling(f(Umax)))
    # FirstPositive <- which(S > 0)[1]
    # S <- S[FirstPositive:length(S)]
    dmin <- S[1]
    dmax <- S[length(S)]
    m <- dmax - dmin + 1

    ud <- rep(NA, m)
    for (i in 1:m){
      ud[i] <- f.inverse(S[i])
    }
    ud[length(ud)] <- Umax
  }


  if (type == "Bayesian"){
    # hyperparameters
    if (length(param) == 1){
      k <- param
      l <- k/lambda0
    } else {
      k <- param[1] # shape
      l <- param[2] # rate
    }

    # find dmin and dmax
    post <- NULL
    for (d in 1:n){
      post[d] <- 1 - pgamma(lambda0, shape = (k + d), rate = l)
    }
    dmin = min(which(post >= W))

    post1 <- NULL
    for (d in 1:n){
      post1[d] <- 1 - pgamma(lambda0, shape = (k + d), rate = (l + Umax))
    }
    dmax = min(which(post1 >= W))
    S = seq(dmin, dmax)
    ud <- rep(NA, length.out = length(S))
    for (j in 1:(length(S) - 1)){
      inner <- function(ud){
        1 - pgamma(lambda0, shape = (k + S[j]), rate = (l + ud))
      }
      if( j == 1){
        ud[j] <- uniroot(function(ud){inner(ud) - W}, lower = 0, upper = Umax,
                         extendInt = "no", trace = 2)$root
      } else {
        ud[j] <- uniroot(function(ud){inner(ud) - W}, lower = ud[j-1], upper = Umax,
                         extendInt = "no", trace = 2)$root
      }
    }
    ud[length(ud)] <- Umax
  }

  return(list(tau = tau, S = S, ud = ud))
}
