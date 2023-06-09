% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.rule.s.R
\name{plot.rule.s}
\alias{plot.rule.s}
\title{Plot Stopping Rule (Survival Data)}
\usage{
\method{plot}{rule.s}(x, ...)
}
\arguments{
\item{x}{A 'rule.s' object calculated by \code{calc.rule.s()} function}

\item{...}{Other parameters passed to the \code{plot} function.}
}
\description{
Display a stopping rule graphically as a curve for time-to-event data
}
\examples{
pocock.rule <- calc.rule.s(n = 30, tau = 100, p0 = 0.1, type = "Pocock", alpha = 0.05)
plot(pocock.rule, col = "red")

}
