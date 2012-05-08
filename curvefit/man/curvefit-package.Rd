\name{curvefit-package}
\alias{curvefit-package}
\alias{curvefit}
\docType{package}
\title{Bayesian curve fitting}
\description{
  Fits data using Bayesian techniques.


}
\details{
\tabular{ll}{
Package: \tab curvefit\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-04-22\cr
License: \tab GPLv3\cr
LazyLoad: \tab yes\cr
}

The function curvefit fits a piecewise polynomial to data. The
function `curvefit.at' is used to evaluate the fitted function at
specified values.

}
\usage{
curvefit(formula, data, ...)

\method{curvefit}{default}(formula, data, prior_mean,
max_knots=length(data), c_param=.4, poly_deg=2, knot_continuity=1,
mse_relative=10^-3, mse_absolute=10^3, burnin=10, diagnostics=FALSE, ...)
\method{curvefit}{print}(fit, ...)
\method{curvefit}{at}(fit, x, ...)
}
\arguments{
  \item{formula}{ a formula specifying the predictor and response
	variable. }
  
  \item{data}{ the data to be fit. If the formula specified is "y ~ x",
	then "data" must contain columns "y" and "x". }

  \item{prior_mean}{ the prior mean of the number of knots in the fitted
	model. Must be a nonnegative integer. }

  \item{max_knots}{ the maximum number of knot locations in the
	model. Increase to increase the granularity of fits. }

  \item{c_param}{ a parameter that controls the ratio between 'birth'
	and 'death' steps and 'move' steps. Must be in (0, .5). }

  \item{poly_deg}{ the degree of fitted polynomials. Must be a
	positive integer. }

  \item{knot_continuity}{ controls the smoothness of the resulting
	fit. If zero, the result is possibly noncontinuous. If positive, the
	fit will have "knot_continuity - 1" continuous derivatives. }

  \item{mse_relative}{ the relative mean-squared error required to
	halt. Decrease to improve fit, at the expense of running time. }

  \item{mse_absolute}{ the absolute mean-squared error required to
	halt. Decrease to improve fit, at the expense of running time. }

  \item{burnin}{ minimum number of iterations we require before testing
  for absolute and relative halting conditions. }

  \item{diagnostics}{ if TRUE, prints out diagnostics as the fit
	runs. Can be useful to see how fast a fit is converging. }

  \item{fit}{ a curvefit object. }

  \item{x}{ a vector of values to evaluate the fit at. }

  \item{\dots}{ ignored. }
}
\value{
  An object of class \code{curvefit}, a list including elements
  \item{coef}{ a named vector of fit coefficients. }
  \item{niters}{ the number of iterations required to terminate. }
  \item{mse}{ the mean-squared error of the fitted model. }
  \item{knots}{ a vector of interior knot indices in the fitted
	model. The value of knot points can be found by accessing the
	corresponding index of the "locations" parameter. }
  \item{locations}{ a vector of candidate knot locations. }
  \item{l}{ the degree of polynomials used. }
  \item{l0}{ the degree of continuity of the fit. If 0, the fit is not
	guaranteed to be continuous. If positive, the fit has l0 - 1
	continuous derivatives. }
}
\author{
Patrick Steele

Maintainer: Patrick Steele, <prs233@cornell.edu>

}
\references{
"Automatic Bayesian curve fitting", Denison et al., 1998
}

\keyword{ regression }

\examples{
# data(EuStockMarkets

# x = 1:200
# y = EuStockMarkets[x, 1]

# fit = curvefit(y ~ x, list(y=y, x=x), prior_mean=10)

# plot(x, fit.at(fit, x))
}
