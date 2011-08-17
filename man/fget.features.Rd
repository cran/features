\name{fget}
\alias{fget}
\alias{fget.features}
\title{
Utility function to extract features from an object of \code{class "features"}.
}
\description{
This is a utility function that operates only on objects of \code{class "features"} that are returned by the function \code{features()}.  Some of the key extracted features include mean value, first and second derivatives, critical points 
(i.e. local maxima and minima), curvature of function at critical points, wiggliness of the function, noise in data, and outliers in data.}
\usage{
fget(x)
}

\arguments{
\item{x}{An object of \code{class "features" } } 
}


\value{
A list with a number of extracted features of the underlying smooth function:
\item{f}{: a numeric vector containing 10 basic features of the function (in order): mean, min, max, std.dev, noise, signal-to-noise ratio, minimum and maximum of first derivative, wiggliness, and number of critical points}

\item{cpts}{: Locations of the critical points (e.g. points where the first derivative is zero)}

\item{curvature}{: Value of second derivative at the critical points}

\item{outliers}{: Locations of outlying data points}
}

\seealso{
\code{\link{features}}, \code{\link{plot.features}}, \code{\link[graphics]{plot}}
}


\examples{
# Estimating the smooth and the derivatives of a noisy and discretely sampled function. 
n <- 200
x <- sort(runif(n))
y <- exp(-0.2 * sin(2*pi*x)) + rnorm(n, sd=0.05)

ans <- features(x, y)
fget(ans)

}
