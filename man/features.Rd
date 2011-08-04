\name{features}
\alias{features}
\title{
Extracting features of a discretely-sampled functional data.
}
\description{
Discretely-sampled function is first smoothed.  Features of the smoothed function are then extracted.  Some of the key features include mean value, first and second derivatives, critical points 
(i.e. local maxima and minima), curvature of cunction at critical points, wiggliness of the function, noise in data, and outliers in data.}
\usage{
features(x, y, smoother=c("glkerns", "smooth.spline"), fits.return=TRUE, 
	control = list(), ...)
}

\arguments{
\item{x}{a vector of independent variable (e.g. time) at which the function is sampled.}
\item{y}{a vector of response values to be smoothed }
\item{smoother}{a character string specifying the name of the smoothing algorithm; currently only 2 smoothers are implemented: smoothing spline and a kernel smoother with a global plug-in bandwidth}
\item{fits.return}{a logical variable specifying whether or not to return the smoothed function and its first 2 derivatives.}
\item{control}{A list of control parameters.  See *Details* for more information.}
\item{...}{Additional arguments to be passed to the smoothers.  An important optional parameter is the \code{bandwidth} that controls the smoothness of the fit.  See the help for `glkerns' for mode details.}
}

\details{
Argument \code{control} is a list specifing any changes to default values of control parameters.  Note that the names of these must be specified completely.  Partial matching will not work.  

Default values of \code{control} are:  ctrl <- list(npts=100, plot.it=FALSE, c.outlier=3, decim.out=2)

\code{npts}: an integer. Number of points to use in estimating smoothed function and for computing features

\code{plot.it}: a logical variable indicating whether the smoothed function should be plotted

\code{c.outlier}: a constant denoting number of standard deviations away from smooth fit for determining whether a point is an outlier. Default is TRUE.

\code{decim.out}: number of decimals to display in the features output. 
}

\value{
A list with a number of extracted features of the underlying smooth function:
\item{f}{: a numeric vector containing 10 basic features of the function: mean, min, max, std.dev, noise, signa-to-noise ratio, minimum and maximum of first derivative, wiggliness, and number of critical points}

\item{cpts}{: Locations of the critical points (e.g. points where the first derivative is zero)}

\item{curvature}{: Value of second derivative at the critical points}

\item{outliers}{: Locations of outlying data points}

\item{fits}{an attribute of the object returned by \code{features} when \code{return.fits=TRUE}.  A list with 4 vectors: abscissae, smoothed function, first derivative, and second derivative}
}


\examples{
# Estimating the smooth and the derivatives of a noisy and discretely sampled function. 
n <- 200
x <- sort(runif(n))
y <- exp(-0.2 * sin(2*pi*x)) + rnorm(n, sd=0.05)

ans <- features(x, y, fits.return=TRUE, control=list(plot.it=TRUE))
ans$f
ans$f["fwiggle"]
ans$out

fits <- attr(ans, "fits")
# fits$fn is the function value at locations fits$x
# fits$d1 is the first derivative at locations fits$x
# fits$d2 is the second derivative at locations fits$x

par(ask=TRUE)
plot(fits$x, fits$d1, type="l", xlab="x", ylab="First derivative")
yexact <- -0.2 * 2*pi * cos(2*pi*x) * exp(-0.2 * sin(2*pi*x)) 
lines(x, yexact, col=2)
legend(x=0, y=0.98*max(fits$d1), legend=c("Estimated", "Exact"), lty=1, col=1:2)

ans2 <- features(x, y, smoother="smooth.spline", fits.return=TRUE, 
	control=list(plot.it=TRUE))
fits2 <- attr(ans2, "fits")

par(ask=TRUE)
plot(fits2$x, fits2$d1, type="l", xlab="x", ylab="First derivative")
yexact <- -0.2 * 2*pi * cos(2*pi*x) * exp(-0.2 * sin(2*pi*x)) 
lines(x, yexact, col=2)
legend(x=0, y=0.98*max(fits2$d1), legend=c("Estimated", "Exact"), lty=1, col=1:2)

# Passing arguments to smoothers; especially important parameter is "bandwidth"

ans3 <- features(x, y, fits.return=TRUE, control=list(plot.it=TRUE), bandwidth=0.2)
fits3 <- attr(ans3, "fits")

par(ask=TRUE)
plot(fits3$x, fits3$d1, type="l", xlab="x", ylab="First derivative")
yexact <- -0.2 * 2*pi * cos(2*pi*x) * exp(-0.2 * sin(2*pi*x)) 
lines(x, yexact, col=2)
legend(x=0, y=0.98*max(fits3$d1), legend=c("Estimated", "Exact"), lty=1, col=1:2)
}
