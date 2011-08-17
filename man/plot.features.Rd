\name{plot}
\alias{plot}
\alias{plot.features}
\title{
Plotting the smoothed fucntion and its derivatives.
}
\description{Plots of the smoothed function and its first 2 derivatives.}

\usage{
\method{plot}{features}(x, ...)
}

\arguments{
\item{x}{an object of \code{class(features)} that is returned by the \code{features()} function.}
\item{...}{Additional arguments to be passed to the \code{plot()} function}
}


\value{
A plot is shown in a 2 x 2 layout, where the top 2 frames are the same and they depict the raw data and the smoothed function.  
The bottom left panel shows the smoothed first derivative, and the bottom right panel depicts the smoothed second derivative.
}


\examples{
# Estimating the smooth and the derivatives of a noisy and discretely sampled function. 
n <- 200
x <- sort(runif(n))
y <- exp(-0.2 * sin(2*pi*x)) + rnorm(n, sd=0.05)

ans <- features(x, y)
plot(ans)

ans.sm <- features(x, y, smoother="smooth.spline")
plot(ans.sm)

}
