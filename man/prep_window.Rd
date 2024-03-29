% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{prep_window}
\alias{prep_window}
\title{Define available plotting region in Jupyter notebook (and set default graphical pars)}
\usage{
prep_window(rc = c(1, 1), w = 4, h = 4, res = 200, pch = 20, ...)
}
\arguments{
\item{rc}{Vector setting the number of rows and columns. Default is c(1,1)}

\item{w}{Scalar setting the width of each of the (r,c) plotting regions. Default is 4}

\item{h}{Scalar setting the height of each of the (r,c) plotting regions. Default is 4}

\item{res}{Scalar setting the resolution of the plots. Default is 200.}

\item{pch}{Set default plotting character. Default is 20 (small filled circles)}

\item{...}{Additional graphical parameters to be passed to par()}
}
\description{
Define available plotting region in Jupyter notebook (and set default graphical pars)
}
