% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colour.R
\name{nr_col_pal}
\alias{nr_col_pal}
\title{Colour palettes}
\usage{
nr_col_pal(palette_name, n)
}
\arguments{
\item{palette_name}{One of \code{"colblind"} (max 9), \code{"bigpal"} (max 29)}

\item{n}{Number of colours to return}
}
\value{
A vector of n colours
}
\description{
Colour palettes
}
\examples{
mycol <- nr_col_pal("colblind", 9)
col_show(mycol)
}
\seealso{
\code{\link{col_show}}, \code{\link[grDevices]{colorRamp}}
}
