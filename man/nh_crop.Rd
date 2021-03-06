% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nh_crop.R
\name{nh_crop}
\alias{nh_crop}
\title{Crop extra rows/columns from all sides of a raster}
\usage{
nh_crop(rast)
}
\arguments{
\item{rast}{input raster}
}
\value{
RasterLayer
}
\description{
This function will reduce the extent of a raster, by removing rows/columns
from all sides if they do not have any non-NA cells with a value greater than zero.
One buffer row/column is left on each side.
}
\examples{
\dontrun{
rast <- rast <- raster("ambymabe_16May2018.tif")
rast.crop <- nh_crop(rast)

# should be TRUE
sum(values(rast), na.rm=T) == sum(values(rast.crop), na.rm=T)
}
}
\author{
David Bucklin
}
