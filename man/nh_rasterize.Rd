% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nh_rasterize.R
\name{nh_rasterize}
\alias{nh_rasterize}
\title{Convert vector format SDM predictions to raster format}
\usage{
nh_rasterize(
  spf,
  rast,
  pred.vals,
  buffer = NULL,
  priority = NULL,
  touches = TRUE,
  rast.out = NULL,
  ...
)
}
\arguments{
\item{spf}{input spatial features (model predictions; sp or sf spatial object)}

\item{rast}{input raster template}

\item{pred.vals}{column name in spf holding feature prediction values}

\item{buffer}{numeric (single value or vector matching length of spf); spatial buffer around spf to include in burn-in}

\item{priority}{column name in spf holding priority values for feature rasterization (higher values have priority)}

\item{touches}{from \code{terra::rasterize}: If TRUE, all cells touched by lines or polygons are affected, not just those on the line render path, 
or whose center point is within the polygon.}

\item{rast.out}{Optional output raster file name (with file extension)}

\item{...}{Additional arguments to terra::writeRaster (e.g. overwrite)}
}
\value{
SpatRaster
}
\description{
Requires spatial features with a prediction values attribute, and a 
template raster. The buffer is optional; when provided, the
buffer distance can be provided as a single numeric value or a 
vector matching length of \code{spf}, for a variable buffer by feature.
}
\details{
Buffer units should be given in the units of \code{rast}. Make 
sure to keep in mind the resolution of the raster when choosing 
a buffer. If rasterizing lines, and a one-cell width is desired,
do not use a buffer.

A vector can be given to `priority` for sorting prior to rasterization, 
where higher values have priority. When \code{priority = NULL} (default),
priority will be defined as the \code{pred.vals}.

If \code{rast.out} is not specified, the raster will remain in temp folder.
}
\examples{
\dontrun{
spf <- st_read("_data/species/acipoxyr/outputs/model_predictions/acipoxyr_20180105_133929_results.shp")
rast <- terra::rast("_data/species/ambymabe/outputs/model_predictions/ambymabe_20171018_130837.tif")

# rasterize
rast_poly <- nh_rasterize(spf, rast, pred.vals = "prbblty", priority = "strord", 
   buffer = spf$strord*15, touches=F, rast.out = "C:/David/scratch/nh_rasterize_poly.tif", 
   overwrite = T)
rast_line <- nh_rasterize(spf, rast, pred.vals = "prbblty", priority = "strord", touches=T, 
   rast.out = "C:/David/scratch/nh_rasterize_line.tif", overwrite = T)
}
}
\author{
David Bucklin
}
