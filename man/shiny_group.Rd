% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_group.R
\name{shiny_group}
\alias{shiny_group}
\title{Interactively view (in a leaflet map) and select grouping distance}
\usage{
shiny_group(spf)
}
\arguments{
\item{spf}{input spatial features (sp or sf spatial object)}
}
\description{
Interactively view (in a leaflet map) and select grouping distance
}
\examples{
\dontrun{
spf <- readOGR("D:/SDM/Tobacco/inputs/species/ambymabe/polygon_data", "ambymabe_expl")
shiny_group(spf)
}
}
