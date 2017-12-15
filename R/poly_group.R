# poly_group

#' Group spatial polygons using a defined separation distance
#'
#' Input polygons are allocated into new groups depending
#' on whether they lie within the seperation distance of
#' another input polygon.
#'
#' The grouping is done by input ID, so multi-polygons
#' are considered as one input polygon (they will not be
#' be ungrouped.) Grouping is greedy - for a given input
#' polygon, it will belong to the group containing all of
#' its neighbor polygons, plus all of their neighboring
#' polygons, etc.
#'
#' A column 'group' will be added to the output SpatialPolygonsDataFrame.
#' Specifying \code{union = TRUE} will output one (multi)polygon per group,
#' meaning original attributes will be discarded.
#'
#' @param spP input SpatialPolygonsDataFrame
#' @param sep.dist separation distance (in projection units) with which to define groups
#' @param union whether to union output groups in the resulting SPDF
#'
#' @author David Bucklin
#'
#' @importFrom sp merge SpatialPolygonsDataFrame
#' @importFrom maptools unionSpatialPolygons
#' @importFrom rgeos gDistance
#'
#' @export

poly_group <- function(spP, sep.dist, union = FALSE) {

  spP$TEMPRN <- row.names(spP)
  gd <- as.data.frame(gDistance(spP, spP, byid = TRUE))
  # remove comparisons to each other
  for (n in row.names(gd)) gd[n,n] <- NA

  # set up group list
  grps <- as.list(rep(NA, length(colnames(gd))))
  names(grps) <- colnames(gd)

  # group
  for (n in colnames(gd)) {
    if (is.na(grps[n])) grps[n] <- n
    d <- gd[n]
    if (min(d, na.rm = T) < sep.dist) {
      subg <- subset(d, d[n] < sep.dist)
      if (all(is.na(grps[row.names(subg)]))) {
        grps[row.names(subg)] <- n
      } else {
        grps[n] <- min(unlist(grps[row.names(subg)]), na.rm = T)
        grps[row.names(subg)] <- min(unlist(grps[row.names(subg)]), na.rm = T)
      }
    }
  }

  polyg <- data.frame(TEMPRN = names(grps), group = unlist(grps))
  spP <-merge(spP, polyg, by = "TEMPRN")
  spP$group <- as.numeric(spP$group)
  spP$TEMPRN <- NULL

  if (union) {
    spP <- unionSpatialPolygons(spP, spP$group)
    spP <- SpatialPolygonsDataFrame(spP, data = data.frame(group=row.names(spP), row.names = row.names(spP)), match.ID = TRUE)
  }

  return(spP)
}
