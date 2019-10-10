# nh_linefill

#' For line-based SDM predictions, identify areas between predicted suitable
#' lines
#' 
#' Provided a network of lines (e.g., reaches in a hydrological network), and 
#' an attribute that identifies those lines representing suitable habitat, this
#' function identifies fill-in lines between suitable lines, given certain limits
#' (total distance and/or total reaches). It does not alter original line segments.
#' 
#' Line directionality is used to identify how they are connected in the network,
#' so start/end nodes must align. These nodes will be calculated; alternatively, 
#' the names \code{startNode} and \code{endNode} can be provided as columns in \code{spf} and 
#' will be used instead.
#' 
#' Values for \code{max.dist} should be in the units of the CRS, except when
#' using a lat/lon based CRS, in which case the value should be in meters.
#' 
#'
#' @param spf input spatial features (model predictions; sp or sf spatial object)
#' @param field name of field in \code{spf} (binary 1/0), where lines will be filled between features with value == 1
#' @param max.dist maximum distance between selected lines to fill in 
#' @param max.line maximum number of lines between selected lines to fill in
#' 
#' @return sf object including suitable lines and filled-in lines
#' 
#' @author David Bucklin
#'
#' @import sf
#'
#' @export
#'
#' @examples
#' \dontrun{
#' spf <- st_read("SDM_results.shp")
#' fill <- nh_linefill(spf, "thresh", max.dist = 5000, max.line = 5)
#' }

nh_linefill <- function(spf, field, max.dist = NA, max.line = NA) {
  
  if (is.na(max.dist) & is.na(max.line)) stop("One or both of `max.dist` and `max.line` must be specified (not NA).")
  
  spf1 <- tospf(spf)
  sp <- spf1[[1]]
  spf <- spf1[[2]]
  rm(spf1)
  
  # get names from original data, calculate new fields
  nm <- names(spf)
  nm <- nm[nm != "geometry"]
  spf$nh_linefill_id <- as.numeric(row.names(spf))
  if (!is.na(max.dist)) spf$nh_linefill_leng_m <- as.numeric(st_length(spf))
  spf$nh_linefill_field <- eval(parse(text = paste0("spf$", field)))
  spf$nh_linefill_field[is.na(spf$nh_linefill_field)] <- 0

  # make start/endpoints
  if (all(c("startNode", "endNode") %in% names(spf))) {
    message("Using existing node columns...")
  } else {
    spf <- nodes(spf)
  }
  
  # expression for loop
  if (is.na(max.line)) {
    bool <- "any(rleng)"
  } else if (is.na(max.dist)) {
    bool <- "r <= max.line"
  } else {
    bool <- "any(rleng) & r <= max.line"
  }
  
  # this uses distance threshold to add reaches
  ls <- spf$nh_linefill_id[spf$nh_linefill_field==1]
  add <- c()
  # max.r <- NA # reaches downstream allowed
  message("Identifying lines to fill in...")
  for (i in ls) {
    r <- 1
    to <- st_drop_geometry(spf[spf$nh_linefill_id==i,])
    rt <- split(to, seq(nrow(to)))
    rleng <- TRUE
    while (eval(parse(text=bool))) {
      ct <- 0
      rt2 <- list()
      for (n in 1:length(rt)) {
        x <- rt[[n]]
        to <- spf[spf$startNode %in% x$endNode[nrow(x)],]
        to <- st_drop_geometry(to)
        if (nrow(to) > 0) {
          for (nt in 1:nrow(to)) {
            ct <- ct + 1
            rt2[[ct]] <- rbind(x, to[nt,])
          }
        } else {
          ct <- ct + 1
          rt2[[ct]] <- x
        }
      }
      rt <- rt2
      
      # if rows are not added, remove sequence from consideration
      rt <- rt[!unlist(lapply(rt, nrow)) == r]
      if (length(rt) == 0) break
      
      # if another reach is hit, add to a vector `rls`
      rls <- lapply(rt, function(x) {
        c1 <- x$nh_linefill_id[2:nrow(x)]
        if (any(c1 %in% ls)) c1 else NA 
      })
      # check length of reaches
      if (!is.na(max.dist)) {
        rleng <- unlist(lapply(rt, function(x) {
          leng <- sum(x$nh_linefill_leng_m[2:nrow(x)])
          leng <= max.dist 
        }))
        rls <- unlist(rls[rleng])
      } else {
        rls <- unlist(rls)
      }
      rls <- rls[!is.na(rls)]
      
      # if rls is not all NA loop ends. New ids (if any) get added to `add`
      if (length(rls) > 0) {
        add <- c(add, rls[!rls %in% ls])
        break
      }
      # continue to next reach
      r <- r+1
    }
  }
  
  add <- unique(add)
  spf.out <- spf[spf$nh_linefill_id %in% c(ls,add),][nm]
  
  if (!is.null(spf) && sp) return(as(spf.out,Class = "Spatial")) else return(spf.out)
}
