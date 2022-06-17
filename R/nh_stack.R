# nh_stack

#' Stack multiple binary SDM rasters into one raster layer, using terra
#' 
#' Takes a \code{rastfiles} of binary raster's filenames (values 0/1), 
#' a template raster covering the extent desired for the stack, 
#' and optionally the \code{codes} (names) to use
#' for the rasters. Returns a single raster layer (factor type with an
#' attribute table). 
#' 
#' The raster attribute table (accessed using \code{levels()}) 
#' has four columns:  `ID`: the unique integer value in the raster; `VALUE`:
#' the internal nh_stack unique value for that `ID`, `ALLCODES`: the
#' identity of species codes, pasted in a character vector seperated by `;`,
#' `ALLCODES_CT`: the number of unique codes for that value.
#' 
#' If \code{codes} is not given, the raster layer name will be used as the 
#' layer's code. If these are not unique, the internal nh_stack unique value will
#' be pasted to the end of the original code, and a message will be printed.
#' 
#' All rasters must have the same projection and resolution, though
#' they can have different extents - the processing extent is defined by
#' \code{rast}. To summarize all rasters across their
#' entire extents, \code{rast} should essentially be the union (mosaic) 
#' of all raster extents.
#' 
#' When \code{return.table} is TRUE, a list with two objects (1), the 
#' stack raster, and (2) a summary table of included rasters is returned (with
#' internal nh_stack unique values, species codes, and file names). This table
#' is required for resampling the stack (i.e. with \code{nh_stack_resample}).
#' 
#' @param rastfiles raster file names; character vector
#' @param rast raster template dataset
#' @param codes species codes; character vector. If given, must match length of \code{rastfiles}
#' @param return.table Whether to return a table with nh_stack unique values, species codes, and filenames
#' @param clip.feat sf data with masking features for rastfiles. Must have column named 'code', matching \code{codes}
#' 
#' @return SpatRaster
#' 
#' @author David Bucklin
#'
#' @import terra
#' @importFrom stringi stri_rand_strings stri_length stri_sub
#' @importFrom sf st_transform st_crs
#' @export
#'
#' @examples
#' \dontrun{
#' setwd("D:/PSH")
#' rast <- rast("inputs/_masks/data_mask.tif")
#' list <- sort(list.files(paste0("proj_psh/thumb"), full.names = T, pattern = ".tif$"))
#' rastfiles <- c(list[sample(1:length(list), size=5)], "proj_psh/thumb/helevirg_10Aug2018.tif", 
#'    "proj_psh/thumb/myotsoda_10May2018.tif")
#' codes <- unlist(lapply(basename(rastfiles), function(x) strsplit(x, "_")[[1]][1]))
#'
#' stack <- nh_stack(rastfiles, rast, codes = NULL)
#' cats(stack[[1]])
#' }

nh_stack <- function(rastfiles, rast, codes = NULL, return.table = TRUE, clip.feat = NULL) {
  # Works with terra.
  
  list <- rastfiles
  if (!is.null(codes) & length(codes) != length(list)) stop("`codes` length must match `rastfiles` length.")
  
  # load clip.feat
  if (!is.null(clip.feat)) {
    if (!"code" %in% names(clip.feat)) stop("`clip.feat` must have the column 'code' for matching raster codes.")
    message('Prepping clipping features...')
    clip <- st_transform(clip.feat, sf::st_crs(rast))
  }
  
  message("Prepping template raster...")
  if (ncell(rast) > 2.147e+9) stop("Raster template dimensions too large, will need to process in blocks.")
  rcont <- setValues(rast(rast), 1:ncell(rast))
  # dataType(rcont) <- "INT4S"
  
  # make unique codes
  sz <- length(list)
  
  len <- 0
  cd <- character(0)
  while (length(cd) < sz) {
    len <- len + 1
    cd <- unique(stri_rand_strings(sz * 10, len, pattern = "[a-z0-9]"))
  }
  # get codes for length you need
  stk_order <- data.frame(uval = sample(cd, size = sz, replace = F), code_nm = NA, file = list)
  # this is a running list of codes for each cell
  bigd <- rep("",ncell(rcont))
  # loop over layers
  message("Stacking layers...")
  for (i in 1:nrow(stk_order)) { 
    # get nh_stack uval
    spcd <- stk_order$uval[i]
    # get/crop raster
    r.orig <- rast(stk_order$file[i])
    r <- crop(r.orig, rcont)
    # assign long code name
    if (!is.null(codes)) {
      names(r) <- codes[i]
    } else {
      names(r) <- r.orig@ptr$get_sourcenames()
    }
    if (names(r) %in% stk_order$code_nm) {
      message("Non-unique code '", names(r), "' changed to '",paste0(names(r), "_", spcd),"'.")
      names(r) <- paste0(names(r), "_", spcd)
    }
    message(paste0("\nWorking on raster ", names(r), "..."), appendLF = F)
    stk_order$code_nm[stk_order$uval == spcd] <- names(r)
    
    # Clip/mask procedure
    if (!is.null(clip.feat)) {
      c1 <- clip[clip$code==names(r),]
      if (nrow(c1) == 1) {
        message('Clipping...', appendLF = F)
        r <- mask(crop(r, c1), vect(c1), touches = F)  # Terra accepts vectors for masking
      }
    }
    # crop/mask continuous raster, get values
    cr <- crop(rcont, r, datatype = "INT4S")
    rval <- values(r)
    v <- values(cr)[!is.na(rval) & rval==1]
    
    # paste to bigd
    if (length(v) > 0) bigd[v] <- paste0(bigd[v], spcd) 
  }
  
  # Get unique values
  uvals <- data.frame(CODE=unique(bigd))
  uvals$VALUE <- 1:nrow(uvals)  # VALUE must be numeric to work with ArcGIS (will create attribute table)
  bigd2 <- match(bigd, uvals$CODE)
  
  message("\nCalculating stack attributes...")
  r2 <- setValues(rast(rcont), bigd2)
  
  # parse
  pars <- lapply(uvals$CODE, 
                 function(x) if (stri_length(x) > len) stri_sub(x, from = seq(1, stri_length(x), by = len), length = len) else x)
  parssp <- unlist(lapply(pars, function(x) paste(stk_order$code_nm[stk_order$uval %in% x], collapse = ";")))
  parsct <- unlist(lapply(pars, length))
  
  uvals$ALLCODES <- parssp
  uvals$ALLCODES_CT <- parsct
  uvals$ALLCODES_CT[uvals$CODE == ""] <- 0
  
  levels(r2) <- uvals[c("VALUE", "CODE", "ALLCODES", "ALLCODES_CT")]
  names(r2) <- "nh_stack"
  
  if (return.table) {
    names(stk_order) <- c("nh_stack_uval","CODE","FILENAME")
    return(list(r2, stk_order))
  } else {
    return(r2)
  }
}

# nh_stack_resample

#' Resample a raster from nh_stack to a lower (coarser) resolution
#' 
#' Takes an output raster from nh_stack, and returns a lower-resolution
#' version, with recalculated species assemblages for the larger cells.
#' New values are "aggregated" by \code{fact}, the number of cells
#' to aggregate in the x/y dimensions (see \code{?terra::aggregate}).
#' 
#' You can also provide polygons (sp or sf-class) to \code{spf}, over which
#' to aggregate species assemblages. The polygons intersecting areas with
#' data in \code{rast} are returned, with columns identifying species codes 
#' and counts. Polygons will be returned in their original projection, but 
#' processing internally is done in the raster's projection.
#' 
#' @param rast raster output from nh_stack
#' @param lookup lookup table from nh_stack
#' @param fact aggregation factor, in number of cells (see \code{?terra::aggregate})
#' @param spf Optional vector spatial features to use for aggregation (sp or sf-class polygons). If supplied, \code{fact} will be ignored
#' 
#' @return SpatRaster or spf (with attributes added)
#' 
#' @author David Bucklin
#'
#' @import terra
#' @importFrom stringi stri_rand_strings stri_length stri_sub
#' @importFrom methods as
#' @importFrom sf st_crs st_as_sfc st_intersects st_bbox st_intersects st_transform
#' @export
#'
#' @examples
#' \dontrun{
#' # stack <- nh_stack(list, rast, return.table = TRUE)
#' 
#' # resample from 30m to 990m (~1km) resolution
#' stack1km <- nh_stack_resample(stack[[1]], stack[[2]], fact = 33)
#' cats(stack1km)
#' 
#' # view species count raster (index=3 is the 'ALLCODES_CT' column)
#' ct <- as.numeric(stack1km, index=3)
#' plot(ct)
#' }

nh_stack_resample <- function(rast, lookup, fact = 10, spf = NULL) {
  # split length
  len <- unique(nchar(lookup$nh_stack_uval))
  # get levels
  lev <- cats(rast)[[1]]
  # lev <- levels(rast)[[1]]
  names(lev)[1:2] <- c("ID", "category")
  
  # aggregate
  spp <- c()
  
  if (is.null(spf)) {
    # check fact
    if (all(fact < 2)) {
      stop("'fact' must be a one or two-length integer greater than 1.")
    } else {
      a.fact <- round(fact)
    }
    
    message("Aggregating stack...")
    r1 <- aggregate(rast, fact = a.fact, 
                    fun = function(x, ...) {
                      uv <- unique(x)
                      if (all(is.na(uv))) {
                        sp <- NA
                      } else {
                        uv <- uv[!is.na(uv)]
                        cats <- paste(lev$category[lev$ID %in% uv], collapse = "")
                        if (cats != "") {
                          sp <- sort(unique(stri_sub(cats, seq(1, stri_length(cats), by = len), length = len)))
                          sp <- paste(sp, collapse = "")
                        } else {
                          sp <- ""
                        }
                      }
                      if (!is.na(sp)) {
                        if (!sp %in% spp) spp <<- c(spp, sp)
                        val <- match(sp, spp)
                        return(val)
                      } else{
                        return(NA)
                      }
                    })
    # Make levels data frame
    rcats <- data.frame(VALUE = 1:length(spp), CODE = spp)
    levels(r1) <- rcats
    uvals <- cats(r1)[[1]]
  } else {
    proj <- st_crs(spf)
    # handle sp/sf class
    spf1 <- tospf(spf, rast)
    sp <- spf1[[1]]
    spf <- spf1[[2]]
    rm(spf1)
    
    # get intersecting polys
    extr <- st_as_sfc(st_bbox(rast))
    spf <- spf[unlist(lapply(st_intersects(spf, extr), any)),]
    
    # add ID
    spf$burnval <- 1:length(spf$geometry)
    zonr <- gRasterize(spf, rast, value = "burnval")
    
    message("Aggregating stack by polygons...")
    zon <- as.data.frame(zonal(rast, zonr,
                               fun = function(x, ...) {
                                 uv <- unique(x)
                                 if (all(is.na(uv))) {
                                   sp <- NA
                                 } else {
                                   uv <- uv[!is.na(uv)]
                                   cats <- paste(lev$category[lev$ID %in% uv], collapse = "")
                                   if (cats != "") {
                                     sp <- sort(unique(stri_sub(cats, seq(1, stri_length(cats),by = len), length = len)))
                                     sp <- paste(sp, collapse = "")
                                   } else {
                                     sp <- ""
                                   }
                                 }
                                 spp <<- c(spp, sp)
                                 return(1)
                               }))
    zon$code <- spp
    polys2 <- merge(spf, zon, by.x = "burnval", by.y = "burnval")
    polys2$burnval <- NULL
    polys2$nh_stack <- NULL
    
    uvals <- data.frame(CODE = unique(polys2$code)[!is.na(unique(polys2$code))])
  }
  
  message("Calculating stack attributes...")
  # species lookup
  stk_order <- lookup
  # parse
  pars <- lapply(uvals$CODE, 
                 function(x) if (stri_length(x) > len) stri_sub(x, from = seq(1, stri_length(x), by = len), length = len) else x)
  parssp <- unlist(lapply(pars, function(x) paste(stk_order$CODE[stk_order$nh_stack_uval %in% x], collapse = ";")))
  parsct <- unlist(lapply(pars, length))
  
  uvals$ALLCODES <- parssp
  uvals$ALLCODES_CT <- parsct
  uvals$ALLCODES_CT[uvals$CODE == ""] <- 0
  
  if (is.null(spf)) {
    levels(r1) <- uvals
    names(r1) <- "nh_stack_resample"
  } else {
    r1 <- merge(polys2, uvals, by.x = "code", by.y = "CODE")
    r1 <- st_transform(r1, proj)
    if (sp) return(as(r1,Class = "Spatial")) else return(r1)
  }
  return(r1)
}
