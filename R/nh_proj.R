# nh_newproj

#' Initiate new SDM project
#' 
#' Creates directory structure (if it doesn't exist) for SDMs and 
#' new SDM project SQLite database. Folder structure is initiated in current
#' working directory, unless \code{folder} is specified.
#' 
#' More details
#' 
#' @param proj.name The name of the project. Appended to a new project database file.
#' @param folder The SDM folder in which to initiate the project. Defaults to `getwd()`.
#' 
#' @return TRUE on succesful project creation
#'
#' @author David Bucklin
#'
#' @import RSQLite
#'
#'
#' @examples
#' \dontrun{
#' setwd("D:/")
#' nh_newproj("invert")
#' }

nh_newproj <- function(proj.name, folder = "nhSDM") {
  
  s <- folder
  if (!dir.exists(s)) dir.create(s)
  # create folders
  dbfile <- paste0(s, "/databases/", "sdm_tracking_", proj.name, ".sqlite")
  if (file.exists(dbfile)) stop("Database file '", dbfile, "' already exists. Choose a different name or use nh_updateproj() to update database.")
  
  folders <- c("databases", "env_vars/background", "env_vars/raster", "env_vars/tabular",
               "other_spatial/raster", "other_spatial/feature",
               "species/_template_copyme/inputs/background", "species/_template_copyme/inputs/presence",
               "species/_template_copyme/inputs/scripts","species/_template_copyme/inputs/model_input",
               "species/_template_copyme/outputs")
  
  for (i in folders) {
    if (!dir.exists(paste0(s, "/" , i))) dir.create(paste0(s, "/" , i), recursive = T)
  }
  
  sql <- paste(readLines("E:/git/nhSDM/data/sql/sqlite_template_db_nodata.sql"), collapse = "")
  sql <- paste(strsplit(sql, ";")[[1]], ";",sep = "")
  tryCatch(expr = {
      conn <- dbConnect(SQLite(), dbfile)
      dbBegin(conn)
      for (i in sql) dbExecute(conn, i)
    },
    finally = {
      dbCommit(conn)
      dbDisconnect(conn)
    })
  message("New project database '", dbfile , "' created.")
}