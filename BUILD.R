setwd("E:/git/nhSDM")
rm(list=ls())
library(devtools)

# use_build_ignore("BUILD.R")
# use_build_ignore("nhSDM.pdf")
use_build_ignore("working")

document()
install()
check()
rm(list=ls())

###
build_pdf <- function() {
  pkg <- as.package(".")
  path <- paste(dirname(pkg$path), noquote(pkg$package),
                sep = "/")
  prev <- ifelse(TRUE, "", " --no-preview")
  over <- ifelse(TRUE, " --force", "")
  cmd <- paste("R CMD Rd2pdf ", shQuote(pkg$path), " -o ",
               path, "/", noquote(pkg$package), ".pdf", prev, over,
               sep = "")
  system(cmd)
}
build_pdf()
