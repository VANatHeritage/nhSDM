setwd("D:/git/nhSDM")
rm(list=ls())
library(devtools)

# use_build_ignore("BUILD.R")
# use_build_ignore("nhSDM.pdf")
use_build_ignore("working")
use_build_ignore("recycling")

document()
check()
install()
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
build_pdf()  # fixed error using: https://felixfan.github.io/Font-ts1-zi4r-at-540-not-found/
