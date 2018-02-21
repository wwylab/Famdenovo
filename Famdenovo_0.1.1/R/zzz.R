.onLoad <- function(libname, pkgname) {
  data("LFSpro.cancer.type", package="Famdenovo" , envir=parent.env(environment()))
  data("LFSpenet.2010", package="Famdenovo", envir=parent.env(environment()))
}