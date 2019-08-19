.onLoad <- function(libname = find.package("TraitEnvMLMWA"), pkgname = "TraitEnvMLMWA"){
  # CRAN Note avoidance
  if (getRversion() >= "2.15.1")
    utils::globalVariables(  c("CWM", "SNC", "E", "R", "K",  "N", "N2", "b_eslope", "b_star", "b_tslope", "c_star",
        "meanE", "meanT", "wmeanCWM", "wmeanSNC","Env","Traits","association"))
  invisible()
}
