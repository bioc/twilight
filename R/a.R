.First.lib <- function(libname,pkgname){
  library.dynam("twilight",pkgname,libname)

  library(stats)
  library(splines)
  
  if(.Platform$OS.type=="windows" && require(Biobase) && interactive() && .Platform$Gui=="Rgui"){
    addPDF2Vig("twilight")
  }
  
}
 
