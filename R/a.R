.First.lib <- function(libname,pkgname,where){
  library.dynam(pkgname,pkgname,libname)

  library(stats)
  library(splines)
  
  if(.Platform$OS.type=="windows" && require(Biobase) && interactive() && .Platform$GUI=="Rgui"){
    addVigs2WinMenu("twilight")
  }
  
}
 
