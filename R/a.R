.First.lib <- function(libname,pkgname,where){
  library.dynam("twilight",pkgname,libname)
  
  if(.Platform$OS.type=="windows" && require(Biobase) && interactive() && .Platform$GUI=="Rgui"){
    addVigs2WinMenu("twilight")
  }
  
}
 
