.onUnload <- function (libpath) {
  library.dynam.unload("covRIT", libpath)
}