#' lmebreed class
#'
#'

setClass("lmebreed", representation = list(relfac = "list"),
         contains = "merMod")
setClass("glmerlmebreed", representation = list(resp="glmResp"),
         contains = "lmebreed")
setClass("lmerlmebreed", representation = list(resp="lmerResp"),
         contains = "lmebreed")




