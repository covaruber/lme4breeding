#' lmebreed class
#'
#'

# setClass("lmebreed", representation = list(relfac = "list"),
#          contains = "merMod")
# setClass("glmerlmebreed", representation = list(resp="glmResp"),
#          contains = "lmebreed")
# setClass("lmerlmebreed", representation = list(resp="lmerResp"),
#          contains = "lmebreed")

setClass("lmebreed", slots = c(relfac = "list", udu = "list"),
         contains = c("merMod"))
setClass("glmerMod", slots = c(resp="glmResp"),
         contains = "lmebreed")
setClass("lmerMod", slots = c(resp="lmerResp"),
         contains = "lmebreed")

# setClass("lmebreed", representation = list(relfac = "list", udu="list"),
#          contains = "merMod")
# setClass("glmerlmebreed", representation = list(resp="glmResp"),
#          contains = "lmebreed")
# setClass("lmerlmebreed", representation = list(resp="lmerResp"),
#          contains = "lmebreed")
