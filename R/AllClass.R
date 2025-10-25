#' lmeb class
#'
#'

# setClass("lmeb", representation = list(relfac = "list"),
#          contains = "merMod")
# setClass("glmerlmebreed", representation = list(resp="glmResp"),
#          contains = "lmeb")
# setClass("lmerlmebreed", representation = list(resp="lmerResp"),
#          contains = "lmeb")

setClass("lmeb", slots = c(relfac = "list", udu = "list"),
         contains = c("merMod"))
setClass("glmerMod", slots = c(resp="glmResp"),
         contains = "lmeb")
setClass("lmerMod", slots = c(resp="lmerResp"),
         contains = "lmeb")

# setClass("lmeb", representation = list(relfac = "list", udu="list"),
#          contains = "merMod")
# setClass("glmerlmebreed", representation = list(resp="glmResp"),
#          contains = "lmeb")
# setClass("lmerlmebreed", representation = list(resp="lmerResp"),
#          contains = "lmeb")
