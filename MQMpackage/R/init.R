######################################################################
#
# init.R
#
# copyright (c) 2009, Danny Arends
# written Jan, 2009
# .First.lib is executed when the package is loaded with library(MQMpackage)
#
######################################################################

.First.lib <- function(lib, pkg) library.dynam("MQMpackage", pkg, lib)

# end of init.R
