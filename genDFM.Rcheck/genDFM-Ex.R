pkgname <- "genDFM"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('genDFM')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("genDFM-package")
### * genDFM-package

flush(stderr()); flush(stdout())

### Name: genDFM-package
### Title: A short title line describing what the package does
### Aliases: genDFM-package genDFM
### Keywords: package

### ** Examples

  ## Not run: 
##D      ## Optional simple examples of the most important functions
##D      ## These can be in \dontrun{} and \donttest{} blocks.   
##D   
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
