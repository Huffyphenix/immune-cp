if(!require("GEOquery")) {
  BiocInstaller::biocLite('GEOquery')
  library(GEOquery)
} else {
  library(GEOquery)
}


gse_accs <- c("GSE2435", "GSE31040", "GSE100020", "GSE90444")
geo_path <- "/home/cliu18/liucj/projects/6.autophagy/GEO"

tryCatch({
  gset <- getGEO(gse_acc[1], GSEMatrix = TRUE, destdir = geo_path)
}, error = function(error) {
  # Try downloading again from scratch...
  tryCatch({
    gset <- getGEO(argv$accession, GSEMatrix = TRUE)
  }, error = function(e) {
    cat("ERROR: Unable to generate gset", file = stderr())
    cat(e, file = stderr())
    quit(save = "no",
         status = 2,
         runLast = FALSE)
  })
})
