## md5 of the package source (R/, src/, inst/stan/, DESCRIPTION, NAMESPACE) in a checkout: used to tell whether the installed ssEI is current
pkg_hash <- function(repo = ".") {
  fs <- c(file.path(repo, c("DESCRIPTION", "NAMESPACE")), list.files(file.path(repo, c("R", "src", "inst/stan")), recursive = TRUE, full.names = TRUE))
  fs <- sort(fs[file.exists(fs) & !grepl("\\.(o|so|dll)$|stanExports_.*\\.(cc|h)$", fs)])
  tf <- tempfile(); on.exit(unlink(tf)); writeLines(paste(basename(dirname(fs)), basename(fs), tools::md5sum(fs)), tf)
  unname(tools::md5sum(tf)) }
