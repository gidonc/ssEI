## md5 of the package source (R/, src/, inst/stan/, DESCRIPTION, NAMESPACE) in a checkout: used to tell whether the installed ssEI is current.
## Paths are taken relative to the checkout, so the result does not depend on how the checkout's path is written.
pkg_hash <- function(repo = ".") {
  repo <- normalizePath(repo)
  rel <- c("DESCRIPTION", "NAMESPACE", unlist(lapply(c("R", "src", "inst/stan"), function(d) file.path(d, list.files(file.path(repo, d), recursive = TRUE)))))
  rel <- sort(rel[file.exists(file.path(repo, rel)) & !grepl("\\.(o|so|dll)$|stanExports_.*\\.(cc|h)$", rel)])
  tf <- tempfile(); on.exit(unlink(tf)); writeLines(paste(rel, tools::md5sum(file.path(repo, rel))), tf)
  unname(tools::md5sum(tf)) }
