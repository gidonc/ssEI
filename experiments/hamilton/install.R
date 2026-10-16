## Install ssEI from this checkout into the Hamilton library if the package source changed since the last install.
## Run in RStudio Server with the repository as the working directory, after a fresh R session (Session > Restart R):
##   source("experiments/hamilton/install.R")
## Library: SSEI_LIB, default /nobackup/<user>/R-libs. Compiling the Stan model takes about 5-10 minutes.
source("experiments/hamilton/pkg_hash.R")
lib <- Sys.getenv("SSEI_LIB", file.path("/nobackup", Sys.info()[["user"]], "R-libs")); dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))
need <- setdiff(c("remotes", "quadprog", "coda", "eiPack", "rstan", "rstantools", "StanHeaders", "BH", "RcppEigen", "RcppParallel"), rownames(installed.packages()))
if (length(need)) { message("installing: ", paste(need, collapse = ", ")); install.packages(need, lib = lib) }
h <- pkg_hash("."); stamp <- file.path(lib, "ssEI", "SOURCE_HASH")
old <- if (file.exists(stamp)) readLines(stamp) else ""
if (identical(old, h) && requireNamespace("ssEI", lib.loc = lib, quietly = TRUE)) {
  message("ssEI is up to date with this checkout (", substr(h, 1, 8), "); nothing to do")
} else {
  if ("ssEI" %in% loadedNamespaces()) stop("ssEI is loaded in this session: restart R (Session > Restart R) and source this file again")
  message("installing ssEI from this checkout (source ", substr(h, 1, 8), ") into ", lib, " ...")
  remotes::install_local(".", lib = lib, force = TRUE, upgrade = "never", dependencies = FALSE)
  writeLines(h, stamp); message("done: ", R.version.string, ", source ", substr(h, 1, 8))
}
