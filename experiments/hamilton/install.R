## Install ssEI from this checkout into the Hamilton library if the package source changed since the last install.
## Either, from the Hamilton shell (runs as a Slurm job, nothing to watch):   bash experiments/hamilton/sweep.sh install
## or in RStudio with the repository as the working directory:                  source("experiments/hamilton/install.R")
## Compiler output goes to a log file (<library>/ssEI_install.log), not the console: it is thousands of lines.
## Library: SSEI_LIB, default /nobackup/<user>/R-libs. Compiling the Stan model takes about 5-10 minutes and 4+ GB of memory.
source("experiments/hamilton/pkg_hash.R")
lib <- Sys.getenv("SSEI_LIB", file.path("/nobackup", Sys.info()[["user"]], "R-libs")); dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths())); logf <- file.path(lib, "ssEI_install.log")
need <- setdiff(c("quadprog", "coda", "eiPack", "rstan", "rstantools", "StanHeaders", "BH", "RcppEigen", "RcppParallel"), rownames(installed.packages()))
if (length(need)) { message("installing: ", paste(need, collapse = ", ")); install.packages(need, lib = lib, quiet = TRUE) }
h <- pkg_hash("."); stamp <- file.path(lib, "ssEI", "SOURCE_HASH")
if (file.exists(stamp) && identical(readLines(stamp), h) && requireNamespace("ssEI", lib.loc = lib, quietly = TRUE)) {
  message("ssEI is up to date with this checkout (source ", substr(h, 1, 8), "); nothing to do")
} else {
  if ("ssEI" %in% loadedNamespaces()) stop("ssEI is loaded in this session: restart R (Session > Restart R) and source this file again")
  message("building and installing ssEI (source ", substr(h, 1, 8), ") into ", lib, "; compiler output -> ", logf, " ...")
  repo <- normalizePath("."); bd <- tempfile("ssei_build"); dir.create(bd); Rbin <- file.path(R.home("bin"), "R")
  old <- setwd(bd); on.exit(setwd(old), add = TRUE)
  st1 <- system2(Rbin, c("CMD", "build", "--no-build-vignettes", "--no-manual", shQuote(repo)), stdout = logf, stderr = logf)
  tarball <- list.files(bd, pattern = "^ssEI_.*\\.tar\\.gz$", full.names = TRUE)
  st2 <- if (st1 == 0 && length(tarball) == 1) system2(Rbin, c("CMD", "INSTALL", "--no-multiarch", paste0("--library=", shQuote(lib)), shQuote(tarball)),
                                                       stdout = paste0(logf, ".install"), stderr = paste0(logf, ".install")) else 1
  if (file.exists(paste0(logf, ".install"))) { cat(readLines(paste0(logf, ".install")), file = logf, sep = "\n", append = TRUE); unlink(paste0(logf, ".install")) }
  setwd(old)
  if (st2 != 0) { lg <- readLines(logf); message("INSTALL FAILED; last lines of ", logf, ":"); writeLines(tail(lg, 25)); stop("ssEI was not installed") }
  writeLines(h, stamp); message("done: ", R.version.string, ", source ", substr(h, 1, 8), " (full log: ", logf, ")")
}
