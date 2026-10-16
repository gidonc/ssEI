## Pre-flight for sweep.sh: ssEI loads, and was installed from the current package source. Args: repository path.
repo <- commandArgs(TRUE)[1]; source(file.path(repo, "experiments/hamilton/pkg_hash.R"))
ok <- suppressWarnings(suppressMessages(requireNamespace("ssEI", quietly = TRUE)))
if (!ok) { cat("ssEI does not load with R", getRversion(), "and library", .libPaths()[1], "\n  -> in RStudio: source(\"experiments/hamilton/install.R\")\n"); quit(status = 1) }
for (p in c("eiPack", "quadprog", "coda")) if (!requireNamespace(p, quietly = TRUE)) { cat("missing package", p, " -> source(\"experiments/hamilton/install.R\")\n"); quit(status = 1) }
stamp <- file.path(dirname(find.package("ssEI")), "ssEI", "SOURCE_HASH"); h <- pkg_hash(repo)
if (!file.exists(stamp) || !identical(readLines(stamp), h)) {
  cat("the installed ssEI does not match this checkout's package source\n  -> in RStudio: Session > Restart R, then source(\"experiments/hamilton/install.R\")\n"); quit(status = 1) }
cat("ssEI OK (R", as.character(getRversion()), ", source", substr(h, 1, 8), ")\n")
