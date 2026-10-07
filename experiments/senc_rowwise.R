## senc (eiPack), soft multinomial, row-by-row margin model.
##
## Reproduces the configuration used for the senc results in the findings notes:
## rows conditioned on the data, expected table allocated row by row, whitened
## column-margin parameters, volume scaling, interior parameters anchored on E_rc
## and scaled by the hierarchy.
##
## The helpers (experiments/rowwise_helpers.R) build the extra model inputs; they are meant to move into the package.
##
## Usage, from the package directory: source("experiments/senc_rowwise.R"), after installing ssEI and eiPack.
## Set ALLOC <- 2L for the odds-ratio logit allocation instead of adjusted row (3).

source("experiments/rowwise_helpers.R")
ALLOC <- 3L
SEED <- 1234; CHAINS <- 4; ITER <- 1000; WARMUP <- 500

## ---------------------------------------------------------------------------
## Data and basis
## ---------------------------------------------------------------------------
if (!exists("senc")) data(senc, package = "eiPack")
rown <- c("white", "black", "natam"); coln <- c("dem", "rep", "non")
rm <- senc[rown]; cm <- senc[coln]
J <- nrow(rm); R <- 3; C <- 3
kc <- array(0, c(J, R, C))                                # the true cells, used only to score the result
cellnm <- outer(c("wh", "bl", "natam"), coln, paste0)
for (r in 1:R) for (c in 1:C) kc[, r, c] <- senc[[cellnm[r, c]]]
dimnames(kc) <- list(NULL, rown, coln)
stopifnot(all(apply(kc, c(1, 2), sum) == as.matrix(rm)), all(apply(kc, c(1, 3), sum) == as.matrix(cm)))

## within each row: {dem, rep} against {non}, then dem against rep; row margins first
Vr <- rw_sbp_to_v(cbind(c(1, 1, -1), c(1, -1, 0)))
V_within <- matrix(0, R * C, R * (C - 1))
for (r in 1:R) V_within[(r - 1) * C + 1:C, (r - 1) * (C - 1) + 1:(C - 1)] <- Vr
V_ilr <- cbind(kronecker(make_helmert_basis(R), matrix(1 / sqrt(C), C, 1)), V_within)
stopifnot(all.equal(crossprod(V_ilr), diag(R * C - 1), tolerance = 1e-8))

## ---------------------------------------------------------------------------
## Fit and report (error index, coverage, E_rc and sigma summaries, sampler diagnostics)
## ---------------------------------------------------------------------------
fit <- rw_fit(kc, V_ilr, alloc = ALLOC, rem_col = c(1L, 1L, 1L),
              chains = CHAINS, iter = ITER, warmup = WARMUP, seed = SEED)
