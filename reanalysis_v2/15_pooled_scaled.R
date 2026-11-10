# =====================================================================
# 15  Step 11 again, with the SAME model, draws and likelihood, but with parameter scaling for the optimiser
# =====================================================================
# Why: the parameters of step 11 differ greatly in size (means and spreads of WTAs in tens or hundreds of 100-NOK units, log-cost around -2,
# Cholesky terms in between).  BFGS then takes badly scaled steps and can stop on a flat part of the likelihood ('relative function convergence').
# Apollo's estimate_settings 'scaling' divides each parameter by a factor INSIDE the optimiser only; the model, the likelihood, the draws and the
# returned estimates (in the original units) are unchanged.  Here: factor 10 for every parameter whose starting value exceeds 5 in absolute value, 1 for the others.
# Jobs: U3 | scale3 | cov3 | attr3 | sq3 | all3  (restarts from the best first-round solution 'g', as the jobs ending in 2 in step 11).
SCf <- function(b) setNames(ifelse(abs(b) > 5, 10, 1), names(b))
src <- readLines(file.path(V2, "11_pooled_mixed.R"), encoding = "UTF-8")
src <- sub('restart <- grepl("2$", job); base <- sub("2$", "", job)', 'restart <- grepl("[23]$", job); base <- sub("[23]$", "", job)', src, fixed = TRUE)
i1 <- grep("prob_P, st, fixed", src, fixed = TRUE); i2 <- grep("prob_P, mP$estimate, fixed", src, fixed = TRUE)
stopifnot(length(i1) == 1, length(i2) == 1)
src[i1] <- sub("estimate_settings = list(maxIterations = 400)", "estimate_settings = list(maxIterations = 400, scaling = SCf(st[setdiff(names(st), fixed)]))", src[i1], fixed = TRUE)
src[i2] <- sub("estimate_settings = list(maxIterations = 400)", "estimate_settings = list(maxIterations = 400, scaling = SCf(mP$estimate[setdiff(names(mP$estimate), fixed)]))", src[i2], fixed = TRUE)
stopifnot(sum(grepl("scaling = SCf", src)) == 2)
eval(parse(text = src), envir = globalenv())
