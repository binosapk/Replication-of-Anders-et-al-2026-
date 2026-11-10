# =====================================================================
# 16  Draw-seed check of the pooled mixed logit: the SAME model, number of draws (250 MLHS), start values and scaled optimiser as step 15,
#     but with a different seed for the simulation draws (Apollo's default seed is 13; here 2027).
# =====================================================================
# Why: log-likelihoods and tests of a mixed logit depend slightly on the random draws. If the unrestricted model (U) and the model with everything equal (all)
# give the same conclusions with other draws, the likelihood-ratio test of step 13 does not rest on one lucky set of draws.
# Jobs: U4 (unrestricted) and all4 (everything equal); both restart from the best first-round solution 'g', as in steps 11 and 15.
SEED <- 2027
SCf <- function(b) setNames(ifelse(abs(b) > 5, 10, 1), names(b))
ftxt <- deparse(fit_model)
stopifnot(sum(grepl("noDiagnostics = TRUE", ftxt, fixed = TRUE)) == 1)
ftxt <- sub("noDiagnostics = TRUE", sprintf("noDiagnostics = TRUE, seed = %d", SEED), ftxt, fixed = TRUE)
fit_model <- eval(parse(text = ftxt))
src <- readLines(file.path(V2, "11_pooled_mixed.R"), encoding = "UTF-8")
src <- sub('restart <- grepl("2$", job); base <- sub("2$", "", job)', 'restart <- grepl("[234]$", job); base <- sub("[234]$", "", job)', src, fixed = TRUE)
i1 <- grep("prob_P, st, fixed", src, fixed = TRUE); i2 <- grep("prob_P, mP$estimate, fixed", src, fixed = TRUE)
stopifnot(length(i1) == 1, length(i2) == 1)
src[i1] <- sub("estimate_settings = list(maxIterations = 400)", "estimate_settings = list(maxIterations = 400, scaling = SCf(st[setdiff(names(st), fixed)]))", src[i1], fixed = TRUE)
src[i2] <- sub("estimate_settings = list(maxIterations = 400)", "estimate_settings = list(maxIterations = 400, scaling = SCf(mP$estimate[setdiff(names(mP$estimate), fixed)]))", src[i2], fixed = TRUE)
stopifnot(sum(grepl("scaling = SCf", src)) == 2)
cat(sprintf("Draw seed %d (Apollo default 13)\n", SEED))
eval(parse(text = src), envir = globalenv())
