# =====================================================================
# 08  Bayesian (hierarchical) estimation of the mixed logit with a FULL 12 x 12 covariance matrix
# =====================================================================
# Leeds [A] 3A.  Footnote 16 of the article: the fully correlated model "did not converge" with simulated maximum likelihood.
# Bayesian estimation (Gibbs + Metropolis-Hastings, package RSGHB through Apollo) does not need to maximise a likelihood, so a full
# covariance matrix is feasible.  Remember (Leeds 1B, slide 3): "HB is NOT a model" - the model is the same mixed logit; only the estimator differs.
# Units: coefficients in 1,000 NOK (so that the Metropolis-Hastings step sizes work; slide 69-70 warn about scaling).
#   V_dev = b_cost * ( sum_k w_k x_k + cost / 1000 ),   V_sq = b_cost * w_sq,   w_k ~ Normal (full covariance),  b_cost ~ positive lognormal
source(file.path(V2, "a3_published.R"))
d <- dat[[SAMPLE]]; d$agegrp <- as.character(d$agegrp)
nm <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view")
database <- d
quiet(apollo_initialise())
apollo_control <- list(modelName = paste0("hb_", SAMPLE), indivID = "id", HB = TRUE, nCores = 1, outputDirectory = outdir)
apollo_beta <- c(setNames(c(3.9, -0.9, -1.4, -0.8, -1.7, -7.4, -1.8, -4.4, -2.6, -0.6, -4.7), nm), b_cost = 1)
apollo_fixed <- c()
apollo_HB <- list(hbDist = c(setNames(rep("N", 11), nm), b_cost = "LN+"), gNCREP = 40000, gNEREP = 20000, gINFOSKIP = 5000, gFULLCV = TRUE, gSeed = 2026)
apollo_inputs <- quiet(apollo_validateInputs())
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs); on.exit(apollo_detach(apollo_beta, apollo_inputs))
  V <- list()
  V[["dev"]] <- b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
                          w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
                          w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 1000)
  V[["sq"]] <- b_cost * w_sq
  P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
  P <- apollo_prepareProb(P, apollo_inputs, functionality); return(P)     # no apollo_panelProd under HB
}
t0 <- Sys.time()
hb <- quiet(apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs))
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("Bayesian mixed logit, %s sample: %d burn-in + %d retained iterations (%.0f min)\n", SAMPLE, 40000, 20000, mins))

rn <- hb$HB_names_random_params
M  <- hb$HB_iterations_means; cn <- colnames(M)
cat("columns of the stored chain of means:", paste(cn, collapse = ", "), "\n")
K  <- hb$HB_iterations_covar
ck <- sub("^.*_", "", ""); invisible(ck)
pick <- function(x) M[, grep(paste0("^", x, "$"), cn)[1]]
post <- sapply(c(nm, "b_cost"), function(v) { j <- which(cn == v | cn == paste0("mu_", v) | cn == paste0(v, "_mean")); if (length(j)) M[, j[1]] else rep(NA, nrow(M)) })
if (anyNA(post)) post <- M[, 1:12, drop = FALSE]
colnames(post) <- c(nm, "b_cost")[seq_len(ncol(post))]
sg <- c(w_sq = 1000, setNames(rep(-1000, 10), nm[-1]))
tab <- data.frame(attribute = nm, WTA_post_mean = round(sg * colMeans(post[, nm])), post_sd = round(1000 * apply(post[, nm], 2, sd)),
                  q025 = round(apply(sg * post[, nm], 2, function(z) quantile(z, if (TRUE) 0.025))), q975 = round(apply(sg * post[, nm], 2, function(z) quantile(z, 0.975))), row.names = NULL)
# for the negative-signed attributes the percentile order flips; recompute properly
lo <- sapply(nm, function(v) min(quantile(sg[[v]] * post[, v], c(.025, .975)))); hi <- sapply(nm, function(v) max(quantile(sg[[v]] * post[, v], c(.025, .975))))
tab$q025 <- round(lo); tab$q975 <- round(hi)
a <- A3[[SAMPLE]]; sg_a <- c(100, rep(-100, 10)); tab$Anders_published <- round(sg_a * a$mean[1:11])
cat("\nPosterior mean of the MEAN WTA (NOK per household-year), posterior sd and 95% credible interval; Anders's published value for comparison\n")
print(tab, row.names = FALSE)

# population standard deviations and correlations of the random coefficients
cv <- apply(K, c(1, 2), mean); sdv <- sqrt(diag(cv)); names(sdv) <- rn
cat("\nPosterior mean standard deviation across respondents (NOK):\n"); print(round(1000 * sdv[1:11]))
cr <- cov2cor(cv); dimnames(cr) <- list(rn, rn)
cat("\nPosterior mean correlations of the status quo constant with the other random coefficients (and with the cost / scale coefficient):\n")
print(round(cr[rn[1], -1], 2))
cat("\nPosterior mean correlation matrix, attribute block (rounded):\n"); print(round(cr[rn[1:11], rn[1:11]], 1))
g <- hb$HB_Geweke_test_means
cat("\nConvergence: Geweke z of the chains of the means (|z| < 2 is fine):\n"); print(round(g, 2))

