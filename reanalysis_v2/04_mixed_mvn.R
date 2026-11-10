# =====================================================================
# 04  Mixed logit with a FULL covariance matrix (Cholesky), random cost and exogenous shifters  -  estimated in STAGES
# =====================================================================
# Leeds [A] 1E: (i) fixing the cost coefficient imposes homogeneous scale; (ii) uncorrelated coefficients overstate heterogeneity
# if scale varies; (iii) a Cholesky factor gives the full covariance of any set of Normals.  Anders restricts correlation to levels of the same
# attribute and writes (footnote 16) that a fully correlated model did not converge.
# Five random dimensions in WTA space with a FULL 5 x 5 Cholesky factor:  status quo, no protection, < 1 km, visible (Normal) and the log of the
# cost / scale coefficient (Normal on the log scale).  Other coefficients are fixed.  The status quo constant depends on exogenous characteristics.
# A first attempt that started the cost dimension at the published lognormal values ran away to an invalid solution
# (log-likelihood -2232 vs -2065 published, mean of log cost 50).  Leeds [A] 1B (slide 93-96, starting values and local optima): build up in STAGES:
#   stage 1  four Normals, full 4 x 4 Cholesky, cost fixed (b_cost estimated),             covariates on the status quo
#   stage 2  add the cost dimension, starting at the stage-1 solution with a small, free log-cost spread
prob_M <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs); on.exit(apollo_detach(apollo_beta, apollo_inputs))
  sq <- b_sq + g_age * age_c + g_fem * female + g_univ * univ + g_inc * inc_c + g_incm * inc_missing + g_exp * exposed + g_org * org
  V <- list()
  V[["dev"]] <- b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
                          b_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + b_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
                          w_d3 * (dist_alt1 == 3) + b_view * (view_alt1 == 2) + cost_alt1 / 100)
  V[["sq"]] <- b_cost * sq
  P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
  P <- apollo_panelProd(P, apollo_inputs, functionality); P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality); return(P)
}
# stage 1: cost enters as a fixed parameter exp(m_cost) (no draw), stage 2: with a draw on the 5th dimension
rc_M1 <- function(apollo_beta, apollo_inputs) {
  r <- list()
  r[["b_sq"]]   <- m_sq   + l11 * xi1
  r[["b_nop"]]  <- m_nop  + l21 * xi1 + l22 * xi2
  r[["b_d1"]]   <- m_d1   + l31 * xi1 + l32 * xi2 + l33 * xi3
  r[["b_view"]] <- m_view + l41 * xi1 + l42 * xi2 + l43 * xi3 + l44 * xi4
  r[["b_cost"]] <- exp(m_cost) + 0 * xi5
  return(r)
}
rc_M2 <- function(apollo_beta, apollo_inputs) {
  r <- list()
  r[["b_sq"]]   <- m_sq   + l11 * xi1
  r[["b_nop"]]  <- m_nop  + l21 * xi1 + l22 * xi2
  r[["b_d1"]]   <- m_d1   + l31 * xi1 + l32 * xi2 + l33 * xi3
  r[["b_view"]] <- m_view + l41 * xi1 + l42 * xi2 + l43 * xi3 + l44 * xi4
  r[["b_cost"]] <- exp(m_cost + l51 * xi1 + l52 * xi2 + l53 * xi3 + l54 * xi4 + l55 * xi5)
  return(r)
}
dr_M <- list(interDrawsType = "mlhs", interNDraws = 250, interNormDraws = paste0("xi", 1:5))
source(file.path(V2, "a3_published.R"))
a <- A3[[SAMPLE]]; av <- setNames(a$mean, a$par)
cn <- c("l11", "l21", "l22", "l31", "l32", "l33", "l41", "l42", "l43", "l44", "l51", "l52", "l53", "l54", "l55")
sf <- c("w_el38", "w_el57", "w_pl3", "w_pl6", "w_wil", "w_d2", "w_d3"); gg <- c("g_age", "g_fem", "g_univ", "g_inc", "g_incm", "g_exp", "g_org")
# sensible starts from the MNL baseline (step 02), NOT from the published lognormal values
m0 <- (function() { e <- new.env(); load(file.path(V2, "state_02_all.RData"), envir = e); e$lad[[SAMPLE]]$m0 })()
e0 <- m0$estimate
st1 <- c(m_sq = e0[["w_sq"]], m_nop = e0[["w_nop"]], m_d1 = e0[["w_d1"]], m_view = e0[["w_view"]], m_cost = log(e0[["b_cost"]]),
         setNames(e0[paste0("w_", c("el38", "el57", "pl3", "pl6", "wil", "d2", "d3"))], sf),
         setNames(rep(0, 15), cn), setNames(rep(0, 7), gg))
st1[c("l11", "l22", "l33", "l44")] <- c(30, 15, 10, 10)
f1 <- c(cn[11:15])                                                        # stage 1: no cost dimension
t0 <- Sys.time()
M1 <- fit_model(paste0("mvn1_", SAMPLE), dat[[SAMPLE]], prob_M, st1, fixed = f1, draws = dr_M, randcoeff = rc_M1,
                estimate_settings = list(maxIterations = 200))
t1 <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("Stage 1 (%s): four correlated Normals, fixed cost: LL %.1f, %d parameters, BIC %.1f (%.0f min; %s)\n", SAMPLE, M1$maximum, nfree(M1), bic_of(M1), t1, M1$message))
st2 <- M1$estimate; st2[c("l51", "l52", "l53", "l54")] <- 0; st2[["l55"]] <- 0.3
t0 <- Sys.time()
M2 <- fit_model(paste0("mvn2_", SAMPLE), dat[[SAMPLE]], prob_M, st2, draws = dr_M, randcoeff = rc_M2, estimate_settings = list(maxIterations = 200))
t2 <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("Stage 2 (%s): + random cost / scale dimension: LL %.1f, %d parameters, BIC %.1f (%.0f min; %s)\n", SAMPLE, M2$maximum, nfree(M2), bic_of(M2), t2, M2$message))
cat(sprintf("Reference: Anders's published model  LL %.1f, 30 parameters, BIC %.0f\n", a$LL[1], a$BIC[1]))
cat(sprintf("LR test stage 2 vs stage 1 (random cost): %.1f on %d df, p = %s\n", 2 * (M2$maximum - M1$maximum), nfree(M2) - nfree(M1), format(lrp(M2, M1), digits = 2)))

mM <- if (M2$maximum > M1$maximum + 1) M2 else M1       # keep the better model, but report which
cat("Model reported below:", if (identical(mM, M2)) "stage 2 (random cost)" else "stage 1 (cost fixed)", "\n")
e <- mM$estimate; sr <- rse(mM)
L <- matrix(0, 5, 5)
idx <- c("l11", "l21", "l31", "l41", "l51", "l22", "l32", "l42", "l52", "l33", "l43", "l53", "l44", "l54", "l55")
L[lower.tri(L, diag = TRUE)] <- e[idx]; Sig <- L %*% t(L)
nm <- c("status quo", "no protection", "< 1 km", "visible", "log cost (scale)"); dimnames(Sig) <- list(nm, nm)
cat("\nMeans of the random coefficients (NOK; compensation = sign-flipped) and standard deviations implied by the Cholesky factor\n")
mn <- c(sq = 100 * e[["m_sq"]], nop = -100 * e[["m_nop"]], d1 = -100 * e[["m_d1"]], view = -100 * e[["m_view"]])
print(data.frame(coefficient = nm[1:4], mean_NOK = round(mn), robust_se = round(100 * sr[c("m_sq", "m_nop", "m_d1", "m_view")]),
                 sd_NOK = round(100 * sqrt(diag(Sig)[1:4])), row.names = NULL), row.names = FALSE)
cat(sprintf("log cost: mean %.2f (robust s.e. %.2f), sd %.2f\n", e[["m_cost"]], sr[["m_cost"]], sqrt(Sig[5, 5])))
cat("\nImplied correlation matrix of the random coefficients\n")
cc <- cov2cor(Sig + diag(1e-12, 5)); print(round(cc, 2))
cat("\nFixed coefficients, WTA in NOK (robust t)\n")
print(data.frame(attribute = c("production 38k", "production 57k", "plants moderate", "plants many", "wilderness only", "1-3 km", "3-6 km"),
                 WTA = round(-100 * e[sf]), t = round(e[sf] / sr[sf], 1), row.names = NULL), row.names = FALSE)
cat("\nStatus quo shifters in NOK (robust t)\n")
print(data.frame(shifter = c("age (+10 y)", "women", "university degree", "income (+1 M NOK)", "income missing", "exposed", "env. organisation"),
                 effect = round(100 * e[gg]), t = round(e[gg] / sr[gg], 1), row.names = NULL), row.names = FALSE)
