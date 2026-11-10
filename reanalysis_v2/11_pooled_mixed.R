# =====================================================================
# 11  Solar vs wind in ONE mixed logit: formal likelihood-ratio tests of scale, taste heterogeneity and means
# =====================================================================
# Article: "further analysis can include testing whether the differences in preferences between the technologies are statistically significant".
# Leeds [W] 4B (pooled estimation, scale), [A] 1E (taste vs scale heterogeneity), [A] 1B (starting values; nested tests).
# The two samples are different respondents, so they are STACKED into one likelihood (id is unique).  The pooled model is parameterised as
#   solar value = theta,      wind value = theta + x          (x = wind minus solar difference)
# for all 34 parameters of the full-covariance mixed logit of step 04 (5 correlated dimensions: status quo, no protection, < 1 km, visible,
# log cost = SCALE).  With every x free this is exactly the two separate models, estimated on the SAME draws; a hypothesis is imposed by fixing
# a group of x at zero, so every test is a likelihood-ratio test between nested models (no reliance on the inflated random-coefficient t-ratios).
#   scale : x for the mean and spread of the log-scale and its 4 correlations with the tastes                 (6 df)
#   cov   : x for the Cholesky terms of the 4 taste dimensions = equal taste covariance                         (10 df)
#   attr  : x for the means of no protection, < 1 km, visible and the 7 fixed WTAs                              (10 df)
#   sq    : x for the mean status quo constant                                                                  (1 df)
#   g     : x for the 7 shifters of the status quo                                                              (7 df)
#   all   : every x = 0, i.e. wind and solar respondents share ONE preference distribution                      (34 df)
# Job name = the 2nd argument:  U (unrestricted) | scale | cov | attr | sq | g | all
job <- SAMPLE
ld1 <- function(f) { e <- new.env(); load(file.path(V2, f), envir = e); e }
ms <- ld1("state_04_solar.RData")$mM; mw <- ld1("state_04_wind.RData")$mM
pn <- names(ms$estimate); stopifnot(setequal(pn, names(mw$estimate)), length(pn) == 34)
es <- ms$estimate; ew <- mw$estimate[pn]
chol_cols <- list(c("l11", "l21", "l31", "l41", "l51"), c("l22", "l32", "l42", "l52"), c("l33", "l43", "l53"), c("l44", "l54"), "l55")
norm_chol <- function(e) { for (cl in chol_cols) if (e[[cl[1]]] < 0) e[cl] <- -e[cl]; e }   # the sign of a Cholesky column is arbitrary
es <- norm_chol(es); ew <- norm_chol(ew)
G <- list(scale = c("m_cost", "l55", "l51", "l52", "l53", "l54"),
          cov   = c("l11", "l21", "l22", "l31", "l32", "l33", "l41", "l42", "l43", "l44"),
          attr  = c("m_nop", "m_d1", "m_view", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_wil", "w_d2", "w_d3"),
          sq    = "m_sq", g = c("g_age", "g_fem", "g_univ", "g_inc", "g_incm", "g_exp", "g_org"))
stopifnot(setequal(unlist(G), pn), !anyDuplicated(unlist(G)))
# Jobs whose name ends in "2" are RESTARTS from the best solution found in the first round (job 'g', LL -4557.14).  The first-round models U, scale, attr and cov
# stopped at LL -4568 to -4571, i.e. BELOW the 'g' model that is nested in U: they had stopped early at a poor point and their LR tests are NOT valid.
restart <- grepl("2$", job); base <- sub("2$", "", job)
restr <- switch(base, U = character(0), all = pn, unlist(G[[base]]))
stopifnot(base %in% c("U", "all", names(G)))
if (!restart) {
  start_th <- ifelse(pn %in% restr, (es[pn] + ew[pn]) / 2, es[pn])        # restricted: start between the two separate solutions
  start_x  <- ifelse(pn %in% restr, 0, ew[pn] - es[pn])
  st <- c(start_th, setNames(start_x, paste0("x_", pn))); names(st)[1:34] <- pn
} else {
  gbest <- ld1("state_11_g.RData")$mP$estimate; st <- gbest[c(pn, paste0("x_", pn))]
  for (p in restr) { st[[p]] <- gbest[[p]] + gbest[[paste0("x_", p)]] / 2; st[[paste0("x_", p)]] <- 0 }   # restricted: common value = midpoint of solar and wind
}
fixed <- if (length(restr)) paste0("x_", restr) else character(0)
stopifnot(!anyDuplicated(dat$pooled$id[dat$pooled$task == 1]), all(tapply(dat$pooled$wind, dat$pooled$id, function(v) length(unique(v))) == 1))

rc_P <- function(apollo_beta, apollo_inputs) {
  mk <- function(th) list(
    sq   = th[["m_sq"]] + th[["l11"]] * xi1,
    nop  = th[["m_nop"]] + th[["l21"]] * xi1 + th[["l22"]] * xi2,
    d1   = th[["m_d1"]] + th[["l31"]] * xi1 + th[["l32"]] * xi2 + th[["l33"]] * xi3,
    view = th[["m_view"]] + th[["l41"]] * xi1 + th[["l42"]] * xi2 + th[["l43"]] * xi3 + th[["l44"]] * xi4,
    cost = exp(th[["m_cost"]] + th[["l51"]] * xi1 + th[["l52"]] * xi2 + th[["l53"]] * xi3 + th[["l54"]] * xi4 + th[["l55"]] * xi5))
  E <- environment(); gv <- function(n) get(n, envir = E)
  # Apollo attaches the parameters (apollo_beta[...] is not allowed); its internal probe calls this function with dummy parameters, hence tryCatch
  tryCatch({
    ths <- sapply(pn, gv); thw <- ths + sapply(paste0("x_", pn), gv)
    a <- mk(ths); b <- mk(thw); r <- list()
    for (n in names(a)) { r[[paste0("rs_", n)]] <- a[[n]]; r[[paste0("rw_", n)]] <- b[[n]] }
    r
  }, error = function(e) list())
}
prob_P <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs); on.exit(apollo_detach(apollo_beta, apollo_inputs))
  sel <- function(a, b) (1 - wind) * a + wind * b
  fx  <- function(n) sel(get(n), get(n) + get(paste0("x_", n)))
  bc <- sel(rs_cost, rw_cost); nop <- sel(rs_nop, rw_nop); d1 <- sel(rs_d1, rw_d1); vw <- sel(rs_view, rw_view)
  sq <- sel(rs_sq, rw_sq) + fx("g_age") * age_c + fx("g_fem") * female + fx("g_univ") * univ + fx("g_inc") * inc_c + fx("g_incm") * inc_missing + fx("g_exp") * exposed + fx("g_org") * org
  V <- list()
  V[["dev"]] <- bc * (fx("w_el38") * (el_alt1 == 38) + fx("w_el57") * (el_alt1 == 57) + fx("w_pl3") * (plants_alt1 == 3) + fx("w_pl6") * (plants_alt1 == 6) +
                      nop * (wild_alt1 == 0) + fx("w_wil") * (wild_alt1 == 1) + d1 * (dist_alt1 == 1) + fx("w_d2") * (dist_alt1 == 2) +
                      fx("w_d3") * (dist_alt1 == 3) + vw * (view_alt1 == 2) + cost_alt1 / 100)
  V[["sq"]] <- bc * sq
  P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
  P <- apollo_panelProd(P, apollo_inputs, functionality); P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality); return(P)
}
dr_P <- list(interDrawsType = "mlhs", interNDraws = 250, interNormDraws = paste0("xi", 1:5))
t0 <- Sys.time()
mP <- fit_model(paste0("pooledmix_", job), dat$pooled, prob_P, st, fixed = fixed, draws = dr_P, randcoeff = rc_P, estimate_settings = list(maxIterations = 400))
if (restart) {   # convergence check: restart once from the solution found; keep the better of the two
  ll1 <- mP$maximum
  mP_b <- fit_model(paste0("pooledmix_", job, "b"), dat$pooled, prob_P, mP$estimate, fixed = fixed, draws = dr_P, randcoeff = rc_P, estimate_settings = list(maxIterations = 400))
  cat(sprintf("First fit LL %.2f; restart from it LL %.2f\n", ll1, mP_b$maximum)); if (mP_b$maximum > mP$maximum) mP <- mP_b
}
cat(sprintf("Pooled mixed logit, job '%s': LL %.2f, %d free parameters, BIC %.1f, %d observations (%.0f min; %s)\n", job, mP$maximum, nfree(mP), bic_of(mP), mP$nObs,
            as.numeric(difftime(Sys.time(), t0, units = "mins")), mP$message))
cat(sprintf("Reference: the two separate models of step 04 sum to LL %.2f with %d parameters (different draws, so only roughly comparable)\n", ms$maximum + mw$maximum, nfree(ms) + nfree(mw)))
if (base == "U") {
  e <- mP$estimate; sr <- rse(mP); xn <- paste0("x_", pn)
  cat("\nWind minus solar difference of every parameter (x), robust t  [coefficients in 100 NOK / Cholesky units]\n")
  print(data.frame(parameter = pn, solar = round(e[pn], 3), wind_minus_solar = round(e[xn], 3), robust_t = round(e[xn] / sr[xn], 1), row.names = NULL), row.names = FALSE)
}
